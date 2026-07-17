#include "stageB_bin/DonorPool.hpp"

#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace pi0 {
namespace {

/// SplitMix64's finalising mix. Used to turn (seed, bin) into that bin's engine
/// seed.
///
/// WHY NOT JUST `mt19937_64(seed + bin)`: adjacent mt19937_64 seeds produce
/// correlated early output, so bin 0 and bin 1 would start with visibly related
/// streams -- and the pool bins are exactly the axis along which we least want a
/// correlation, since a Q^2-x_B neighbour is a physically similar bin. SplitMix64
/// is the standard, and cheap, decorrelating step.
///
/// THIS FUNCTION IS PART OF THE FROZEN ARTEFACT. Changing a constant here changes
/// every pool this project has ever built -- silently, and completely.
std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/// Unbiased uniform draw in `[0, n)`.
///
/// std::uniform_int_distribution IS DELIBERATELY NOT USED. The standard does not
/// specify its algorithm, so libc++ and libstdc++ can and do hand back different
/// values from the same engine in the same state. Using it would make the pool
/// reproducible per-toolchain -- i.e. not reproducible, while looking exactly as
/// if it were, which is worse than not trying. This is the classic rejection
/// method over the engine's full 64-bit range and is fixed for all time.
///
/// `(0 - n) % n` is 2^64 mod n in unsigned arithmetic: the size of the partial
/// block at the bottom of the range. Rejecting it leaves a whole number of blocks
/// of size n, so `u % n` is exactly uniform. Written that way rather than as
/// `2^64 % n` because 2^64 is not representable; the wraparound is intended and
/// is defined for unsigned types.
///
/// mt19937_64 is specified to produce exactly [0, 2^64), so no masking is needed.
std::uint64_t rand_below(std::mt19937_64& eng, std::uint64_t n) {
    if (n <= 1) return 0;
    const std::uint64_t threshold = (0ULL - n) % n;
    for (;;) {
        const std::uint64_t u = eng();
        if (u >= threshold) return u % n;
    }
}

}  // namespace

DonorPool::DonorPool(const PoolGridCuts& pool_grid, std::uint64_t seed,
                     std::size_t donors_per_bin)
    : m_grid(pool_grid), m_donors_per_bin(donors_per_bin) {
    m_grid.validate();

    if (m_donors_per_bin == 0) {
        throw std::invalid_argument(
            "DonorPool: donors_per_bin is 0, so every bin would hold nothing and the mixed "
            "background would be empty. It is a knob, but it is not that kind of knob.");
    }

    const std::size_t n = m_grid.n_bins();
    m_donors.resize(n);
    m_seen.assign(n, 0);

    // One engine per bin, each seeded from (seed, bin). A single shared engine
    // would make every bin's reservoir depend on how the bins interleaved in the
    // file -- reproducible, but needlessly coupled: re-running over the same
    // events in a different ORDER would rebuild a different pool even though the
    // input set is identical. Per-bin streams reduce that dependence to the order
    // WITHIN a bin, which is the irreducible part of reservoir sampling.
    m_engines.reserve(n);
    for (std::size_t b = 0; b < n; ++b) {
        m_engines.emplace_back(splitmix64(seed ^ splitmix64(b)));
    }

    // Reserving up front costs n_bins * donors_per_bin * 20 bytes -- 896 KB for
    // the shipped 224 x 200 -- and buys a build phase with no reallocation, so
    // the reservoirs never move and the pre-pass never pauses to copy.
    for (auto& bin : m_donors) bin.reserve(m_donors_per_bin);
}

void DonorPool::require_bin(int pool_bin) const {
    if (pool_bin < 0 || static_cast<std::size_t>(pool_bin) >= m_donors.size()) {
        throw std::out_of_range("DonorPool: pool bin " + std::to_string(pool_bin) +
                                " is not in [0, " + std::to_string(m_donors.size()) +
                                "). Note pool_bin() returns -1 for an event off the grid or with "
                                "no multiplicity class; -1 is not a bin -- test for it.");
    }
}

void DonorPool::decompose(int pool_bin, std::size_t& i_q2, std::size_t& i_xb,
                          std::size_t& i_mult) const {
    const std::size_t n_mult = m_grid.n_photons_classes.size();
    const std::size_t n_xb = static_cast<std::size_t>(m_grid.xb.nbins());
    const std::size_t b = static_cast<std::size_t>(pool_bin);
    i_mult = b % n_mult;
    i_xb = (b / n_mult) % n_xb;
    i_q2 = b / (n_mult * n_xb);
}

int DonorPool::pool_bin(double q2_gev2, double xb, std::size_t n_photons) const {
    const int i_q2 = m_grid.q2.find(q2_gev2);
    if (i_q2 < 0) return -1;
    const int i_xb = m_grid.xb.find(xb);
    if (i_xb < 0) return -1;
    const int i_mult = m_grid.mult_class(n_photons);
    if (i_mult < 0) return -1;  // includes n_photons == 0

    const int n_xb = m_grid.xb.nbins();
    const int n_mult = static_cast<int>(m_grid.n_photons_classes.size());
    return (i_q2 * n_xb + i_xb) * n_mult + i_mult;
}

void DonorPool::offer(int pool_bin, const std::vector<DonorPhoton>& event_photons) {
    if (m_frozen) {
        throw std::logic_error(
            "DonorPool::offer: the pool is frozen. A frozen pool is the whole design -- mixing is "
            "a pure function of (event, pool) only because the pool stopped changing before the "
            "first mixed pair existed. Offering after freeze() would resurrect the old rolling "
            "FIFO's time dependence one event at a time.");
    }
    require_bin(pool_bin);

    // THE TRIPWIRE. Not a proof of the e-gamma contract -- nothing here can prove
    // it -- but it catches the specific way it gets broken: filtering the photons
    // and not the count, or the count and not the photons. If the size handed over
    // does not fall in the class the bin encodes, the two disagree, and the pool
    // would be keyed on a multiplicity its own contents contradict.
    std::size_t i_q2 = 0, i_xb = 0, i_mult = 0;
    decompose(pool_bin, i_q2, i_xb, i_mult);
    const MultClass& cls = m_grid.n_photons_classes[i_mult];
    if (!cls.contains(event_photons.size())) {
        std::ostringstream os;
        os << "DonorPool::offer: bin " << pool_bin << " encodes multiplicity class '" << cls.label
           << "' (min " << cls.min << ", max ";
        if (cls.max.has_value()) {
            os << *cls.max;
        } else {
            os << "open";
        }
        os << "), but " << event_photons.size()
           << " photons were offered. The count used to choose the bin and the photons actually "
              "handed over disagree. The usual cause: the multiplicity was computed BEFORE the "
              "e-gamma cut and the photons were filtered AFTER it (or the reverse). Both must use "
              "the e-gamma-PASSING photons -- see DonorPool.hpp.";
        throw std::invalid_argument(os.str());
    }

    auto& reservoir = m_donors[static_cast<std::size_t>(pool_bin)];
    auto& eng = m_engines[static_cast<std::size_t>(pool_bin)];
    auto& seen = m_seen[static_cast<std::size_t>(pool_bin)];

    // Algorithm R, per photon. Note the engine is consulted ONLY once the
    // reservoir is full: while filling there is no choice to make and so no draw
    // to take. That is not an optimisation -- it fixes the stream position as a
    // function of the offer sequence alone, which is what makes two runs over the
    // same file land on the same pool.
    for (const DonorPhoton& p : event_photons) {
        ++seen;
        if (reservoir.size() < m_donors_per_bin) {
            reservoir.push_back(p);
            continue;
        }
        const std::uint64_t j = rand_below(eng, seen);
        if (j < m_donors_per_bin) {
            reservoir[static_cast<std::size_t>(j)] = p;
        }
    }
}

void DonorPool::freeze() {
    if (m_frozen) return;  // idempotent: a second freeze is confusion, not damage
    for (auto& bin : m_donors) bin.shrink_to_fit();
    m_frozen = true;
}

const std::vector<DonorPhoton>& DonorPool::donors(int pool_bin) const {
    if (!m_frozen) {
        throw std::logic_error(
            "DonorPool::donors: the pool is not frozen. Reading a reservoir mid-build hands out a "
            "list that is still moving; call freeze() first.");
    }
    require_bin(pool_bin);
    return m_donors[static_cast<std::size_t>(pool_bin)];
}

std::uint64_t DonorPool::n_offered(int pool_bin) const {
    require_bin(pool_bin);
    return m_seen[static_cast<std::size_t>(pool_bin)];
}

std::size_t DonorPool::n_filled() const {
    std::size_t n = 0;
    for (const auto& bin : m_donors) {
        if (!bin.empty()) ++n;
    }
    return n;
}

void DonorPool::report_underfilled(std::ostream& os) const {
    if (!m_frozen) {
        throw std::logic_error(
            "DonorPool::report_underfilled: the pool is not frozen. Occupancy taken mid-build "
            "describes a moving target; call freeze() first.");
    }

    const auto name = [&](int b) {
        std::size_t i_q2 = 0, i_xb = 0, i_mult = 0;
        decompose(b, i_q2, i_xb, i_mult);
        const auto& q2e = m_grid.q2.edges;
        const auto& xbe = m_grid.xb.edges;
        std::ostringstream s;
        s.precision(4);
        s << "bin " << b << "  Q2 [" << q2e[i_q2] << ", " << q2e[i_q2 + 1] << ")  xB [" << xbe[i_xb]
          << ", " << xbe[i_xb + 1] << ")  n_gamma '" << m_grid.n_photons_classes[i_mult].label
          << "'";
        return s.str();
    };

    std::vector<int> empty, thin;
    for (std::size_t b = 0; b < m_donors.size(); ++b) {
        if (m_donors[b].empty()) {
            empty.push_back(static_cast<int>(b));
        } else if (m_donors[b].size() < m_donors_per_bin) {
            thin.push_back(static_cast<int>(b));
        }
    }

    os << "DonorPool: " << m_donors.size() << " bins, " << n_filled() << " filled, depth "
       << m_donors_per_bin << " donor photons/bin.\n";

    if (empty.empty() && thin.empty()) {
        os << "DonorPool: all " << m_donors.size() << " bins reached depth.\n";
        return;
    }

    // Empty first and separately: no donors at all means NO mixed background in
    // that bin -- a different failure from a sparse one, and not one to read past
    // in a list of near-misses.
    if (!empty.empty()) {
        os << "DonorPool: WARNING -- " << empty.size()
           << " bin(s) received NO photons at all. There is no mixed background in them:\n";
        for (const int b : empty) os << "    EMPTY       " << name(b) << "\n";
    }
    if (!thin.empty()) {
        os << "DonorPool: WARNING -- " << thin.size() << " bin(s) never reached depth "
           << m_donors_per_bin << ". Their background shape is statistics-limited:\n";
        for (const int b : thin) {
            os << "    UNDERFILLED " << name(b) << "  held "
               << m_donors[static_cast<std::size_t>(b)].size() << "/" << m_donors_per_bin
               << ", offered " << m_seen[static_cast<std::size_t>(b)] << "\n";
        }
    }
}

}  // namespace pi0
