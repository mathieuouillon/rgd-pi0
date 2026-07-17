#include "core/Pairing.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

#include "core/Constants.hpp"
#include "core/Kinematics.hpp"

namespace pi0 {
namespace {

/// The massless four-vector sum of two photons.
struct PairVec {
    double px{}, py{}, pz{}, e{};
    double mgg{};       ///< GeV, clamped at 0
    double p_mag{};     ///< |p_1 + p_2|, GeV
    double theta_deg{}; ///< opening angle between the two photons, DEGREES
};

[[nodiscard]] PairVec combine(const Photon& a, const Photon& b) {
    PairVec v{};
    v.px = a.px + b.px;
    v.py = a.py + b.py;
    v.pz = a.pz + b.pz;
    v.e = a.e() + b.e();

    v.p_mag = std::sqrt(v.px * v.px + v.py * v.py + v.pz * v.pz);

    // m^2 = E^2 - |p|^2. Algebraically this equals 2 E1 E2 (1 - cos theta_12)
    // and so is >= 0, but E and |p| are built by different sqrt chains, so for a
    // collinear pair they agree only to within rounding and the difference lands
    // on either side of 0. This clamp is NOT dead defensive code: over 5e6
    // exactly-collinear pairs at 0.05-6 GeV, m2 came out NEGATIVE 50.0% of the
    // time, down to -9.8e-14. Every one of those would be a NaN out of sqrt, and
    // a NaN mgg silently fails every downstream comparison (including the cuts
    // below, so the pair would vanish with no trace of why).
    const double m2 = v.e * v.e - v.p_mag * v.p_mag;
    v.mgg = (m2 > 0.0) ? std::sqrt(m2) : 0.0;

    v.theta_deg = kin::angle_between_deg(a.px, a.py, a.pz, b.px, b.py, b.pz);
    return v;
}

}  // namespace

double min_opening_angle_deg(double pair_momentum_gev, const PairingCuts& cuts) {
    return cuts.open_a_deg * std::exp(-cuts.open_b_inv_gev * pair_momentum_gev) + cuts.open_offset_deg;
}

std::optional<GGPair> admissible_pair(const Photon& a, const Photon& b, const PairingCuts& cuts) {
    const PairVec v = combine(a, b);

    // --- admissibility ------------------------------------------------------
    // THE ONLY COPY OF THESE THREE TESTS IN THE PROJECT. find_gg_pairs() calls
    // this rather than repeating them, and stageB_bin's mixed pass calls it too
    // -- which is what makes "the mixed pairs use the same cuts as the same-event
    // pairs" a fact about the code rather than a claim in a comment.
    if (!(v.mgg > cuts.min_mgg_gev)) return std::nullopt;
    if (!(v.theta_deg > min_opening_angle_deg(v.p_mag, cuts))) return std::nullopt;
    if (!(std::abs(v.mgg - kPi0MassGeV) < cuts.mass_window_gev)) return std::nullopt;

    GGPair p{};
    // i and j index nothing here: this function has no list. find_gg_pairs()
    // overwrites them; a mixed-pair caller must not read them.
    p.i = 0;
    p.j = 0;
    p.mgg = v.mgg;
    p.px = v.px;
    p.py = v.py;
    p.pz = v.pz;
    p.e = v.e;
    return p;
}

std::vector<GGPair> find_gg_pairs(const std::vector<Photon>& photons, const PairingCuts& cuts) {
    std::vector<GGPair> pairs;

    const std::size_t n = photons.size();
    if (n < 2) return pairs;

    std::vector<bool> used(n, false);
    pairs.reserve(n / 2);

    // Greedy exclusive pairing: each sweep picks the single best admissible
    // pair over the whole remaining photon list, consumes it, and re-sweeps.
    // At most floor(n/2) sweeps can succeed, so the loop always terminates.
    for (std::size_t sweep = 0; sweep < n / 2; ++sweep) {
        bool found = false;
        double best_dm = std::numeric_limits<double>::max();
        std::size_t best_i = 0;
        std::size_t best_j = 0;
        GGPair best_vec{};

        for (std::size_t i = 0; i < n; ++i) {
            if (used[i]) continue;
            for (std::size_t j = i + 1; j < n; ++j) {
                if (used[j]) continue;

                // The admissibility tests live in admissible_pair(), which is
                // this function's own predicate given a name so that the mixed
                // pass can call it too. See Pairing.hpp.
                const std::optional<GGPair> v = admissible_pair(photons[i], photons[j], cuts);
                if (!v.has_value()) continue;

                const double dm = std::abs(v->mgg - kPi0MassGeV);

                // Strict `<` keeps the first-found pair on a tie, and the loops
                // run in lexicographic (i, j) order, so ties resolve to the
                // lowest indices deterministically.
                if (dm < best_dm) {
                    best_dm = dm;
                    best_i = i;
                    best_j = j;
                    best_vec = *v;
                    found = true;
                }
            }
        }

        if (!found) break;

        // admissible_pair() left i and j at 0 because it had no list to index.
        // This function does, so it fills them in.
        best_vec.i = best_i;
        best_vec.j = best_j;
        pairs.push_back(best_vec);

        used[best_i] = true;
        used[best_j] = true;
    }

    return pairs;
}

}  // namespace pi0
