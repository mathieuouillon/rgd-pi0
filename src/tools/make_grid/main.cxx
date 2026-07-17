// ---------------------------------------------------------------------------
// make_grid -- Stage A slim files -> the two factorized equal-statistics grids.
//
// Reads one or more Stage A slim ROOT files, computes quantile (equal-
// statistics) edges for
//
//   Grid A = (Q^2, x_B)   over EVENTS,   default 8 x 7 = 56 cells
//   Grid B = (z, p_T^2)   over PI0,      default 5 x 5 = 25 cells
//
// and writes them to config/binning/grid_A_q2_xb.json and grid_B_z_pt2.json,
// each with a provenance block recording exactly what produced it.
//
// usage:
//   make_grid --input <slim.root> [--input <more.root> ...] --config <cuts.json>
//             --out-a <grid_A.json> --out-b <grid_B.json>
//             [--na 8x7] [--nb 5x5] [--max-events N]
//             [--cap-z <v>] [--cap-pt2 <v>] [--threads N]
//
// ---------------------------------------------------------------------------
// WHY THIS PROGRAM EXISTS AT ALL
// ---------------------------------------------------------------------------
// The superseded analysis built an adaptive kd-tree from an UNSEEDED,
// thread-timing-dependent reservoir sample. Two passes over the same data gave
// different edges. The edges were written to /work and never archived, so the
// binning of the published production is recoverable only from that disk, and
// only approximately from the outputs (by chaining box centres) if it is lost.
//
// So the edges are computed ONCE, here, by a deterministic function of the
// input files; written to version-controlled JSON; and hashed into every
// downstream output. Nothing about this program is allowed to depend on thread
// scheduling, iteration order or an RNG -- and nothing does: it is a sort and a
// linear interpolation. Run it twice on the same inputs and you get the same
// file, byte for byte apart from the timestamp.
//
// ---------------------------------------------------------------------------
// "EQUAL STATISTICS" MEANS EQUAL PER AXIS, NOT EQUAL PER CELL
// ---------------------------------------------------------------------------
// This is the single most misreadable thing about the output and it is stated
// here, in cuts.json, in the emitted JSON, and in the printed occupancy table.
//
// Each axis is cut at quantiles of its OWN MARGINAL distribution. Each Q^2 row
// therefore holds ~1/n_q2 of the events and each x_B column ~1/n_xb of them --
// but Q^2 and x_B are strongly correlated, so the 56 CELLS of the product are
// nowhere near equal-occupancy, and the low-x_B/high-Q^2 corner is close to
// empty. The kd-tree this replaces DID produce equal-occupancy leaves. That
// property is genuinely lost, and it is the price of having global, quotable,
// reproducible edges (note sec:binning-future).
//
// The occupancy table this program prints is not decoration. It is how a human
// sees the non-uniformity rather than taking a comment's word for it.
//
// ---------------------------------------------------------------------------
// WHY GRID B IS BUILT FROM PI0 AND NOT FROM EVENTS
// ---------------------------------------------------------------------------
// Grid B bins pi0, so its edges must be quantiles of the pi0 distribution. That
// means this program has to reproduce the pi0 reconstruction: pair the photons
// with pi0::find_gg_pairs, apply the e-gamma cut, compute z and p_T^2 with
// pi0::kin::compute_sidis. It calls the SAME functions Stage B will call, from
// the SAME config, so the grid is defined on exactly the objects it will later
// bin. A grid computed on some proxy quantity would be a grid that fits nothing.
//
// The one thing it does NOT reproduce is the mass window's effect on the
// SAMPLE: the pairs here are pi0 CANDIDATES over the full wide m_gg window
// (m_gg < 0.335, essentially every pair -- see cuts.json's
// pairing._mass_window_comment), so the (z, p_T^2) distribution these edges are
// drawn from is largely COMBINATORIAL, not signal. That is the right choice --
// the histograms Stage B fills are filled with candidates too, so the binning
// should follow the candidates -- but it is a real property of the result and
// it is recorded in the emitted JSON rather than left to be discovered.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <TFile.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <nlohmann/json.hpp>

#include "config/Cuts.hpp"
#include "core/Binning.hpp"
#include "core/Constants.hpp"
#include "core/Kinematics.hpp"
#include "core/Pairing.hpp"
#include "core/Types.hpp"
#include "util/Progress.hpp"
#include "util/Sha256.hpp"

// nlohmann/json is used here to WRITE a grid file, never to read a cut.
// Cuts::load() remains the only reader of cuts.json (see the top-level
// meson.build's note on why that dependency is private to Cuts.cpp).

#ifndef PI0_CODE_VERSION
#error "PI0_CODE_VERSION is not defined. It must be injected by meson (see src/tools/make_grid/meson.build)."
#endif

namespace {

using nlohmann::json;

// ===========================================================================
// quantiles
// ===========================================================================
//
// ONE convention, used for every axis of both grids, named in the output:
// Hyndman-Fan TYPE 7 -- the same thing numpy calls interpolation='linear' and
// the default of numpy.quantile, R's quantile(), and the scan the old analysis
// used. Given a sorted sample x[0..n-1] and p in [0, 1]:
//
//     h = (n - 1) * p
//     Q = x[floor(h)] + (h - floor(h)) * (x[floor(h) + 1] - x[floor(h)])
//
// p = 0 gives x[0] and p = 1 gives x[n-1] exactly, so the first and last edges
// ARE the sample min and max with no special-casing -- which is the required
// behaviour, and the reason type 7 is a convenient choice rather than merely a
// conventional one.
//
// WHY NAMING THE CONVENTION MATTERS ENOUGH TO SPELL IT OUT: there are nine
// standard quantile definitions and they disagree at exactly the sample sizes
// where it is hardest to notice. Two of them differ from type 7 by ~1/n in
// probability, which at n = 10^6 is invisible and at n = 800 is a bin edge
// moving by a percent. An edge array with no stated convention cannot be
// checked against an independent reimplementation, and "recoverable forever" is
// this program's entire purpose.

/// Type-7 quantile of an ALREADY-SORTED sample.
///
/// \param sorted  non-empty, ascending.
/// \param p       in [0, 1]; clamped rather than trusted.
[[nodiscard]] double quantile_type7(const std::vector<double>& sorted, double p) {
    if (sorted.empty()) throw std::runtime_error("quantile of an empty sample");
    if (sorted.size() == 1) return sorted.front();
    p = std::min(1.0, std::max(0.0, p));

    const double h = static_cast<double>(sorted.size() - 1) * p;
    // floor before the cast: h is non-negative here, but a truncating cast is
    // the wrong operation to reach for by habit and this is the one line where
    // getting it wrong would be silent.
    const double h_lo = std::floor(h);
    const std::size_t lo = static_cast<std::size_t>(h_lo);
    if (lo + 1 >= sorted.size()) return sorted.back();  // p == 1 exactly
    return sorted[lo] + (h - h_lo) * (sorted[lo + 1] - sorted[lo]);
}

/// `n_bins + 1` equal-statistics edges over `sample`, ascending.
///
/// The sample is taken BY VALUE and sorted here: the caller's vector is large
/// and it is not this function's business to reorder it.
///
/// \throws std::runtime_error if the sample is smaller than n_bins (there is
///         nothing to divide), or if the resulting edges are not STRICTLY
///         increasing.
///
/// The strictly-increasing check is not paranoia. A sample with a large atom --
/// a spike of identical values, e.g. a variable that saturates, or a p_T^2 that
/// is exactly 0 for many entries -- produces repeated quantiles and therefore a
/// ZERO-WIDTH bin. Every downstream index lookup on such a grid is ambiguous,
/// and the failure is invisible: the bin simply takes no entries and its
/// abscissa is 0/0. Better to refuse here and make the operator choose fewer
/// bins than to write a grid that silently has a hole in it.
[[nodiscard]] std::vector<double> equal_statistics_edges(std::vector<double> sample, int n_bins,
                                                         const std::string& axis_name) {
    if (n_bins < 1) throw std::runtime_error(axis_name + ": n_bins must be >= 1");
    if (sample.size() < static_cast<std::size_t>(n_bins)) {
        throw std::runtime_error(axis_name + ": only " + std::to_string(sample.size()) + " entries for " +
                                 std::to_string(n_bins) +
                                 " bins. An equal-statistics grid needs at least one entry per bin, and in "
                                 "practice needs orders of magnitude more to mean anything.");
    }
    std::sort(sample.begin(), sample.end());

    std::vector<double> edges;
    edges.reserve(static_cast<std::size_t>(n_bins) + 1);
    for (int i = 0; i <= n_bins; ++i) {
        edges.push_back(quantile_type7(sample, static_cast<double>(i) / static_cast<double>(n_bins)));
    }

    for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
        if (!(edges[i] < edges[i + 1])) {
            std::ostringstream os;
            os << axis_name << ": edges " << i << " and " << i + 1 << " are both " << std::setprecision(17)
               << edges[i]
               << ", i.e. a zero-width bin. More than " << 100.0 / n_bins
               << "% of the sample sits at one value, so it cannot be split into " << n_bins
               << " equal-statistics bins. Use fewer bins on this axis, or bin a different variable.";
            throw std::runtime_error(os.str());
        }
    }
    return edges;
}

/// Index of `x` in `edges`, or nullopt if it is outside.
///
/// Convention, applied identically everywhere in this project: every bin is
/// [lo, hi) EXCEPT the last, which is [lo, hi] so that the sample maximum --
/// which IS the top edge, by construction of the type-7 quantile at p = 1 --
/// lands in the top bin rather than in the overflow. Downstream code must use
/// this same rule; it is stated in the emitted JSON for that reason.
[[nodiscard]] std::optional<std::size_t> find_bin(const std::vector<double>& edges, double x) {
    if (edges.size() < 2) return std::nullopt;
    if (x < edges.front() || x > edges.back()) return std::nullopt;
    if (x == edges.back()) return edges.size() - 2;  // the closed top end
    const auto it = std::upper_bound(edges.begin(), edges.end(), x);
    if (it == edges.begin()) return std::nullopt;
    return static_cast<std::size_t>(it - edges.begin()) - 1;
}

// ===========================================================================
// CLI
// ===========================================================================

struct Args {
    std::vector<std::string> inputs;
    std::string config;
    std::string out_a;
    std::string out_b;
    std::optional<int> na_q2, na_xb;  ///< --na WxH, overriding cuts.json's /binning
    std::optional<int> nb_z, nb_pt2;  ///< --nb WxH, ditto
    std::optional<long long> max_events;
    std::optional<double> cap_z;
    std::optional<double> cap_pt2;
    /// Worker threads for the sample-collection read. DEFAULT 1, which is
    /// byte-for-byte the old sequential program. THE EDGES DO NOT DEPEND ON THIS:
    /// each axis's edges are quantiles of the pooled, SORTED sample, so the order
    /// (and therefore the thread count) in which the samples are collected cannot
    /// change them. Threads only decide how many input files are read at once.
    /// Ignored -- forced back to 1 -- when --max-events is given; see parse_args.
    unsigned int threads{1};
};

void print_usage(const char* argv0) {
    std::cerr
        << "usage: " << argv0 << " --input <slim.root> [--input <more.root> ...]\n"
        << "                      --config <cuts.json>\n"
        << "                      --out-a <grid_A.json> --out-b <grid_B.json>\n"
        << "                      [--na 8x7] [--nb 5x5] [--max-events N]\n"
        << "                      [--cap-z <v>] [--cap-pt2 <v>] [--threads N]\n"
        << "\n"
        << "  --input       a Stage A slim file. Repeat for more; the grid is computed over the\n"
        << "                union of all of them. Order does not affect the edges (they are\n"
        << "                quantiles of a sorted pooled sample), only the provenance listing.\n"
        << "  --config      the cut configuration. The grid SHAPE, the pairing cuts and the z\n"
        << "                window all come from here. Its sha256 is stamped into both outputs\n"
        << "                and checked against the one each input recorded at skim time.\n"
        << "  --out-a       Grid A = (Q^2, x_B), per event.\n"
        << "  --out-b       Grid B = (z, p_T^2), per pi0.\n"
        << "  --na WxH      override binning.grid_a from the config, e.g. 8x7.\n"
        << "  --nb WxH      override binning.grid_b from the config, e.g. 5x5.\n"
        << "  --max-events  stop after N events TOTAL across all inputs (default: all). The\n"
        << "                sample is then a PREFIX of the input list, not a random subset, so\n"
        << "                the edges it gives are not the edges of the full dataset. Stamped\n"
        << "                into the output so a truncated grid cannot pass for a real one.\n"
        << "  --cap-z       clamp Grid B's top z edge to this and DISCARD pi0 above it.\n"
        << "  --cap-pt2     the same for p_T^2, in GeV^2.\n"
        << "  --threads     worker threads for reading the inputs (default 1). THE EDGES DO NOT\n"
        << "                DEPEND ON THIS: each axis's edges are quantiles of the pooled, sorted\n"
        << "                sample, so the order the samples are collected in -- and hence the\n"
        << "                thread count -- cannot change them. The inputs are read one whole file\n"
        << "                per worker, so N > 1 helps only with more than one --input (production\n"
        << "                is thousands of slims per target). IGNORED under --max-events, which is\n"
        << "                read strictly single-threaded so the truncation stays a deterministic\n"
        << "                prefix of the input list. Progress goes to stderr; stdout stays the\n"
        << "                occupancy table and edges.\n"
        << "\n"
        << "THE CAPS ARE OFF BY DEFAULT AND SHOULD USUALLY STAY OFF. See the long comment in\n"
        << "this program's source, or config/cuts.json's /binning/_cap_comment: with the\n"
        << "count-weighted abscissae Stage B accumulates, a wide outer bin is honestly\n"
        << "POSITIONED even though it is coarse, so a cap trades real pi0 for a cosmetic gain.\n"
        << "Both caps are recorded in the output along with how many pi0 they threw away.\n";
}

/// Parse "8x7" into a pair. Accepts 'x' or 'X'.
[[nodiscard]] std::pair<int, int> parse_shape(const std::string& flag, const std::string& v) {
    const std::size_t x = v.find_first_of("xX");
    if (x == std::string::npos || x == 0 || x + 1 == v.size()) {
        throw std::runtime_error(flag + " wants <bins>x<bins>, e.g. 8x7, got \"" + v + "\"");
    }
    const auto one = [&](const std::string& s) {
        try {
            std::size_t pos = 0;
            const int n = std::stoi(s, &pos);
            if (pos != s.size()) throw std::invalid_argument("trailing characters");
            if (n < 1) throw std::invalid_argument("not positive");
            return n;
        } catch (const std::exception&) {
            throw std::runtime_error(flag + " wants two positive integers, e.g. 8x7, got \"" + v + "\"");
        }
    };
    return {one(v.substr(0, x)), one(v.substr(x + 1))};
}

[[nodiscard]] double parse_double(const std::string& flag, const std::string& v) {
    try {
        std::size_t pos = 0;
        const double d = std::stod(v, &pos);
        if (pos != v.size()) throw std::invalid_argument("trailing characters");
        if (!std::isfinite(d)) throw std::invalid_argument("not finite");
        return d;
    } catch (const std::exception&) {
        throw std::runtime_error(flag + " wants a finite number, got \"" + v + "\"");
    }
}

[[nodiscard]] Args parse_args(int argc, char* argv[]) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(flag + " requires a value");
            return argv[++i];
        };
        if (flag == "--input") a.inputs.push_back(value());
        else if (flag == "--config") a.config = value();
        else if (flag == "--out-a") a.out_a = value();
        else if (flag == "--out-b") a.out_b = value();
        else if (flag == "--na") {
            const auto [w, h] = parse_shape(flag, value());
            a.na_q2 = w;
            a.na_xb = h;
        } else if (flag == "--nb") {
            const auto [w, h] = parse_shape(flag, value());
            a.nb_z = w;
            a.nb_pt2 = h;
        } else if (flag == "--cap-z") {
            a.cap_z = parse_double(flag, value());
        } else if (flag == "--cap-pt2") {
            a.cap_pt2 = parse_double(flag, value());
        } else if (flag == "--max-events") {
            const std::string v = value();
            long long n = 0;
            try {
                // The '-' check is not decoration: without it "-1" parses as a
                // negative that then compares as "never reached", i.e. a typo
                // that silently means the opposite of what it says.
                if (v.find('-') != std::string::npos) throw std::invalid_argument("negative");
                std::size_t pos = 0;
                n = std::stoll(v, &pos);
                if (pos != v.size()) throw std::invalid_argument("trailing characters");
            } catch (const std::exception&) {
                throw std::runtime_error("--max-events wants a positive integer, got \"" + v + "\"");
            }
            if (n == 0) {
                throw std::runtime_error(
                    "--max-events 0 would read nothing. Omit the flag to use every event.");
            }
            a.max_events = n;
        } else if (flag == "--threads") {
            const std::string v = value();
            unsigned long long n = 0;
            try {
                // Same '-' guard as --max-events: stoull turns "-1" into a huge
                // positive that would spawn an absurd pool, so a typo must be
                // refused rather than silently doing the opposite.
                if (v.find('-') != std::string::npos) throw std::invalid_argument("negative");
                std::size_t pos = 0;
                n = std::stoull(v, &pos);
                if (pos != v.size()) throw std::invalid_argument("trailing characters");
            } catch (const std::exception&) {
                throw std::runtime_error("--threads wants a positive integer, got \"" + v + "\"");
            }
            if (n == 0) {
                throw std::runtime_error("--threads 0 is meaningless; the minimum is 1 (the default).");
            }
            if (n > std::numeric_limits<unsigned int>::max()) {
                throw std::runtime_error("--threads " + v + " is absurd; the maximum is " +
                                         std::to_string(std::numeric_limits<unsigned int>::max()) + ".");
            }
            a.threads = static_cast<unsigned int>(n);
        } else {
            throw std::runtime_error("unknown argument \"" + flag + "\"");
        }
    }

    if (a.inputs.empty()) throw std::runtime_error("at least one --input is required");
    if (a.config.empty()) throw std::runtime_error("--config is required");
    if (a.out_a.empty()) throw std::runtime_error("--out-a is required");
    if (a.out_b.empty()) throw std::runtime_error("--out-b is required");
    if (a.out_a == a.out_b) throw std::runtime_error("--out-a and --out-b are the same file");
    return a;
}

// ===========================================================================
// reading a slim file
// ===========================================================================

/// One input's contribution, plus what its provenance said about it.
struct InputSummary {
    std::string path;
    long long n_events_read{};
    long long n_pi0{};
    /// Selected Stage A provenance, verbatim. Empty string = the key was absent.
    std::string stagea_config_sha, stagea_run, stagea_target, stagea_fallback, stagea_code_version;
};

/// One TNamed's title from the file's provenance directory, or "" if absent.
[[nodiscard]] std::string read_provenance(TFile& f, const char* key) {
    auto* dir = f.GetDirectory("provenance");
    if (dir == nullptr) return {};
    auto* obj = dir->Get<TNamed>(key);
    return obj ? std::string(obj->GetTitle()) : std::string{};
}

/// Accumulated samples. Deliberately four flat vectors rather than a vector of
/// structs: they are sorted independently (each axis is its own marginal), and
/// at production scale these are the program's whole memory footprint.
struct Samples {
    std::vector<double> q2, xb;    ///< per event
    std::vector<double> z, pt2;    ///< per pi0, index-aligned with each other
};

/// Read one slim file, appending to `s`.
///
/// \param budget  remaining --max-events allowance; decremented. Negative means
///                unlimited.
///
/// APPENDS TO ITS OWN `s`. Under multi-threading each worker owns a distinct
/// Samples, filled from one whole file, and the caller concatenates them in
/// input-index order afterwards. Because the edges are quantiles of the pooled
/// SORTED sample, that concatenation gives an identical pooled sample -- and thus
/// identical edges -- for any thread count, with no compensated summation and no
/// fixed partitioning of the kind Stage B needs. This is the whole reason
/// make_grid can be parallelised more simply than Stage B.
void read_slim(const std::string& path, const pi0::Cuts& cuts, Samples& s, InputSummary& sum, long long& budget) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (f == nullptr || f->IsZombie()) throw std::runtime_error("cannot open input: " + path);

    sum.path = path;
    sum.stagea_config_sha = read_provenance(*f, "config.sha256");
    sum.stagea_run = read_provenance(*f, "run");
    sum.stagea_target = read_provenance(*f, "target");
    sum.stagea_fallback = read_provenance(*f, "gbt.fallback_used");
    sum.stagea_code_version = read_provenance(*f, "stageA_skim.code_version");

    auto* tree = f->Get<TTree>("events");
    if (tree == nullptr) {
        throw std::runtime_error("no TTree named \"events\" in " + path +
                                 " -- is this a Stage A slim file? (make_grid reads the slim schema, "
                                 "not raw HIPO.)");
    }

    TTreeReader reader(tree);
    // The Stage A schema. A missing branch is caught by TTreeReader's status
    // check below rather than by reading garbage.
    TTreeReaderValue<double> r_q2(reader, "q2");
    TTreeReaderValue<double> r_xb(reader, "xb");
    TTreeReaderValue<double> r_nu(reader, "nu");
    TTreeReaderValue<double> r_ex(reader, "ex");
    TTreeReaderValue<double> r_ey(reader, "ey");
    TTreeReaderValue<double> r_ez(reader, "ez");
    TTreeReaderValue<double> r_ee(reader, "ee");
    TTreeReaderArray<double> r_gpx(reader, "gpx");
    TTreeReaderArray<double> r_gpy(reader, "gpy");
    TTreeReaderArray<double> r_gpz(reader, "gpz");
    TTreeReaderArray<double> r_gang(reader, "g_e_gamma_deg");

    std::vector<pi0::Photon> photons;

    // Whether WE ended the loop or the reader did. The status check below is
    // only meaningful in the second case -- see it for why.
    bool stopped_on_budget = false;

    while (reader.Next()) {
        if (budget == 0) {
            stopped_on_budget = true;
            break;
        }
        if (budget > 0) --budget;
        ++sum.n_events_read;

        // ---- Grid A: per EVENT --------------------------------------------
        // q2 and xb are taken from the file, NOT recomputed. Stage A already
        // computed them and the contract says so; recomputing here would be a
        // second implementation of the same formula that could drift from the
        // one that actually applied the DIS cuts.
        s.q2.push_back(*r_q2);
        s.xb.push_back(*r_xb);

        // ---- Grid B: per PI0 -----------------------------------------------
        // The e-gamma cut is applied HERE, on the precomputed angle, because
        // Stage A deliberately did not apply it (it only stored the angle) and
        // find_gg_pairs deliberately does not either -- it is a single-photon
        // property, so it belongs upstream of pairing. See Pairing.hpp.
        photons.clear();
        const std::size_t n_g = std::min({r_gpx.GetSize(), r_gpy.GetSize(), r_gpz.GetSize(), r_gang.GetSize()});
        for (std::size_t i = 0; i < n_g; ++i) {
            if (!(r_gang[i] > cuts.pairing.e_gamma_min_angle_deg)) continue;  // strict, per cuts.json
            photons.push_back(pi0::Photon{r_gpx[i], r_gpy[i], r_gpz[i], r_gang[i]});
        }
        if (photons.size() < 2) continue;

        // The DIS kinematics come from the file. compute_sidis reads only
        // dis.nu, but the whole struct is filled from the branches rather than
        // part-filled, so nothing here depends on which field it happens to use.
        pi0::DisKin dis{};
        dis.q2 = *r_q2;
        dis.nu = *r_nu;
        dis.xb = *r_xb;

        for (const pi0::GGPair& pair : pi0::find_gg_pairs(photons, cuts.pairing)) {
            const pi0::SidisKin sk = pi0::kin::compute_sidis(pair.px, pair.py, pair.pz, pair.e, dis, *r_ex, *r_ey,
                                                             *r_ez, *r_ee, cuts.beam.energy_gev);
            // 0 < z < 1, strict at both ends (cuts.json pairing._z_comment).
            // NOT optional here: the mass window is so wide that the candidate
            // sample is largely combinatorial, and a combinatorial pair can
            // carry more energy than nu -- z > 1, unphysical. Letting one into
            // the sample would drag the top z edge past the kinematic limit,
            // which is precisely the defect this whole exercise is about.
            if (!(sk.z > cuts.sidis.z_min && sk.z < cuts.sidis.z_max)) continue;
            s.z.push_back(sk.z);
            s.pt2.push_back(sk.pt2);
            ++sum.n_pi0;
        }
    }

    // Only ask the reader why the loop ended if the READER ended it. The status
    // describes the read that ended the loop; break on the budget and the last
    // Next() succeeded, so it is kEntryValid -- a healthy read this check used to
    // report as "failed reading ... Check the branch list matches the Stage A
    // schema", blaming the schema for a file that is fine. It made --max-events
    // unusable on any file holding more events than the budget, which is every
    // file the flag is for. The old form guessed that a clean break "leaves
    // kEntryNotFound behind on some ROOT versions"; it does not on ROOT 6.40, and
    // the guess was never needed -- we know who ended the loop.
    //
    // kEntryBeyondEnd is the ONLY clean end: per TTreeReader.h it is "last entry
    // loop has reached its end", while kEntryNotFound is "the tree entry number
    // does not exist" -- a real failure alongside kEntryChainSetupError and
    // kEntryDictionaryError. Excusing it here, now that `stopped_on_budget` says
    // what it was standing in for, would mask the error it names.
    if (!stopped_on_budget && reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd) {
        // A missing branch shows up as kEntryDictionaryError / kEntryChainSetupError.
        throw std::runtime_error("failed reading \"events\" from " + path +
                                 " (TTreeReader status " + std::to_string(static_cast<int>(reader.GetEntryStatus())) +
                                 "). Check the branch list matches the Stage A schema.");
    }
}

// ===========================================================================
// output
// ===========================================================================

[[nodiscard]] std::string utc_now() {
    const std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    std::array<char, 32> buf{};
    std::strftime(buf.data(), buf.size(), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return std::string(buf.data());
}

/// A cap as it ended up being applied, for the JSON and for the printout.
struct CapRecord {
    bool applied{};
    double value{};
    long long n_discarded{};
    std::string axis;
};

json cap_to_json(const CapRecord& c) {
    json j;
    j["applied"] = c.applied;
    if (c.applied) {
        j["value"] = c.value;
        j["n_pi0_discarded"] = c.n_discarded;
        j["_comment"] =
            "A cap WAS applied to this axis. The top edge below is the cap value exactly, not the "
            "sample maximum, and the " + std::to_string(c.n_discarded) +
            " pi0 above it were DISCARDED -- they are in no bin of this grid and Stage B will not "
            "count them. Note the discard is per-PI0 and therefore affects BOTH axes of Grid B: a "
            "pi0 dropped by --cap-z is gone from the p_T^2 marginal too, so every edge on this "
            "grid, not just this axis's top one, was computed from the surviving sample.";
    } else {
        j["_comment"] =
            "No cap. The top edge is the sample maximum, so the outermost bin runs to the "
            "kinematic limit of the data and is far wider than its neighbours. THIS IS THE "
            "DEFAULT AND IS DELIBERATE: Stage B accumulates count-weighted sums per (4D bin, "
            "m_gg bin), so the outermost bin is reported at the mean of the pi0 actually in it, "
            "not at the midpoint of its box. A wide bin that is honestly positioned is coarse; a "
            "wide bin reported at its geometric centre is wrong, and that is what the superseded "
            "analysis did (note sec:binning-caveat: 6.9% of its bins implied nu > E_beam).";
    }
    return j;
}

/// One axis as it will be written.
struct AxisOut {
    std::string name;            ///< MUST be what Binning::load expects: q2/xb, z/pt2.
    std::vector<double> edges;
    std::string comment;
    json cap;
};

/// Write one grid file.
///
/// THE SCHEMA IS NOT THIS PROGRAM'S TO CHOOSE. src/core/Binning.hpp's
/// Binning::load() is the sole reader of these files, and it is strict on
/// purpose: an `axes` ARRAY of exactly two objects, each {name, edges}, in an
/// order that is load-bearing (axes[0] is the SLOW axis of the row-major flat
/// index), with the names checked against the grid -- "q2","xb" for A and
/// "z","pt2" for B -- plus a `provenance` object carrying at least a string
/// `source` and an integer `n_events`. Anything else it refuses.
///
/// So this function writes that shape and main() then LOADS BOTH FILES BACK
/// through Binning::load() before claiming success. A generator whose output
/// its own loader rejects is worse than useless: it produces a plausible file
/// that fails at the start of the next stage, hours later, with the scan gone.
///
/// Extra keys are free: nlohmann ignores what the reader does not ask for, and
/// _comment / quantile_convention / cap are exactly the things a human opening
/// this file in two years needs and the loader does not.
void write_grid(const std::string& path, const std::string& grid_label, const std::string& comment,
                const std::vector<AxisOut>& axes, const json& provenance) {
    json j;
    j["_comment"] = comment;
    j["_schema_version"] = "1.0";
    j["_grid"] = grid_label;
    j["_generated"] = utc_now();
    j["_bin_convention"] =
        "Every bin is half-open [lo, hi) EXCEPT that a value exactly on the top edge lands in the "
        "LAST bin rather than in the overflow. This is Grid1D::find()'s documented rule and it is "
        "not negotiable here: the top edge is the kinematic limit (and, with no cap, is exactly "
        "the sample maximum), so returning 'outside' for it would silently discard the events that "
        "reach it -- which is what the superseded find_1d_bin() did to every p_T^2 sitting on the "
        "top edge. Any decoder written against this file must implement the same rule.";
    j["_quantile_convention"] =
        "Hyndman-Fan type 7, i.e. numpy.quantile's default interpolation='linear', and R's "
        "quantile() default. For a sorted sample x[0..n-1] and probability p: h = (n-1)*p, "
        "Q = x[floor(h)] + (h - floor(h)) * (x[floor(h)+1] - x[floor(h)]). Consequently p=0 gives "
        "the sample minimum and p=1 the sample maximum exactly, so the first and last edges ARE "
        "the sample min and max unless a cap moved the last one. Stated explicitly, and in the "
        "file rather than only in the source, because the nine standard quantile definitions "
        "differ from one another by ~1/n in probability -- invisible at n = 1e6, a percent-level "
        "edge shift at n = 1e3 -- so an edge array whose convention is unstated cannot be checked "
        "against an independent reimplementation. 'Recoverable forever' means recoverable by "
        "someone who does not have this source tree.";
    j["_equal_statistics_is_marginal"] =
        "EQUAL STATISTICS IS PER AXIS, NOT PER CELL. Each axis below is cut at quantiles of its "
        "OWN marginal, so each bin of each axis holds ~1/nbins of the sample. The 2D CELLS of the "
        "product do NOT hold equal counts: the two axes are correlated, so the cell occupancy is "
        "strongly non-uniform and some corner cells are nearly empty. The superseded kd-tree DID "
        "equalise its leaves; that is the adaptivity traded away for global, quotable, "
        "reproducible edges (note sec:binning-future). make_grid prints the full occupancy table "
        "when it writes this file.";

    json axes_json = json::array();
    for (const AxisOut& a : axes) {
        json ax;
        // `name` and `edges` are the two keys Binning::load actually reads.
        // ORDER IS LOAD-BEARING: axes[0] is the slow axis. The caller passes
        // them in index order and the loader re-checks the names anyway, so a
        // swap here is caught rather than silently mislabelling every result.
        ax["name"] = a.name;
        ax["_comment"] = a.comment;
        ax["edges"] = a.edges;
        ax["n_bins"] = a.edges.size() - 1;
        ax["cap"] = a.cap;
        axes_json.push_back(std::move(ax));
    }
    j["axes"] = axes_json;
    j["provenance"] = provenance;

    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write: " + path);
    // Full double precision. The superseded analysis wrote its display copy of
    // the binning at the iostream default of SIX significant digits, which made
    // the only human-readable record of it display-only -- and the 17-digit copy
    // was never archived. A grid file that cannot round-trip a double does not
    // define the binning it claims to. nlohmann serialises doubles with a
    // shortest-round-trip algorithm, so this is belt and braces on the stream
    // rather than the mechanism, but it costs nothing and the failure it guards
    // against is unrecoverable.
    out << std::setprecision(17) << j.dump(2) << '\n';
    if (!out) throw std::runtime_error("failed writing: " + path);
}

// ===========================================================================
// printing
// ===========================================================================

void print_edges(const std::string& axis, const std::vector<double>& e) {
    std::cout << "  " << std::left << std::setw(8) << axis << " [";
    for (std::size_t i = 0; i < e.size(); ++i) {
        std::cout << (i ? ", " : "") << std::fixed << std::setprecision(4) << e[i];
    }
    std::cout << "]   (" << e.size() - 1 << " bins)\n" << std::defaultfloat;
}

/// The 2D occupancy table, plus the marginals.
///
/// The marginals are the point: they are what "equal statistics" actually
/// promises, and printing them next to the manifestly non-uniform cells is the
/// cheapest way to stop somebody reading the cell counts as a bug.
void print_occupancy(const std::string& title, const std::string& xname, const std::vector<double>& xe,
                     const std::string& yname, const std::vector<double>& ye, const std::vector<double>& xs,
                     const std::vector<double>& ys) {
    const std::size_t nx = xe.size() - 1, ny = ye.size() - 1;
    std::vector<long long> cells(nx * ny, 0), rowsum(nx, 0), colsum(ny, 0);
    long long total = 0, outside = 0;

    for (std::size_t k = 0; k < xs.size(); ++k) {
        const auto ix = find_bin(xe, xs[k]);
        const auto iy = find_bin(ye, ys[k]);
        if (!ix || !iy) {
            ++outside;
            continue;
        }
        ++cells[*ix * ny + *iy];
        ++rowsum[*ix];
        ++colsum[*iy];
        ++total;
    }

    std::cout << "\n--- " << title << " occupancy ---\n";
    std::cout << "  rows = " << xname << " bins, columns = " << yname << " bins\n\n";

    std::cout << "        ";
    for (std::size_t j = 0; j < ny; ++j) std::cout << std::right << std::setw(8) << (yname + std::to_string(j));
    std::cout << std::setw(10) << "ROW SUM" << '\n';

    for (std::size_t i = 0; i < nx; ++i) {
        std::cout << "  " << std::left << std::setw(6) << (xname + std::to_string(i));
        for (std::size_t j = 0; j < ny; ++j) std::cout << std::right << std::setw(8) << cells[i * ny + j];
        std::cout << std::setw(10) << rowsum[i] << '\n';
    }

    std::cout << "  " << std::left << std::setw(6) << "COLSUM";
    for (std::size_t j = 0; j < ny; ++j) std::cout << std::right << std::setw(8) << colsum[j];
    std::cout << std::setw(10) << total << '\n';

    // The honest summary of what was and was not achieved.
    long long lo = std::numeric_limits<long long>::max(), hi = 0;
    std::size_t empty = 0;
    for (long long c : cells) {
        lo = std::min(lo, c);
        hi = std::max(hi, c);
        if (c == 0) ++empty;
    }
    const double ideal = nx * ny > 0 ? static_cast<double>(total) / static_cast<double>(nx * ny) : 0.0;
    std::cout << "\n  cells: " << nx * ny << "   empty: " << empty << "   min " << lo << "   max " << hi
              << "   flat-grid ideal " << std::fixed << std::setprecision(1) << ideal << std::defaultfloat << '\n';
    if (outside > 0) {
        std::cout << "  outside the grid: " << outside
                  << "   (only possible when a cap discarded them; otherwise the edges span the sample)\n";
    }
    std::cout << "  The ROW SUM and COLSUM columns are what 'equal statistics' means here: each is\n"
              << "  flat to within rounding. THE CELLS ARE NOT FLAT AND ARE NOT MEANT TO BE -- the\n"
              << "  two axes are correlated, and a factorized grid cannot follow that. The kd-tree\n"
              << "  this replaces did equalise the cells; that is the adaptivity traded away for\n"
              << "  reproducible global edges (note sec:binning-future).\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    Args args;
    try {
        args = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n\n";
        print_usage(argv[0]);
        return 1;
    }

    try {
        const pi0::Cuts cuts = pi0::Cuts::load(args.config);
        const std::string config_sha = pi0::util::sha256_file(args.config);

        // Shape: config first, command line overrides. Which one won is recorded
        // in the output -- a grid whose shape came from an ad-hoc flag and a
        // grid whose shape came from the version-controlled config are not the
        // same artefact, and only one of them is reproducible from the repo.
        const int n_q2 = args.na_q2.value_or(cuts.binning.n_q2);
        const int n_xb = args.na_xb.value_or(cuts.binning.n_xb);
        const int n_z = args.nb_z.value_or(cuts.binning.n_z);
        const int n_pt2 = args.nb_pt2.value_or(cuts.binning.n_pt2);
        const bool na_overridden = args.na_q2.has_value();
        const bool nb_overridden = args.nb_z.has_value();

        std::cout << "make_grid\n"
                  << "  config     : " << args.config << "  (sha256 " << config_sha.substr(0, 16) << "...)\n"
                  << "  inputs     : " << args.inputs.size() << '\n'
                  << "  Grid A     : " << n_q2 << " x " << n_xb << " (Q2, xB)"
                  << (na_overridden ? "   [--na override]" : "   [from config]") << '\n'
                  << "  Grid B     : " << n_z << " x " << n_pt2 << " (z, pT2)"
                  << (nb_overridden ? "   [--nb override]" : "   [from config]") << '\n'
                  << "  max events : "
                  << (args.max_events.has_value()
                          ? std::to_string(*args.max_events) + "   *** TRUNCATED (--max-events) ***"
                          : std::string("all"))
                  << '\n'
                  << "  cap z      : " << (args.cap_z ? std::to_string(*args.cap_z) : std::string("none")) << '\n'
                  << "  cap pT2    : " << (args.cap_pt2 ? std::to_string(*args.cap_pt2) : std::string("none")) << '\n';

        // ---- read ----------------------------------------------------------
        // Multi-threading is OPT-IN (--threads, default 1) and READS ONE WHOLE FILE
        // PER WORKER. It is FORCED OFF under --max-events: that flag means "N events
        // TOTAL across all inputs", and the only cheap way to keep that a
        // DETERMINISTIC prefix of the input list is a single shared budget consumed
        // in input order -- the sequential loop below. Racing workers on one atomic
        // budget would still stop at N events, but WHICH N (and so the sample, and
        // so the edges) would depend on who won the races. Truncated grids are for
        // smoke tests, where determinism matters more than the seconds MT would
        // save, so this trade is free. See the --threads usage text.
        const bool use_threads = args.threads > 1 && !args.max_events.has_value();

        // Each input's samples land in its OWN slot; concatenated in input-index
        // order afterwards, so the pooled sample -- and every edge -- is identical
        // for any thread count (edges are quantiles of the sorted pool; order does
        // not reach them). `read_ok[i]` records which inputs were actually read: the
        // sequential path may stop early on the budget, leaving a tail unread.
        std::vector<Samples> per_input(args.inputs.size());
        std::vector<InputSummary> per_summary(args.inputs.size());
        std::vector<char> read_ok(args.inputs.size(), 0);
        bool truncated = false;  // set iff the budget stopped us before the last input

        // Progress on STDERR only -- stdout is the result (occupancy table + edges),
        // which may be parsed. The unit is INPUT FILES, not events: make_grid has no
        // cheap up-front event count (unlike Stage B, whose pass 0 hands pass 1 a
        // total for free), so an event-based bar would need a serial pre-pass that
        // opens every input just to sum GetEntries() -- a serial O(n_files) cost
        // that grows with the input and fights the very parallelism this adds, worst
        // on the thousands-of-slims production runs where progress matters most. One
        // tick per completed file needs no pre-pass, is thread-safe, and is smooth
        // exactly when there are many files. add() is called by the reader loops.
        pi0::Progress progress("make_grid", static_cast<std::int64_t>(args.inputs.size()), std::cerr);

        if (!use_threads) {
            // Sequential -- byte-for-byte the old program, and the ONLY path that
            // honours --max-events. One shared budget, consumed in input order.
            long long budget = args.max_events.value_or(-1);
            for (std::size_t i = 0; i < args.inputs.size(); ++i) {
                read_slim(args.inputs[i], cuts, per_input[i], per_summary[i], budget);
                read_ok[i] = 1;
                progress.add();
                if (budget == 0) {
                    truncated = i + 1 < args.inputs.size();
                    break;
                }
            }
        } else {
            // Parallel: file-level, work-stealing over an atomic index so uneven
            // file sizes load-balance (production is thousands of slims per target).
            // No --max-events on this path, so every worker reads its whole file
            // with an unlimited (negative) budget -- nothing is truncated here.
            //
            // ROOT's global state (TFile::Open, TClass, streamers) is made
            // thread-safe only now that more than one file is read at once; the
            // --threads 1 path never touches it and stays the plain sequential flow.
            ROOT::EnableThreadSafety();

            const unsigned int lanes =
                static_cast<unsigned int>(std::min<std::size_t>(args.threads, args.inputs.size()));
            std::atomic<std::size_t> next_input{0};
            std::vector<std::exception_ptr> errs(args.inputs.size());

            const auto worker = [&]() {
                std::size_t i;
                while ((i = next_input.fetch_add(1, std::memory_order_relaxed)) < args.inputs.size()) {
                    try {
                        long long unlimited = -1;  // a per-worker budget; never spent
                        read_slim(args.inputs[i], cuts, per_input[i], per_summary[i], unlimited);
                        read_ok[i] = 1;
                        progress.add();  // thread-safe: atomic count, try_lock redraw
                    } catch (...) {
                        // A worker's throw must not cross the thread boundary; stash
                        // it and rethrow on the main thread after the join.
                        errs[i] = std::current_exception();
                    }
                }
            };

            std::vector<std::thread> pool;
            pool.reserve(lanes > 0 ? lanes - 1 : 0);
            for (unsigned int k = 1; k < lanes; ++k) pool.emplace_back(worker);
            worker();  // this thread is a worker too
            for (std::thread& t : pool) t.join();

            for (const std::exception_ptr& e : errs) {
                if (e) std::rethrow_exception(e);
            }
        }
        progress.finish();

        // Concatenate the per-input samples and gather the per-input summaries in
        // INPUT-INDEX ORDER, so stdout and the pooled sample are deterministic and,
        // for --threads 1, byte-identical to before. Only inputs actually read
        // contribute -- the budget may have left a tail unread.
        Samples s;
        std::vector<InputSummary> summaries;
        for (std::size_t i = 0; i < args.inputs.size(); ++i) {
            if (!read_ok[i]) continue;
            const Samples& p = per_input[i];
            s.q2.insert(s.q2.end(), p.q2.begin(), p.q2.end());
            s.xb.insert(s.xb.end(), p.xb.begin(), p.xb.end());
            s.z.insert(s.z.end(), p.z.begin(), p.z.end());
            s.pt2.insert(s.pt2.end(), p.pt2.begin(), p.pt2.end());
            std::cout << "  read " << std::setw(9) << per_summary[i].n_events_read << " events, " << std::setw(9)
                      << per_summary[i].n_pi0 << " pi0   from " << args.inputs[i] << '\n';
            summaries.push_back(std::move(per_summary[i]));
        }
        if (truncated) {
            std::cout << "  (--max-events reached; remaining inputs not read)\n";
        }

        // ---- the config cross-check ----------------------------------------
        // The grid is defined by the pairing cuts and the z window in
        // args.config, but the pi0 it is computed from were selected by whatever
        // config the SKIM ran under. If those differ, the edges describe a
        // selection that never produced this data. Cheap to check because Stage
        // A stamps its config's sha256 into every file.
        {
            std::vector<std::string> mismatched;
            for (const auto& sm : summaries) {
                if (!sm.stagea_config_sha.empty() && sm.stagea_config_sha != config_sha) mismatched.push_back(sm.path);
            }
            if (!mismatched.empty()) {
                std::cout << "\n  WARNING: " << mismatched.size() << " of " << summaries.size()
                          << " inputs were skimmed under a DIFFERENT config than --config:\n";
                for (const auto& m : mismatched) std::cout << "           " << m << '\n';
                std::cout << "           The photons in those files were selected by cuts that are not the ones\n"
                          << "           this grid is being defined with. The mismatch is recorded in both output\n"
                          << "           files. It is a warning and not a refusal because a cuts.json edit that\n"
                          << "           touches only, say, an extraction knob changes the hash without changing\n"
                          << "           the skim -- but if the /photon or /pairing blocks moved, re-skim.\n";
            }
        }

        long long n_events = 0, n_pi0 = 0;
        for (const auto& sm : summaries) {
            n_events += sm.n_events_read;
            n_pi0 += sm.n_pi0;
        }
        std::cout << "\n  total: " << n_events << " events, " << n_pi0 << " pi0\n";

        if (n_events == 0) throw std::runtime_error("no events were read from any input; there is nothing to bin");
        if (n_pi0 == 0) {
            throw std::runtime_error(
                "no pi0 survived pairing + the e-gamma cut + the z window, so Grid B cannot be computed. "
                "Grid A alone is not a useful output, so nothing is written.");
        }

        // ---- caps ----------------------------------------------------------
        // Applied BEFORE the quantiles, by deleting the overflow from the
        // sample. The discard is per-PI0 and therefore couples the axes: a pi0
        // dropped for its z is gone from the p_T^2 marginal too. Done in one
        // pass over the index-aligned vectors so that coupling is structural
        // rather than something two independent filters might disagree about.
        CapRecord cap_z{false, 0.0, 0, "z"}, cap_pt2{false, 0.0, 0, "pt2"};
        if (args.cap_z || args.cap_pt2) {
            const double zc = args.cap_z.value_or(std::numeric_limits<double>::infinity());
            const double pc = args.cap_pt2.value_or(std::numeric_limits<double>::infinity());
            std::vector<double> kz, kp;
            kz.reserve(s.z.size());
            kp.reserve(s.pt2.size());
            for (std::size_t i = 0; i < s.z.size(); ++i) {
                const bool over_z = s.z[i] > zc;
                const bool over_p = s.pt2[i] > pc;
                if (over_z) ++cap_z.n_discarded;
                if (over_p) ++cap_pt2.n_discarded;
                if (over_z || over_p) continue;
                kz.push_back(s.z[i]);
                kp.push_back(s.pt2[i]);
            }
            s.z.swap(kz);
            s.pt2.swap(kp);
            if (args.cap_z) { cap_z.applied = true; cap_z.value = zc; }
            if (args.cap_pt2) { cap_pt2.applied = true; cap_pt2.value = pc; }

            std::cout << "  cap: kept " << s.z.size() << " of " << n_pi0 << " pi0";
            if (cap_z.applied) std::cout << "   (z > " << cap_z.value << " dropped " << cap_z.n_discarded << ")";
            if (cap_pt2.applied) std::cout << "   (pT2 > " << cap_pt2.value << " dropped " << cap_pt2.n_discarded << ")";
            std::cout << "\n  NOTE the discards above OVERLAP: a pi0 over both caps is counted in both, so\n"
                      << "       they do not sum to the number removed. The kept count is the authority.\n";
            if (s.z.empty()) throw std::runtime_error("the caps discarded every pi0; there is nothing to bin");
        }

        // ---- edges ---------------------------------------------------------
        std::vector<double> e_q2 = equal_statistics_edges(s.q2, n_q2, "Grid A / Q2");
        std::vector<double> e_xb = equal_statistics_edges(s.xb, n_xb, "Grid A / xB");
        std::vector<double> e_z = equal_statistics_edges(s.z, n_z, "Grid B / z");
        std::vector<double> e_pt2 = equal_statistics_edges(s.pt2, n_pt2, "Grid B / pT2");

        // The cap clamps the TOP EDGE to the cap value exactly, rather than
        // leaving it at the highest surviving pi0. Two reasons: the file then
        // says what was asked for instead of an accident of the sample, and the
        // bin is closed at a number a human chose, so two runs on different data
        // with the same cap have the same top edge.
        if (cap_z.applied) e_z.back() = cap_z.value;
        if (cap_pt2.applied) e_pt2.back() = cap_pt2.value;

        std::cout << "\n--- Grid A edges (Q2, xB) -- per EVENT ---\n";
        print_edges("Q2", e_q2);
        print_edges("xB", e_xb);
        std::cout << "\n--- Grid B edges (z, pT2) -- per PI0 ---\n";
        print_edges("z", e_z);
        print_edges("pT2", e_pt2);

        print_occupancy("Grid A (Q2, xB)", "Q2_", e_q2, "xB_", e_xb, s.q2, s.xb);
        print_occupancy("Grid B (z, pT2)", "z_", e_z, "pT2_", e_pt2, s.z, s.pt2);

        // ---- provenance ----------------------------------------------------
        json inputs = json::array();
        for (const auto& sm : summaries) {
            json i;
            i["path"] = sm.path;
            i["n_events_used"] = sm.n_events_read;
            i["n_pi0_used"] = sm.n_pi0;
            i["stageA.config_sha256"] = sm.stagea_config_sha;
            i["stageA.config_sha256_matches_this_run"] =
                sm.stagea_config_sha.empty() ? std::string("unknown (no provenance in the file)")
                : sm.stagea_config_sha == config_sha ? std::string("yes")
                                                     : std::string("NO -- this file was skimmed under other cuts");
            i["stageA.run"] = sm.stagea_run;
            i["stageA.target"] = sm.stagea_target;
            i["stageA.gbt_fallback_used"] = sm.stagea_fallback;
            i["stageA.code_version"] = sm.stagea_code_version;
            inputs.push_back(std::move(i));
        }

        const auto make_prov = [&](const std::string& which) {
            json p;
            // `source` and `n_events` are REQUIRED by Binning::load() and are
            // type-checked there (string, integer). They are first here because
            // they are the two fields the loader will refuse the file without;
            // everything after them is for humans and for re-running this scan.
            //
            // `source` is the scan's identity, in one line, per the contract the
            // placeholder files state: it replaces the literal "placeholder", so
            // a grep for that word across config/binning is exactly the check
            // for "has this grid ever been fitted to data".
            std::ostringstream src;
            src << "make_grid over " << summaries.size() << " Stage A slim file(s), " << n_events << " events, "
                << n_pi0 << " pi0";
            if (!summaries.empty()) {
                src << "; targets=";
                std::vector<std::string> t;
                for (const auto& sm : summaries) {
                    if (!sm.stagea_target.empty() && std::find(t.begin(), t.end(), sm.stagea_target) == t.end()) {
                        t.push_back(sm.stagea_target);
                    }
                }
                for (std::size_t i = 0; i < t.size(); ++i) src << (i ? "," : "") << t[i];
                if (t.empty()) src << "unknown";
                src << "; runs=";
                std::vector<std::string> rr;
                for (const auto& sm : summaries) {
                    if (!sm.stagea_run.empty() && std::find(rr.begin(), rr.end(), sm.stagea_run) == rr.end()) {
                        rr.push_back(sm.stagea_run);
                    }
                }
                for (std::size_t i = 0; i < rr.size(); ++i) src << (i ? "," : "") << rr[i];
                if (rr.empty()) src << "unknown";
            }
            src << "; config.sha256=" << config_sha.substr(0, 16);
            if (args.max_events.has_value()) src << "; TRUNCATED --max-events=" << *args.max_events;
            p["source"] = src.str();
            // The loader wants an INTEGER and wants it to mean "how many events
            // the equal-statistics edges were computed from". For Grid B that is
            // still the event count -- n_pi0 is carried separately below rather
            // than smuggled into this field, because the two grids must answer
            // the same question with the same units or the field is a trap.
            p["n_events"] = n_events;
            p["n_pi0"] = n_pi0;
            p["_comment"] =
                "Everything needed to recompute these edges. The point of this block is that the "
                "binning of a production must be recoverable from version control FOREVER. The "
                "superseded analysis's kd-tree edges were written to /work and never archived: if "
                "that disk dies, the binning of its published results is gone and the numbers "
                "cannot be re-derived, only approximated by chaining box centres out of the "
                "outputs (note sec:binning-reconstruction).";
            p["make_grid.code_version"] = PI0_CODE_VERSION;
            p["make_grid.created_utc"] = utc_now();
            p["grid"] = which;
            p["inputs"] = inputs;
            p["n_input_files"] = static_cast<long long>(summaries.size());
            p["n_events_used"] = n_events;
            p["n_pi0_used"] = n_pi0;
            p["config.path"] = args.config;
            p["config.sha256"] = config_sha;
            p["beam.energy_gev"] = cuts.beam.energy_gev;
            p["shape.source"] = (which == "A" ? na_overridden : nb_overridden)
                                    ? std::string("--na/--nb on the command line, OVERRIDING config binning.grid_" +
                                                  std::string(which == "A" ? "a" : "b") +
                                                  ". This grid's shape is NOT reproducible from the repository "
                                                  "alone -- you need this command line too.")
                                    : std::string("config binning.grid_" + std::string(which == "A" ? "a" : "b"));
            p["max_events"] = args.max_events.has_value()
                                  ? std::to_string(*args.max_events) +
                                        " -- TRUNCATED: the sample is a PREFIX of the input list, not a random "
                                        "subset of it, so these edges are not the edges of the full dataset."
                                  : std::string("all (no --max-events; every event in every input was used)");
            return p;
        };

        json prov_a = make_prov("A");
        prov_a["_selection_comment"] =
            "Grid A is computed from the per-EVENT (q2, xb) branches of the slim files, taken as "
            "written by Stage A rather than recomputed -- a second implementation of the DIS "
            "formulae here could drift from the one that actually applied the DIS cuts. Every "
            "event in the input contributes, including events whose photons form no pi0.";

        json prov_b = make_prov("B");
        prov_b["_selection_comment"] =
            "Grid B is computed from PI0 CANDIDATES, not events and not truth pi0. Each event's "
            "photons are filtered on g_e_gamma_deg > " + std::to_string(cuts.pairing.e_gamma_min_angle_deg) +
            " deg (Stage A stores this angle but deliberately does not cut on it, and "
            "find_gg_pairs deliberately does not either -- it is a single-photon property), paired "
            "by pi0::find_gg_pairs under this config's /pairing block, turned into (z, pT2) by "
            "pi0::kin::compute_sidis against the event's stored nu, and required to satisfy " +
            std::to_string(cuts.sidis.z_min) + " < z < " + std::to_string(cuts.sidis.z_max) +
            ". These are the same functions and the same config Stage B uses, so the grid is "
            "defined on exactly the objects it will later bin. WHAT THIS SAMPLE IS NOT: the m_gg "
            "window is deliberately wide (m_gg < 0.335, essentially every pair), so these edges "
            "are quantiles of a largely COMBINATORIAL sample rather than of the pi0 signal. That "
            "is the correct choice -- the histograms Stage B fills are filled with candidates too "
            "-- but it means the edges follow the background's (z, pT2) distribution wherever the "
            "background dominates.";

        json factorized_comment =
            "FACTORIZED, AND DELIBERATELY SO. Grid B's edges are computed ONCE over ALL pi0 -- they "
            "are NOT recomputed per Grid A cell. The 4D bin is the plain product A x B (56 x 25 = "
            "1400), not a nested tree. This is a deliberate simplification of the superseded "
            "adaptive kd-tree, which split Q2 -> xB -> z -> pT2 in that order and so had a "
            "different z edge set per (Q2, xB) cell and a different pT2 set per (Q2, xB, z) cell "
            "-- 1450 distinct pT2 centres against 290 z centres. What the change BUYS: global "
            "quotable edges in every dimension, and determinism (the kd-tree was built from an "
            "unseeded, thread-timing-dependent reservoir sample, so two passes over the same data "
            "gave different trees). What it COSTS, and this is a real loss, not a formality: pT2 "
            "and z are strongly correlated, and factorized edges cannot follow that correlation "
            "the way the nested tree did.";

        // Axis names are NOT free text: Binning::load() checks axes[0].name and
        // axes[1].name against "q2"/"xb" for A and "z"/"pt2" for B, and refuses
        // a swapped pair -- which would otherwise pass every geometry check and
        // simply label every result with the wrong variable.
        write_grid(args.out_a, "A",
                   "Grid A: equal-statistics (quantile) edges in (Q^2, x_B), computed per EVENT from the "
                   "q2/xb branches of the Stage A slim files. " +
                       factorized_comment.get<std::string>(),
                   {AxisOut{"q2",
                            e_q2,
                            "Q^2 in GeV^2, the SLOW axis of the flat A index (cell = i_q2 * n_xb + i_xb). The "
                            "bottom edge is the sample minimum, which sits at dis.q2_min = 1.0 because that is "
                            "the cut every event passed. The TOP edge is the sample maximum, i.e. the kinematic "
                            "limit reached by this dataset -- so the top bin is far wider than its neighbours "
                            "and MUST NOT be reported at its midpoint (note sec:binning-caveat).",
                            cap_to_json(CapRecord{false, 0, 0, "q2"})},
                    AxisOut{"xb",
                            e_xb,
                            "x_B, dimensionless, the FAST axis of the flat A index. There is no x_B cut in the "
                            "analysis (cuts.json /dis), so unlike the pool grid's hand-chosen [0.1, 0.7] these "
                            "edges span the data itself and truncate nothing at either end.",
                            cap_to_json(CapRecord{false, 0, 0, "xb"})}},
                   prov_a);

        write_grid(args.out_b, "B",
                   "Grid B: equal-statistics (quantile) edges in (z, p_T^2), computed per PI0 CANDIDATE. " +
                       factorized_comment.get<std::string>(),
                   {AxisOut{"z",
                            e_z,
                            "z = E_pi0 / nu, dimensionless, the SLOW axis of the flat B index (cell = i_z * "
                            "n_pt2 + i_pt2). Restricted to " +
                                std::to_string(cuts.sidis.z_min) + " < z < " + std::to_string(cuts.sidis.z_max) +
                                " before the quantiles were taken: the m_gg window is wide enough that the "
                                "candidate sample is largely combinatorial, and a combinatorial pair can carry "
                                "more energy than nu (z > 1, unphysical). One of those in the sample would drag "
                                "the top edge past the kinematic limit, which is the exact defect this grid "
                                "exists to avoid.",
                            cap_to_json(cap_z)},
                    AxisOut{"pt2",
                            e_pt2,
                            "p_T^2 in GeV^2 relative to the virtual photon (Trento frame, pi0::kin::"
                            "compute_sidis), the FAST axis of the flat B index. The top bin is the one the "
                            "note's Cronin-rise caveat is about: in the old production the top p_T^2 box spanned "
                            "[0.18, 3.64] and was reported at 1.92 while its true <p_T^2> was ~0.3.",
                            cap_to_json(cap_pt2)}},
                   prov_b);

        std::cout << "\nwrote " << args.out_a << "\n      " << args.out_b << '\n';

        // ---- load the output back -------------------------------------------
        // The generator proves its own product is loadable, using the SAME
        // loader Stage B will use, before it reports success. A grid file that
        // Binning::load() rejects is a plausible-looking file that fails at the
        // start of the next stage -- hours later, with the scan gone and the
        // farm slot spent. The schema Binning::load() enforces (an axes ARRAY,
        // ordered, names checked, provenance.source a string and
        // provenance.n_events an integer) is narrow enough that writing it
        // blind and hoping is not a reasonable thing to do.
        //
        // This also prints the provenance_hash -- the 16 hex digits that get
        // stamped into every downstream output and are how a result is matched
        // back to the binning that produced it. That is the hole this whole
        // rewrite exists to close: the superseded production's outputs carry
        // nothing that could identify their kd-tree.
        const pi0::Binning check = pi0::Binning::load(args.out_a, args.out_b);
        std::cout << "\n--- readback (pi0::Binning::load, the same loader Stage B uses) ---\n"
                  << "  Grid A cells  : " << check.A.ncells() << "\n"
                  << "  Grid B cells  : " << check.B.ncells() << "\n"
                  << "  4D bins       : " << check.n4d() << "\n"
                  << "  3D bins (pT)  : " << check.n3d() << "\n"
                  << "  provenance_hash: " << check.provenance_hash() << "\n"
                  << "  ^ stamp this into every output made with this binning.\n";
        if (args.max_events.has_value()) {
            std::cout << "\n*** --max-events was given: these edges come from a PREFIX of the inputs and are\n"
                      << "    not the edges of the full dataset. Both files say so in their provenance. ***\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 2;
    }
}
