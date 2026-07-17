// ---------------------------------------------------------------------------
// stageA_skim -- HIPO -> slim ROOT TTree.
//
// One input HIPO file per invocation. Reads the raw banks through RHipoDS,
// applies the electron selection (including the target's vertex window), the
// DIS cuts and the GBT photon selection, and writes one entry per surviving
// event into a TTree named "events".
//
// usage:
//   stageA_skim --input <file.hipo> --output <slim.root> --config <cuts.json>
//               --target <LD2|CxC|Cu|Sn> [--run <N>] [--max-events <N>]
//               [--qa-ntuple <qa.root>]
//
// ---------------------------------------------------------------------------
// WHY THERE IS A SECOND, OPT-IN OUTPUT (--qa-ntuple)
// ---------------------------------------------------------------------------
// The slim schema holds only what SURVIVED selection. It carries no vz, no
// chi2pid, no sampling fraction, no PCAL lv/lw, no DC edge and no GBT score, so
// the figures that show the selection working -- sampling fraction vs p with its
// band, vz against the target window, the GBT score against its threshold --
// cannot be drawn from it: every value those plots need is one the slim dropped,
// and the rows they need are the rejected ones the slim never wrote.
//
// --qa-ntuple writes those PRE-CUT values to a separate file, one row per
// CANDIDATE, recorded as the cuts are applied. It is DIAGNOSTIC ONLY and both
// files' provenance says so: it is per-candidate rather than per-event, its
// photon rows exist only for events that already passed the electron and DIS
// cuts, and nothing in it is a yield.
//
// It RECORDS the selection; it does not re-run it. Every number written below
// comes out of the same call the skim decides on -- pass_electron()'s verdict,
// the score the GBT threshold was compared against, the lv/lw handed to the
// fiducial cut. src/stageB_bin/main.cxx documents what the alternative costs:
// its duplicated e-gamma filter means "the pool is drawn from a different photon
// population than the spectra it models". A QA plot drawn from a second
// implementation illustrates a selection that never ran.
//
// ---------------------------------------------------------------------------
// WHY >= 1 PHOTON AND NOT >= 2
// ---------------------------------------------------------------------------
// A pi0 needs two photons, so it is tempting to require two here. Do not. The
// mixed-event background (stage C) draws its donor photons from single-photon
// events as much as from multi-photon ones; requiring two would silently shrink
// the donor pool and change the background shape, not just the file size.
//
// ---------------------------------------------------------------------------
// WHY THE MODEL IS CHOSEN ONCE, BEFORE THE EVENT LOOP
// ---------------------------------------------------------------------------
// select_model() THROWS for RG-D runs unless photon.allow_rga_fallback is set
// (see config/cuts.json, which explains at length why). That throw is the
// intended behaviour and this program's job is to report it and exit non-zero,
// not to work around it. Doing the lookup before the loop means the refusal
// costs a file open instead of a full pass over the data.
//
// ---------------------------------------------------------------------------
// EVERY COLUMN IS AN RVec
// ---------------------------------------------------------------------------
// RHipoDS is constructed with n_inspect = 0, which is deliberate: a schema scan
// mis-types a multi-row bank as a scalar whenever the scanned events happen to
// hold one row, and drops banks that are absent from the scan window. With
// n_inspect = 0 the schema comes from the bank dictionary instead, so it is
// stable -- at the price that EVERY column is an RVec, including the single-row
// banks. REC_Event_helicity is RVec<short>, not short. RUN_config_run is
// RVec<int>, not int. Read element [0] and guard for empty.
//
// ---------------------------------------------------------------------------
// OPT-IN MULTI-THREADING (--threads N, DEFAULT 1) -- AND WHY ORDER IS SACRED
// ---------------------------------------------------------------------------
// Stage A writes a slim TTree "events". Stage B builds its donor pool by
// RESERVOIR SAMPLING, which is a pure function of the SEQUENCE of entries in that
// slim. So the slim -- and the three QA trees (qa_electron, qa_photon,
// qa_electron_stages) -- must be BYTE-IDENTICAL regardless of --threads: a
// thread-nondeterministic emission order would make the pool, and hence the
// published mixed-event background, nondeterministic. Same philosophy as the
// Stage B pass-1 MT that just landed (see src/stageB_bin/main.cxx): parallelise
// the WORK, emit the OUTPUT in canonical (entry) order.
//
// The part we parallelise is the per-event work -- above all the GBT photon
// scoring (CatBoost). The TTree::Fill is NOT parallelised; it is serialised on
// purpose, because that is what fixes the emission order.
//
//   * --threads 1 is EXACTLY the pre-MT path: one RHipoDS + RDataFrame + Foreach
//     that fills the trees directly, in entry order, with O(1) memory. The farm
//     runs one thread per job across many files, so it is on this path and is
//     entirely unaffected by everything below. Its output is byte-for-byte what
//     the code produced before MT existed (Stage A has no floating-point
//     reduction, unlike Stage B, so it is EXACTLY equal, not equal-to-ulp).
//
//   * --threads N>1 reads the events ONCE, sequentially, into an in-memory buffer
//     (copying each event's raw columns), then SCORES them in parallel: the buffer
//     index range [0, total) is split into a FIXED number of contiguous partitions
//     decided by the EVENT COUNT, never by the thread count (pi0::partition_count /
//     pi0::partition_range from core/Summation.hpp -- the same primitives Stage B
//     uses). Each partition runs process_event() over its slice on a worker thread,
//     appending output rows to per-partition buffers. Because a partition is a
//     contiguous index range, its rows are already in entry order; the partitions
//     are then filled into the trees SINGLE-THREADED, in PARTITION-INDEX ORDER, so
//     the concatenation is the exact entry-order sequence one thread would emit.
//     TTree::Fill is never called off the main thread.
//
// WHY THE READ IS SHARED, NOT PARTITIONED. The obvious design -- one RHipoDS per
// partition, each seeking to its own entry range -- does NOT work for HIPO, and
// this was measured, not assumed. On run 22083 the ENTIRE wall time is a FIXED,
// one-time cost inside RHipoDS: constructing the datasource and asking it for its
// entry ranges decompresses the whole file's record index (GetEntryRanges loads
// every record), and it is SERIAL and PER-INSTANCE. Reading 50 events and reading
// 4000 events from this 9.1 GB file both take ~56 s; the per-event work is in the
// noise. So N independent readers pay that whole-file scan N times over and run N
// times SLOWER, and RDataFrame::Range makes it worse still (it never seeks -- the
// loop calls SetEntry for every entry in [0, end) and Range only gates the action).
// One shared reader pays the scan once; only the cheap, embarrassingly-parallel
// per-event scoring is threaded. On THIS smoke file that means --threads shows no
// wall-clock win (the scan dominates and cannot be parallelised from here); on a
// GBT-bound production run -- the RG-A fallback scoring many photons across
// millions of events -- the parallel scoring is what pays off. Either way the
// OUTPUT is identical, which is the property that actually matters.
//
// Determinism here does NOT depend on the partition COUNT (there is no sum whose
// rounding could shift): ANY contiguous tiling filled in index order gives the
// same bytes. The fixed, data-decided count is kept only to mirror Stage B.
//
// MEMORY, STATED HONESTLY. --threads N>1 buffers the raw events (O(events)) plus
// the output rows, so it is meant for BOUNDED (--max-events) smoke runs. The farm
// runs --threads 1, which streams and buffers nothing (O(1)). Do not point
// --threads N>1 at a full 64M-event file.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TTree.h>

// external/hipo-cpp is a pristine v4.4.1 checkout under a standing "DO NOT
// MODIFY" rule, and its headers trip -Wunused-parameter (bank.h's stubbed-out
// `composite` constructors, reader.h's copy constructor). Those warnings are
// unactionable -- nobody here is allowed to fix them -- and at warning_level=3
// five permanent lines of noise are exactly what stops anyone reading the
// warnings that DO mean something. Suppressed at the include and nowhere else,
// so this file's own code is still held to the strict level. Same reasoning as
// the warning_level=0 override on pi0_photonid in the top-level meson.build.
//
// `#pragma GCC diagnostic` is understood by both gcc and clang; spelling it
// `clang` would break the other one.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "RHipoDS.hxx"
#pragma GCC diagnostic pop

#include "config/Cuts.hpp"
#include "core/Constants.hpp"
#include "core/Kinematics.hpp"
#include "core/Summation.hpp"
#include "photonid/Features.hpp"
#include "photonid/PhotonGBT.hpp"
#include "photonid/RunRangeModelMap.hpp"
#include "selection/ElectronSelection.hpp"
#include "selection/PhotonSelection.hpp"
#include "selection/SamplingFraction.hpp"
#include "util/Progress.hpp"
#include "util/Provenance.hpp"
#include "util/Sha256.hpp"
#include "vertex/VzCorrector.hpp"

// The vz correction is defined by data/Vz/vz_corrector_params.txt. meson injects
// its absolute path; there is deliberately no fallback, because a skim that
// silently ran without the correction would produce a plausible-looking file
// with the wrong target assignment -- which is worse than not running at all.
#ifndef PI0_VZ_PARAMS_FILE
#error "PI0_VZ_PARAMS_FILE is not defined. It must be injected by meson (see src/stageA_skim/meson.build)."
#endif

// Set from meson.project_version(). See the provenance note in main().
#ifndef PI0_CODE_VERSION
#error "PI0_CODE_VERSION is not defined. It must be injected by meson (see src/stageA_skim/meson.build)."
#endif

namespace {

using RVecI = ROOT::VecOps::RVec<int>;
using RVecS = ROOT::VecOps::RVec<short>;
using RVecF = ROOT::VecOps::RVec<float>;

// ===========================================================================
// SHA-256
// ===========================================================================
// Lives in util/Sha256.hpp, not here. It used to be a local class in this
// anonymous namespace; make_grid needs the SAME digest (it stamps the config
// sha256 into the grid JSON and compares it against the one this program wrote
// into each slim file it reads), and a second copy of a hash function is a
// second copy of a constant. See the header for the full argument.

// ===========================================================================
// CLI
// ===========================================================================

/// Default approximate events per partition when --threads > 1. Coarser than
/// Stage B's 200000 because a Stage A event is far heavier (five CatBoost models),
/// so fewer, larger partitions still keep every core busy while bounding the
/// per-wave buffer memory. Only used for the opt-in multi-threaded path; a small
/// --max-events smoke run wants a smaller --partition-target to expose parallelism.
constexpr std::size_t kPartitionTargetDefault = 20000;

struct Args {
    std::string input;
    std::string output;
    std::string config;
    std::string target;
    std::optional<int> run;  ///< overrides RUN_config_run[0] for the model lookup
    /// Stop after this many tag-0 events. nullopt = the whole file.
    ///
    /// `unsigned int` rather than a wider type because that is what
    /// RDataFrame::Range takes; parse_args range-checks into it rather than
    /// letting a huge --max-events wrap around into a small one.
    std::optional<unsigned int> max_events;
    /// Second, DIAGNOSTIC-ONLY output: per-candidate pre-cut values. nullopt =
    /// not requested, and nothing about the run changes.
    std::optional<std::string> qa_ntuple;
    /// Worker threads for the GBT scoring / read. DEFAULT 1. THE OUTPUT DOES NOT
    /// DEPEND ON THIS: the entry range is split into a fixed number of contiguous
    /// partitions decided by the event count, and each partition's rows are filled
    /// into the trees in partition-index (i.e. entry) order, so the slim and the
    /// QA trees are byte-identical for any thread count. Threads only decide how
    /// many partitions score in parallel. --threads 1 is the pre-MT single-pass
    /// path, exactly. See the file header for the memory/I-O caveats of N>1.
    unsigned int threads{1};
    /// Approximate events per partition, a TUNING knob for parallel granularity,
    /// NOT a correctness knob: pi0::partition_count is a pure function of (event
    /// count, this target) and never of --threads, so the output is byte-identical
    /// across thread counts for any fixed value here. Only consulted when
    /// --threads > 1. Lower it to expose more parallelism on a small (--max-events)
    /// input; the default is coarse. See core/Summation.hpp.
    std::size_t partition_target{kPartitionTargetDefault};
};

void print_usage(const char* argv0) {
    std::cerr
        << "usage: " << argv0 << " --input <file.hipo> --output <slim.root>\n"
        << "                        --config <cuts.json> --target <LD2|CxC|Cu|Sn>\n"
        << "                        [--run <N>] [--max-events <N>] [--qa-ntuple <qa.root>]\n"
        << "                        [--threads <N>] [--partition-target <N>]\n"
        << "\n"
        << "  --input       ONE HIPO file. Entry numbering restarts per file; do not chain inputs.\n"
        << "  --output      slim ROOT file to create (overwritten if it exists).\n"
        << "  --config      the cut configuration. Every threshold comes from here.\n"
        << "  --target      RG-D target. Selects the vertex window, and for Cu/Sn IS the\n"
        << "                target assignment -- the two foils are distinguished only by v_z.\n"
        << "  --run         override the run number used to pick the GBT model. Normally the\n"
        << "                run is read from RUN::config; pass this only when you know the\n"
        << "                file's run number is wrong or absent.\n"
        << "  --max-events  stop after the first N tag-0 events (default: all of them).\n"
        << "                For smoke tests. The output is a PREFIX of the file, not a\n"
        << "                sample of it, so it is not a physics-representative subset --\n"
        << "                the truncation is stamped into the provenance and the cutflow\n"
        << "                so a partial run cannot be mistaken for a full one.\n"
        << "  --qa-ntuple   also write a DIAGNOSTIC-ONLY file of per-candidate PRE-CUT values\n"
        << "                (trees \"qa_electron\" and \"qa_photon\"), for plotting the selection\n"
        << "                and what each cut removes. Omitting it changes nothing about the\n"
        << "                slim. Rows are CANDIDATES, not events: take no yield from it.\n"
        << "  --threads     worker threads for the GBT scoring (default 1). THE OUTPUT DOES\n"
        << "                NOT DEPEND ON THIS: the entry range is split into a fixed number of\n"
        << "                contiguous partitions decided by the event count, and each\n"
        << "                partition's rows are filled in entry order, so the slim and the QA\n"
        << "                trees are byte-identical for any thread count. --threads 1 is the\n"
        << "                pre-MT single-pass path. N>1 buffers rows and re-unpacks each\n"
        << "                partition's lead-in, so it is for BOUNDED (--max-events) runs; the\n"
        << "                farm uses 1. Progress goes to stderr; stdout stays the cutflow.\n"
        << "  --partition-target\n"
        << "                approximate events per partition when --threads > 1 (default "
        << kPartitionTargetDefault << ").\n"
        << "                A granularity knob, NOT a correctness knob: the partition count is\n"
        << "                a pure function of (event count, this target) and never of the\n"
        << "                thread count, so the output is byte-identical across --threads for\n"
        << "                any fixed value. Lower it to expose parallelism on a small input.\n";
}

/// \throws std::runtime_error on any malformed or missing argument.
[[nodiscard]] Args parse_args(int argc, char* argv[]) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(flag + " requires a value");
            return argv[++i];
        };
        if (flag == "--input") a.input = value();
        else if (flag == "--output") a.output = value();
        else if (flag == "--config") a.config = value();
        else if (flag == "--target") a.target = value();
        else if (flag == "--run") {
            const std::string v = value();
            try {
                std::size_t pos = 0;
                a.run = std::stoi(v, &pos);
                if (pos != v.size()) throw std::invalid_argument("trailing characters");
            } catch (const std::exception&) {
                throw std::runtime_error("--run wants an integer, got \"" + v + "\"");
            }
        } else if (flag == "--max-events") {
            const std::string v = value();
            unsigned long long n = 0;
            try {
                // stoull, not stoi: "--max-events 5000000000" must be reported
                // as out of range rather than silently become something else.
                // The '-' check is not decoration: stoull happily turns "-1"
                // into ULLONG_MAX, which would then clamp to "all events" -- a
                // typo that silently means the opposite of what it says.
                if (v.find('-') != std::string::npos) throw std::invalid_argument("negative");
                std::size_t pos = 0;
                n = std::stoull(v, &pos);
                if (pos != v.size()) throw std::invalid_argument("trailing characters");
            } catch (const std::exception&) {
                throw std::runtime_error("--max-events wants a positive integer, got \"" + v + "\"");
            }
            if (n == 0) {
                throw std::runtime_error(
                    "--max-events 0 would read nothing and write an empty skim. Omit the flag to "
                    "process the whole file.");
            }
            if (n > std::numeric_limits<unsigned int>::max()) {
                throw std::runtime_error(
                    "--max-events " + v + " exceeds the largest range this program can take (" +
                    std::to_string(std::numeric_limits<unsigned int>::max()) +
                    "). Omit the flag to process the whole file.");
            }
            a.max_events = static_cast<unsigned int>(n);
        } else if (flag == "--qa-ntuple") {
            a.qa_ntuple = value();
        } else if (flag == "--threads") {
            const std::string v = value();
            unsigned long long n = 0;
            try {
                // Same guards as --max-events and as stageB_bin: reject "-1"
                // (stoull turns it into ULLONG_MAX) and trailing junk, so a typo
                // is reported rather than silently reinterpreted.
                if (v.find('-') != std::string::npos) throw std::invalid_argument("negative");
                std::size_t pos = 0;
                n = std::stoull(v, &pos);
                if (pos != v.size()) throw std::invalid_argument("trailing characters");
            } catch (const std::exception&) {
                throw std::runtime_error("--threads wants a positive integer, got \"" + v + "\"");
            }
            if (n == 0) throw std::runtime_error("--threads 0 is meaningless; the minimum is 1 (the default).");
            if (n > std::numeric_limits<unsigned int>::max()) {
                throw std::runtime_error("--threads " + v + " is absurd; the maximum is " +
                                         std::to_string(std::numeric_limits<unsigned int>::max()) + ".");
            }
            a.threads = static_cast<unsigned int>(n);
        } else if (flag == "--partition-target") {
            const std::string v = value();
            unsigned long long n = 0;
            try {
                if (v.find('-') != std::string::npos) throw std::invalid_argument("negative");
                std::size_t pos = 0;
                n = std::stoull(v, &pos);
                if (pos != v.size()) throw std::invalid_argument("trailing characters");
            } catch (const std::exception&) {
                throw std::runtime_error("--partition-target wants a positive integer, got \"" + v + "\"");
            }
            // partition_count() treats 0 as 1, but a user typing 0 has a typo, not
            // an intent -- refuse rather than silently reinterpret. Same as Stage B.
            if (n == 0) throw std::runtime_error("--partition-target 0 is meaningless; it is events per partition.");
            a.partition_target = static_cast<std::size_t>(n);
        } else {
            throw std::runtime_error("unknown argument \"" + flag + "\"");
        }
    }

    if (a.input.empty()) throw std::runtime_error("--input is required");
    if (a.output.empty()) throw std::runtime_error("--output is required");
    if (a.config.empty()) throw std::runtime_error("--config is required");
    if (a.target.empty()) throw std::runtime_error("--target is required");
    // Two RECREATEs on one path leaves whichever closed last, so the slim would
    // be silently replaced by a diagnostic file that looks nothing like it.
    if (a.qa_ntuple.has_value() && *a.qa_ntuple == a.output) {
        throw std::runtime_error(
            "--qa-ntuple and --output name the same file (\"" + a.output +
            "\"). The QA ntuple is a separate diagnostic file, not a variant of the slim.");
    }
    return a;
}

/// Parse --target. Case-insensitive; the spellings match cuts.json's /vertex/targets keys.
[[nodiscard]] pi0::vertex::Target parse_target(const std::string& s) {
    std::string k;
    for (char c : s) k += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (k == "ld2") return pi0::vertex::Target::LD2;
    if (k == "cxc") return pi0::vertex::Target::CxC;
    if (k == "cu") return pi0::vertex::Target::Cu;
    if (k == "sn") return pi0::vertex::Target::Sn;
    throw std::runtime_error("unknown --target \"" + s + "\" (want LD2, CxC, Cu or Sn)");
}

// ===========================================================================
// small geometry helpers
// ===========================================================================

[[nodiscard]] double theta_deg_of(double px, double py, double pz) {
    const double p = std::sqrt(px * px + py * py + pz * pz);
    if (p == 0.0) return 0.0;
    return std::acos(pz / p) * pi0::kRadToDeg;
}

[[nodiscard]] double phi_deg_of(double px, double py) { return std::atan2(py, px) * pi0::kRadToDeg; }

/// REC::Traj `edge` for one particle at one DC layer.
///
/// \return 0.0 when the track has no REC::Traj row for that layer. Zero fails
///         every dc_edge threshold, so a track with no trajectory is rejected
///         rather than admitted -- the safe direction, and visible in the
///         cutflow as a dc_edge failure rather than as a silent absence.
[[nodiscard]] double traj_edge(const RVecS& pindex, const RVecS& detector, const RVecS& layer, const RVecF& edge,
                               std::size_t row, int want_layer) {
    constexpr int kDetectorDc = 6;
    const std::size_t n = std::min({pindex.size(), detector.size(), layer.size(), edge.size()});
    for (std::size_t i = 0; i < n; ++i) {
        if (static_cast<std::size_t>(pindex[i]) == row && detector[i] == kDetectorDc && layer[i] == want_layer) {
            return static_cast<double>(edge[i]);
        }
    }
    return 0.0;
}

// ===========================================================================
// cutflow
// ===========================================================================
//
// A cutflow row's label is RENDERED FROM THE CONFIG, never typed. The old
// analysis printed "Momentum > 0.8 GeV" beside a cut that had always been
// p > 2.0: the label was a string literal that drifted from the constant next
// to it, and every table anyone ever copied out of those logs inherited the
// error. A label computed from the value it describes cannot lie. The six
// electron rows come from electron_cutflow_label(); the rest are formatted here
// from `cuts` for the same reason.

struct CutflowRow {
    std::string label;
    long long survivors{};
};

[[nodiscard]] std::string fmt(double v) {
    std::ostringstream os;
    os << v;
    return os.str();
}

void print_cutflow(const std::vector<CutflowRow>& rows) {
    std::size_t w = 0;
    for (const auto& r : rows) w = std::max(w, r.label.size());

    const long long total = rows.empty() ? 0 : rows.front().survivors;
    std::cout << "\n--- cutflow ---\n";
    for (const auto& r : rows) {
        const double abs_pct = total > 0 ? 100.0 * static_cast<double>(r.survivors) / static_cast<double>(total) : 0.0;
        std::cout << "  " << std::left << std::setw(static_cast<int>(w)) << r.label << "  " << std::right << std::setw(10)
                  << r.survivors << "   " << std::fixed << std::setprecision(2) << std::setw(6) << abs_pct << "% of all\n";
    }
    std::cout << std::defaultfloat;
}

// ===========================================================================
// provenance
// ===========================================================================
//
// NOT OPTIONAL. The analysis this replaces could not reproduce its own
// production: nothing recorded which cuts, which model, or which code had
// produced a given file, so a plot could not be traced back to the thing that
// made it. Everything needed to re-run this invocation goes in the output file.
//
// The Provenance struct and utc_now() live in util/Provenance.hpp, not here.
// They used to be local to this file, which was fine while this was the only
// program that wrote a block. stageB_bin now READS this block back out of the
// slim file and propagates it forward, which makes the layout a contract BETWEEN
// two programs rather than a detail inside one -- and a contract with two
// implementations holds only until somebody edits one of them. Same reasoning,
// and the same day, as util/Sha256.hpp. See that header.

using pi0::util::Provenance;
using pi0::util::utc_now;

// ===========================================================================
// the run number
// ===========================================================================
//
// *** RUN 0 IS NOT A RUN NUMBER -- IT IS AN UNFILLED RUN::config BANK. ***
//
// Measured, not assumed. /Users/mathieuouillon/Documents/tmp/hipo-utils/
// DVKpKmP_006616.hipo is RG-A run 6616 and has 137 tag-0 events. EIGHT of them
// carry RUN_config_run[0] == 0, and THE VERY FIRST EVENT IS ONE OF THEM. (The
// whole file, across all tags, reports run min 0 / max 6616 / mean 6615.57 --
// consistent with ~21 junk rows in 326038.)
//
// So "take RUN_config_run[0] of the first entry" reads a 0, and the model
// lookup then refuses run 0 -- a run the table covers perfectly well -- with a
// message about RG-D. The refusal machinery works; it is being handed garbage.
// Scanning forward for the first REAL run is the fix, and it is not a fudge:
// zero is not in any run range that has ever existed, so a zero carries no
// information except "this bank was not filled for this event".
//
// This does NOT weaken the RG-D refusal. A genuine RG-D run (18305-19131) is
// non-zero, is found by this scan, and is then refused by select_model()
// exactly as intended. Only the "bank not filled" case is skipped past.

/// How far to scan for a filled RUN::config before giving up.
///
/// Not a cut -- nothing physics-related is selected on it, so it is not in
/// cuts.json. It bounds an I/O probe: a file whose first thousand tag-0 events
/// all lack a run number is broken in a way the operator should hear about
/// rather than have guessed around, and --run is the escape hatch.
constexpr ULong64_t kRunProbeMaxEvents = 1000;

/// The file's run number: the first non-zero RUN_config_run[0] in the first
/// kRunProbeMaxEvents tag-0 events.
///
/// Opened as its OWN RHipoDS rather than reusing the main one: the main
/// RDataFrame takes ownership of its datasource by move, and re-running an
/// event loop over a moved-from HIPO datasource is not a behaviour this project
/// has verified. A second file open costs milliseconds and is unambiguous.
///
/// This is a PROBE, not a verification: it reads a prefix of the file, so it
/// cannot prove the file holds a single run. The main event loop checks that --
/// it already reads RUN_config_run for every event, so the check is free there
/// and impossible here.
///
/// \return nullopt if no event in the scanned window carries a non-zero run.
[[nodiscard]] std::optional<int> probe_run_number(const std::string& input) {
    auto ds = std::make_unique<RHipoDS>(input, 0);
    ROOT::RDataFrame df(std::move(ds));
    // RUN_config_run is RVec<int> even though RUN::config is a single-row bank
    // -- see the header comment on n_inspect = 0.
    auto runs = df.Range(0, kRunProbeMaxEvents).Take<RVecI>("RUN_config_run");
    for (const auto& r : *runs) {
        if (!r.empty() && r[0] != 0) return r[0];
    }
    return std::nullopt;
}

// ===========================================================================
// THE EVENT LOGIC, FACTORED OUT SO THERE IS ONE COPY
// ===========================================================================
// process_event() is the SINGLE implementation of the skim's per-event work:
// electron selection, DIS cuts, GBT photon selection, and the QA recording. Both
// the --threads 1 path and each --threads N>1 worker call it, so the selection
// that runs is one function, not a fast copy and a slow copy that could drift --
// the same defect this program's header warns about for the QA values.
//
// It does not touch a TTree. It APPENDS finished rows to per-partition buffers and
// bumps per-partition counters; the caller fills the trees from those buffers,
// SINGLE-THREADED and in entry order. That split is the whole reason the output is
// deterministic under threads: the expensive, order-independent work (scoring) is
// what runs in parallel, and the order-defining work (Fill) never leaves one
// thread. See the file header.

/// One surviving event's slim row: exactly the branches of the "events" tree.
struct SlimRow {
    int b_run{};
    std::int64_t b_event{};
    int b_helicity{};
    double b_q2{}, b_xb{}, b_nu{}, b_w{}, b_y{};
    double b_ex{}, b_ey{}, b_ez{}, b_ee{};
    std::vector<double> b_gpx, b_gpy, b_gpz, b_g_e_gamma_deg;
};

/// One trigger-electron candidate's QA row (pre-cut; rejected ones included).
struct QaElectronRow {
    int run{}, sector{}, failed_at{};
    std::int64_t event{};
    double p{}, theta_deg{}, phi_deg{};
    double vz{}, vz_corrected{}, chi2pid{}, sf{};
    double pcal_e{}, ecin_e{}, ecout_e{}, pcal_lv{}, pcal_lw{};
    double dc_edge_r1{}, dc_edge_r2{}, dc_edge_r3{};
};

/// One photon candidate's QA row (pre-threshold).
struct QaPhotonRow {
    int run{}, prefiltered{}, passed{};
    std::int64_t event{};
    double e{}, theta_deg{}, phi_deg{}, gbt_score{};
    double pcal_e{}, pcal_lv{}, pcal_lw{}, beta{};
};

/// The cutflow / provenance counters, accumulated per partition and merged after.
/// Every field is order-independent -- integer counts and per-key sums -- so the
/// merge is a plain add regardless of partition order (unlike Stage B's sums,
/// which needed compensation; Stage A has none).
struct Counters {
    long long n_all{}, n_has_electron{}, n_electron_pass{};
    long long n_q2{}, n_w{}, n_y{}, n_has_photon{}, n_written{};
    long long n_run_zero{};             ///< events whose RUN::config said 0
    long long n_qa_electron_passed{};   ///< qa_electron rows with failed_at == -1
    std::map<std::string, long long> electron_fail;  ///< stage id -> count
    std::map<int, long long> runs_seen;              ///< non-zero run -> events

    void merge(const Counters& o) {
        n_all += o.n_all;
        n_has_electron += o.n_has_electron;
        n_electron_pass += o.n_electron_pass;
        n_q2 += o.n_q2;
        n_w += o.n_w;
        n_y += o.n_y;
        n_has_photon += o.n_has_photon;
        n_written += o.n_written;
        n_run_zero += o.n_run_zero;
        n_qa_electron_passed += o.n_qa_electron_passed;
        for (const auto& [k, v] : o.electron_fail) electron_fail[k] += v;
        for (const auto& [k, v] : o.runs_seen) runs_seen[k] += v;
    }
};

/// One partition's output: the rows it emitted (in entry order) plus its counters.
struct PartitionBuffers {
    std::vector<SlimRow> slim;
    std::vector<QaElectronRow> qa_e;
    std::vector<QaPhotonRow> qa_g;
    Counters cnt;
};

/// One event's raw columns, COPIED out of RHipoDS's per-event buffers so the event
/// can be scored later, off the reader's thread. Only the --threads N>1 path uses
/// it: the read is sequential (one shared RHipoDS -- see the header on why the read
/// cannot be partitioned), and these buffered events are then scored in parallel.
/// The member names and types mirror process_event()'s parameters exactly.
struct RawEvent {
    RVecI p_pid;
    RVecF p_px, p_py, p_pz, p_vz, p_beta, p_chi2pid;
    RVecS p_status;
    RVecS c_pindex, c_layer, c_sector;
    RVecF c_energy, c_x, c_y, c_z, c_m2u, c_m2v, c_lu, c_lv, c_lw;
    RVecS t_pindex, t_detector, t_layer;
    RVecF t_edge;
    RVecS ev_helicity;
    RVecI rc_run, rc_event;
};

/// Everything process_event() needs that is CONST for the whole run. All of it is
/// safe to share across worker threads: the cuts and the model are immutable, and
/// VzCorrector::correct() is documented const/thread-safe once the variant is set
/// (which happens before any worker starts).
struct EventContext {
    int run{};                                  ///< the FILE's run (not this event's)
    const pi0::Cuts* cuts{};
    const pi0::photonid::FeatureCuts* feat_cuts{};
    pi0::photonid::ModelFn model{};
    const pi0::vertex::VzCorrector* vz{};
    bool vz_corrected{};
    const pi0::vertex::VzTargetCuts* vz_cuts{};
    bool qa_enabled{};
};

/// Process ONE event into `buf`. A faithful move of the former Foreach lambda body
/// -- every value written comes out of the same call the skim decides on -- with
/// the tree Fills replaced by buffer appends and the counters made per-partition.
void process_event(const RVecI& p_pid, const RVecF& p_px, const RVecF& p_py, const RVecF& p_pz, const RVecF& p_vz,
                   const RVecF& p_beta, const RVecF& p_chi2pid, const RVecS& p_status,
                   const RVecS& c_pindex, const RVecS& c_layer, const RVecS& c_sector, const RVecF& c_energy,
                   const RVecF& c_x, const RVecF& c_y, const RVecF& c_z, const RVecF& c_m2u, const RVecF& c_m2v,
                   const RVecF& c_lu, const RVecF& c_lv, const RVecF& c_lw,
                   const RVecS& t_pindex, const RVecS& t_detector, const RVecS& t_layer, const RVecF& t_edge,
                   const RVecS& ev_helicity, const RVecI& rc_run, const RVecI& rc_event,
                   const EventContext& ctx, PartitionBuffers& buf) {
    const pi0::Cuts& cuts = *ctx.cuts;
    ++buf.cnt.n_all;

    // Single-row banks are RVecs too. Guard for empty, always.
    const int this_run = rc_run.empty() ? 0 : rc_run[0];
    const std::int64_t this_event = rc_event.empty() ? 0 : static_cast<std::int64_t>(rc_event[0]);
    if (this_run == 0) ++buf.cnt.n_run_zero; else ++buf.cnt.runs_seen[this_run];

    // HELICITY IS HWP-CORRECTED ALREADY (REC::Event.helicity). 0 == UNDEFINED.
    const int helicity = ev_helicity.empty() ? 0 : static_cast<int>(ev_helicity[0]);

    // ---- trigger electron -------------------------------------
    const auto e_row_opt = pi0::selection::find_trigger_electron(p_pid, p_status, p_px, p_py, p_pz, cuts);
    if (!e_row_opt.has_value()) return;
    ++buf.cnt.n_has_electron;
    const std::size_t e = *e_row_opt;

    const pi0::photonid::CaloMap calo = pi0::photonid::CaloMap::build(
        c_pindex, c_layer, c_sector, c_energy, c_x, c_y, c_z, c_m2u, c_m2v, c_lu, c_lv, c_lw);

    const double ex = p_px[e], ey = p_py[e], ez = p_pz[e];
    const double e_p = std::sqrt(ex * ex + ey * ey + ez * ez);
    const double e_theta_deg = theta_deg_of(ex, ey, ez);
    const double e_phi_deg = phi_deg_of(ex, ey);

    // Calorimeter quantities. A track with no ECAL row gets zeros, which fail the
    // sampling-fraction band and the PCAL fiducial, so it is rejected either way.
    const pi0::photonid::CaloRowData* e_calo = calo.find(e);
    const double e_pcal_e = e_calo ? e_calo->pcal.e : 0.0;
    const double e_ecin_e = e_calo ? e_calo->ecin.e : 0.0;
    const double e_ecout_e = e_calo ? e_calo->ecout.e : 0.0;
    const int e_sector = e_calo ? e_calo->pcal.sector : 0;
    const double e_lv = e_calo ? e_calo->pcal.lv : 0.0;
    const double e_lw = e_calo ? e_calo->pcal.lw : 0.0;
    const double e_sf = e_p > 0.0 ? (e_pcal_e + e_ecin_e + e_ecout_e) / e_p : 0.0;

    // ---- vertex ------------------------------------------------
    // LD2 on the RAW vz; the solid targets on the corrected one (for Cu/Sn the
    // correction IS the target assignment).
    const double e_vz_raw = static_cast<double>(p_vz[e]);
    const double e_vz_used =
        ctx.vz_corrected ? ctx.vz->correct(e_vz_raw, e_p, e_theta_deg, e_phi_deg, e_sector) : e_vz_raw;
    const bool vertex_passed = pi0::vertex::VzCorrector::pass_window(e_vz_used, *ctx.vz_cuts);

    const double edge_r1 = traj_edge(t_pindex, t_detector, t_layer, t_edge, e, 6);
    const double edge_r2 = traj_edge(t_pindex, t_detector, t_layer, t_edge, e, 18);
    const double edge_r3 = traj_edge(t_pindex, t_detector, t_layer, t_edge, e, 36);

    const auto verdict = pi0::selection::pass_electron(static_cast<double>(p_chi2pid[e]), e_p, vertex_passed,
                                                       e_sf, e_sector, e_lv, e_lw, edge_r1, edge_r2, edge_r3,
                                                       cuts);

    // ---- QA row ------------------------------------------------
    // Between the verdict and the return that acts on it, so rejected candidates
    // are here. Every branch is a value pass_electron() was just handed, or the
    // verdict it returned.
    if (ctx.qa_enabled) {
        QaElectronRow r;
        r.run = ctx.run;
        r.event = this_event;
        r.sector = e_sector;
        r.p = e_p;
        r.theta_deg = e_theta_deg;
        r.phi_deg = e_phi_deg;
        r.vz = e_vz_raw;
        // = vz when the target has no correction (LD2); still filled (a NaN would
        // make "corrected" mean two things, and the provenance records which).
        r.vz_corrected = e_vz_used;
        r.chi2pid = static_cast<double>(p_chi2pid[e]);
        r.sf = e_sf;
        r.pcal_e = e_pcal_e;
        r.ecin_e = e_ecin_e;
        r.ecout_e = e_ecout_e;
        r.pcal_lv = e_lv;
        r.pcal_lw = e_lw;
        r.dc_edge_r1 = edge_r1;
        r.dc_edge_r2 = edge_r2;
        r.dc_edge_r3 = edge_r3;
        // The SAME identifier the cutflow counts, mapped to the index of the row it
        // belongs to. -1 means it passed all six.
        r.failed_at = pi0::selection::electron_stage_index(verdict.failed_at);
        if (r.failed_at < 0) ++buf.cnt.n_qa_electron_passed;
        buf.qa_e.push_back(std::move(r));
    }

    if (!verdict.passed) {
        ++buf.cnt.electron_fail[verdict.failed_at];
        return;
    }
    ++buf.cnt.n_electron_pass;

    // ---- DIS ---------------------------------------------------
    // E' = sqrt(p^2 + m_e^2): the electron mass is carried, not dropped.
    const double ee = pi0::energy_from_p(e_p, pi0::kElectronMassGeV);
    const pi0::DisKin dis = pi0::kin::compute_dis(ex, ey, ez, ee, cuts.beam.energy_gev);

    if (!(dis.q2 > cuts.dis.q2_min)) return;
    ++buf.cnt.n_q2;
    if (!(dis.w > cuts.dis.w_min)) return;
    ++buf.cnt.n_w;
    if (!(dis.y < cuts.dis.y_max)) return;
    ++buf.cnt.n_y;

    // ---- photons ------------------------------------------------
    std::vector<double> gpx, gpy, gpz, g_e_gamma_deg;

    const std::size_t n_rows = std::min({p_pid.size(), p_px.size(), p_py.size(), p_pz.size(), p_beta.size()});
    for (std::size_t r = 0; r < n_rows; ++r) {
        if (p_pid[r] != pi0::selection::kPdgPhoton) continue;

        const pi0::photonid::CaloRowData* g_calo = calo.find(r);
        // No ECAL row -> PCAL energy 0 -> the pre-filter's E_PCAL > 0 rejects it and
        // score_fn is never invoked, which is what keeps build_features() from ever
        // seeing a row with no calorimeter data (it throws for one).
        const double g_pcal_e = g_calo ? g_calo->pcal.e : 0.0;
        const double g_lv = g_calo ? g_calo->pcal.lv : 0.0;
        const double g_lw = g_calo ? g_calo->pcal.lw : 0.0;

        const double gx = p_px[r], gy = p_py[r], gz = p_pz[r];

        // MOST CLUSTERS HAVE NO SCORE, AND THAT IS NOT A GAP TO FILL. The pre-filter
        // runs before the callable, so for anything it rejects the GBT is never
        // evaluated and no score exists. The QA row says so with prefiltered == 1
        // and a NaN, rather than a 0 (a legal score).
        double g_score = std::numeric_limits<double>::quiet_NaN();
        bool g_scored = false;

        // The scored overload, not the eager one: the GBT is the expensive part,
        // and the pre-filter exists precisely so it is not evaluated on most
        // clusters. This is what --threads parallelises.
        const bool ok = pi0::selection::pass_photon_scored(
            p_pid[r], gx, gy, gz, g_pcal_e, g_lv, g_lw, static_cast<double>(p_beta[r]),
            [&]() {
                const std::vector<float> feats =
                    pi0::photonid::build_features(r, calo, p_pid, p_px, p_py, p_pz, *ctx.feat_cuts);
                g_score = pi0::photonid::score(feats, ctx.model);
                g_scored = true;
                return g_score;
            },
            cuts);

        // ---- QA row ---------------------------------------------
        // Before the `continue`, so rejected candidates are here.
        if (ctx.qa_enabled) {
            QaPhotonRow r_g;
            r_g.run = ctx.run;
            r_g.event = this_event;
            // The selection's OWN energy and angle, not a second sqrt and acos.
            r_g.e = pi0::selection::photon_energy_gev(gx, gy, gz);
            r_g.theta_deg = pi0::selection::photon_theta_deg(gx, gy, gz);
            // phi is the one branch no cut reads; carried because a fiducial map is
            // drawn in (theta, phi) and this file exists to be plotted.
            r_g.phi_deg = phi_deg_of(gx, gy);
            r_g.gbt_score = g_score;
            r_g.pcal_e = g_pcal_e;
            r_g.pcal_lv = g_lv;
            r_g.pcal_lw = g_lw;
            r_g.beta = static_cast<double>(p_beta[r]);
            r_g.prefiltered = g_scored ? 0 : 1;
            r_g.passed = ok ? 1 : 0;
            buf.qa_g.push_back(std::move(r_g));
        }

        if (!ok) continue;

        gpx.push_back(gx);
        gpy.push_back(gy);
        gpz.push_back(gz);
        // Precomputed here so the pairing stage never needs the electron again.
        g_e_gamma_deg.push_back(pi0::kin::angle_between_deg(gx, gy, gz, ex, ey, ez));
    }

    // >= 1, NOT >= 2 -- single-photon events feed the mixed-event donor pool.
    if (gpx.empty()) return;
    ++buf.cnt.n_has_photon;

    // The FILE's run, not this_run. They differ only when RUN::config was unfilled
    // for this event (this_run == 0), and a 0 in the run branch is a landmine for
    // every downstream per-run normalisation. The single-run check after the loop
    // licenses the substitution; n_run_zero records it.
    SlimRow row;
    row.b_run = ctx.run;
    row.b_event = this_event;
    row.b_helicity = helicity;
    row.b_q2 = dis.q2;
    row.b_xb = dis.xb;
    row.b_nu = dis.nu;
    row.b_w = dis.w;
    row.b_y = dis.y;
    row.b_ex = ex;
    row.b_ey = ey;
    row.b_ez = ez;
    row.b_ee = ee;
    row.b_gpx = std::move(gpx);
    row.b_gpy = std::move(gpy);
    row.b_gpz = std::move(gpz);
    row.b_g_e_gamma_deg = std::move(g_e_gamma_deg);
    buf.slim.push_back(std::move(row));
    ++buf.cnt.n_written;
}

}  // namespace

int main(int argc, char* argv[]) {
    // ---- arguments -------------------------------------------------------
    Args args;
    try {
        args = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n\n";
        print_usage(argv[0]);
        return 1;
    }

    try {
        // ---- configuration ------------------------------------------------
        const pi0::Cuts cuts = pi0::Cuts::load(args.config);
        const pi0::vertex::Target target = parse_target(args.target);

        // One polarity string in cuts.json; two enums, because pi0::vertex and
        // pi0::selection are deliberately independent libraries and neither may
        // depend on the other's types. Both are parsed from the SAME string, so
        // they cannot disagree.
        const pi0::selection::Polarity sf_polarity = pi0::selection::polarity_from_string(cuts.beam.polarity);
        const pi0::vertex::Polarity vz_polarity = (sf_polarity == pi0::selection::Polarity::Outbending)
                                                      ? pi0::vertex::Polarity::Outbending
                                                      : pi0::vertex::Polarity::Inbending;

        const std::string config_sha = pi0::util::sha256_file(args.config);

        // ---- vertex correction --------------------------------------------
        // LD2 is uncorrected by construction (one long cryotarget cell, no foil
        // peaks to separate), so it has no variant and set_variant() would throw.
        pi0::vertex::VzCorrector vz = pi0::vertex::VzCorrector::load(PI0_VZ_PARAMS_FILE);
        const pi0::vertex::VzTargetCuts& vz_cuts = cuts.vertex.for_target(target);
        const bool vz_corrected = vz_cuts.correction_enabled;
        if (vz_corrected) vz.set_variant(target, vz_polarity);

        // ---- the run number and the GBT model -----------------------------
        int run = 0;
        std::string run_source;
        if (args.run.has_value()) {
            run = *args.run;
            run_source = "--run (command line override)";
        } else {
            const std::optional<int> probed = probe_run_number(args.input);
            if (!probed.has_value()) {
                std::cerr << "error: no run number in the first " << kRunProbeMaxEvents << " tag-0 events of "
                          << args.input << "\n"
                          << "       Every one of them has RUN::config absent or reporting run 0, which means the\n"
                          << "       bank was never filled. The GBT model is keyed on the run number, so there is\n"
                          << "       nothing to select on. Pass --run <N> if you know the run independently.\n";
                return 2;
            }
            run = *probed;
            run_source = "first non-zero RUN_config_run[0]";
        }

        pi0::photonid::ModelFn model = nullptr;
        try {
            model = pi0::photonid::select_model(run, cuts.photon.gbt_pass, cuts.photon.allow_rga_fallback);
        } catch (const std::exception& e) {
            // This is the designed refusal, not a crash. cuts.json's
            // photon._allow_rga_fallback_comment explains why RG-D must not be
            // scored by an RG-A inbending model behind the analyst's back.
            std::cerr << "error: no GBT photon model covers this file.\n"
                      << "       " << e.what() << "\n"
                      << "       run                = " << run << "  (from " << run_source << ")\n"
                      << "       gbt_pass           = " << cuts.photon.gbt_pass << "\n"
                      << "       allow_rga_fallback = false\n"
                      << "       This is the intended behaviour, not a bug to route around. Setting\n"
                      << "       photon.allow_rga_fallback = true in " << args.config << " applies the RG-A\n"
                      << "       inbending pass-1 model to this data and stamps the fallback into the\n"
                      << "       output's provenance. Do that only for a study you intend to report as such.\n";
            return 3;
        }
        const std::string model_name(pi0::photonid::model_name(model));
        const bool fallback_used = pi0::photonid::fallback_used(run, cuts.photon.gbt_pass);

        std::cout << "stageA_skim\n"
                  << "  input      : " << args.input << '\n'
                  << "  output     : " << args.output << '\n'
                  << "  config     : " << args.config << "  (sha256 " << config_sha.substr(0, 16) << "...)\n"
                  << "  target     : " << pi0::vertex::to_string(target) << '\n'
                  << "  polarity   : " << cuts.beam.polarity << '\n'
                  << "  run        : " << run << "  (from " << run_source << ")\n"
                  << "  max events : "
                  << (args.max_events.has_value()
                          ? std::to_string(*args.max_events) + "   *** TRUNCATED RUN (--max-events) ***"
                          : std::string("all"))
                  << '\n'
                  << "  GBT model  : " << model_name << (fallback_used ? "   *** RGA FALLBACK IN USE ***" : "") << '\n';
        if (fallback_used) {
            std::cout << "  WARNING: this run matches no trained model. Every photon below is being\n"
                      << "           identified by a model trained on different data. The output's\n"
                      << "           provenance records this; any plot made from it must say so.\n";
        }
        if (args.qa_ntuple.has_value()) {
            std::cout << "  qa ntuple  : " << *args.qa_ntuple << "\n"
                      << "           DIAGNOSTIC ONLY -- rows are CANDIDATES, not events. Nothing in it is a\n"
                      << "           yield. The slim above is unaffected by its presence.\n";
        }

        // ---- output --------------------------------------------------------
        // Created before the loop so a bad path fails now, not after the pass.
        std::unique_ptr<TFile> out(TFile::Open(args.output.c_str(), "RECREATE"));
        if (out == nullptr || out->IsZombie()) throw std::runtime_error("cannot create output file: " + args.output);

        int b_run = 0;
        std::int64_t b_event = 0;
        int b_helicity = 0;
        double b_q2 = 0, b_xb = 0, b_nu = 0, b_w = 0, b_y = 0;
        double b_ex = 0, b_ey = 0, b_ez = 0, b_ee = 0;
        std::vector<double> b_gpx, b_gpy, b_gpz, b_g_e_gamma_deg;

        // A TTree, not an RNTuple: uproot reads TTree reliably today, and the
        // downstream stages are Python.
        auto* tree = new TTree("events", "stage A slim skim: DIS electron + selected photons");
        tree->Branch("run", &b_run, "run/I");
        tree->Branch("event", &b_event, "event/L");
        tree->Branch("helicity", &b_helicity, "helicity/I");
        tree->Branch("q2", &b_q2, "q2/D");
        tree->Branch("xb", &b_xb, "xb/D");
        tree->Branch("nu", &b_nu, "nu/D");
        tree->Branch("w", &b_w, "w/D");
        tree->Branch("y", &b_y, "y/D");
        tree->Branch("ex", &b_ex, "ex/D");
        tree->Branch("ey", &b_ey, "ey/D");
        tree->Branch("ez", &b_ez, "ez/D");
        tree->Branch("ee", &b_ee, "ee/D");
        tree->Branch("gpx", &b_gpx);
        tree->Branch("gpy", &b_gpy);
        tree->Branch("gpz", &b_gpz);
        tree->Branch("g_e_gamma_deg", &b_g_e_gamma_deg);

        // ---- QA ntuple (--qa-ntuple), diagnostic only ------------------------
        // =====================================================================
        // ONE ROW PER CANDIDATE. NOT A YIELD. NOT PER EVENT.
        // =====================================================================
        // Flat scalar branches in TTrees, the same house pattern and for the same
        // reason as stageB_bin's four output trees: uproot reads them into numpy
        // with no ROOT installed and no dictionary.
        //
        // WHAT THESE ROWS ARE, exactly, because a QA file is read by whoever is
        // debugging at the time and the counts must not be mistaken:
        //
        //   qa_electron -- one row per TRIGGER ELECTRON, i.e. per event that had
        //     one at all (find_trigger_electron returns at most one row per
        //     event). The row is filled BEFORE the verdict is acted on, so the
        //     rejected ones are here. `failed_at` is the index of the cut that
        //     rejected it, -1 if it survived all six.
        //
        //   qa_photon -- one row per pid == 22 row in REC::Particle, filled
        //     before the GBT threshold is applied. THESE EXIST ONLY FOR EVENTS
        //     THAT ALREADY PASSED THE ELECTRON AND DIS CUTS: the photon loop is
        //     downstream of both, and hoisting it earlier would mean scoring the
        //     GBT on events the skim discards. So this tree is a diagnostic of
        //     the photon selection GIVEN a DIS event, not of photons in the file.
        //
        // Every value here comes out of the call that decided, never out of a
        // second computation of it. See the header comment.
        // =====================================================================
        std::unique_ptr<TFile> qa_out;
        TTree* qa_e_tree = nullptr;
        TTree* qa_g_tree = nullptr;

        int qe_run = 0, qe_sector = 0, qe_failed_at = 0;
        std::int64_t qe_event = 0;
        double qe_p = 0, qe_theta_deg = 0, qe_phi_deg = 0;
        double qe_vz = 0, qe_vz_corrected = 0, qe_chi2pid = 0, qe_sf = 0;
        double qe_pcal_e = 0, qe_ecin_e = 0, qe_ecout_e = 0, qe_pcal_lv = 0, qe_pcal_lw = 0;
        double qe_dc_edge_r1 = 0, qe_dc_edge_r2 = 0, qe_dc_edge_r3 = 0;

        int qg_run = 0, qg_passed = 0, qg_prefiltered = 0;
        std::int64_t qg_event = 0;
        double qg_e = 0, qg_theta_deg = 0, qg_phi_deg = 0, qg_gbt_score = 0;
        double qg_pcal_e = 0, qg_pcal_lv = 0, qg_pcal_lw = 0, qg_beta = 0;

        // qa_electron rows with failed_at == -1 are counted (in Counters, per
        // partition, by process_event) independently of n_electron_pass, so the two
        // can be checked against each other after the loop. See the check there.

        if (args.qa_ntuple.has_value()) {
            // Opened here, before the loop, for the same reason the slim is: a
            // bad path must fail now rather than after the pass.
            qa_out.reset(TFile::Open(args.qa_ntuple->c_str(), "RECREATE"));
            if (qa_out == nullptr || qa_out->IsZombie()) {
                throw std::runtime_error("cannot create QA ntuple file: " + *args.qa_ntuple);
            }
            // TTree attaches itself to gDirectory. TFile::Open has already left the
            // QA file current, but this program now has two files open and the
            // slim's `tree` belongs to the other one, so the cd() is stated rather
            // than inherited.
            qa_out->cd();

            qa_e_tree = new TTree("qa_electron",
                                  "DIAGNOSTIC ONLY. One row per trigger-electron candidate, PRE-CUT: filled "
                                  "before the verdict is acted on, so rejected candidates are here. failed_at "
                                  "indexes qa_electron_stages; -1 = passed all six");
            qa_e_tree->Branch("run", &qe_run, "run/I");
            qa_e_tree->Branch("event", &qe_event, "event/L");
            qa_e_tree->Branch("sector", &qe_sector, "sector/I");
            qa_e_tree->Branch("p", &qe_p, "p/D");
            qa_e_tree->Branch("theta_deg", &qe_theta_deg, "theta_deg/D");
            qa_e_tree->Branch("phi_deg", &qe_phi_deg, "phi_deg/D");
            qa_e_tree->Branch("vz", &qe_vz, "vz/D");
            qa_e_tree->Branch("vz_corrected", &qe_vz_corrected, "vz_corrected/D");
            qa_e_tree->Branch("chi2pid", &qe_chi2pid, "chi2pid/D");
            qa_e_tree->Branch("sf", &qe_sf, "sf/D");
            qa_e_tree->Branch("pcal_e", &qe_pcal_e, "pcal_e/D");
            qa_e_tree->Branch("ecin_e", &qe_ecin_e, "ecin_e/D");
            qa_e_tree->Branch("ecout_e", &qe_ecout_e, "ecout_e/D");
            qa_e_tree->Branch("pcal_lv", &qe_pcal_lv, "pcal_lv/D");
            qa_e_tree->Branch("pcal_lw", &qe_pcal_lw, "pcal_lw/D");
            qa_e_tree->Branch("dc_edge_r1", &qe_dc_edge_r1, "dc_edge_r1/D");
            qa_e_tree->Branch("dc_edge_r2", &qe_dc_edge_r2, "dc_edge_r2/D");
            qa_e_tree->Branch("dc_edge_r3", &qe_dc_edge_r3, "dc_edge_r3/D");
            qa_e_tree->Branch("failed_at", &qe_failed_at, "failed_at/I");

            qa_g_tree = new TTree("qa_photon",
                                  "DIAGNOSTIC ONLY. One row per pid == 22 candidate, PRE-THRESHOLD. Rows exist "
                                  "ONLY for events that already passed the electron and DIS cuts. gbt_score is "
                                  "NaN where prefiltered == 1: the GBT was never evaluated, so there is no score");
            qa_g_tree->Branch("run", &qg_run, "run/I");
            qa_g_tree->Branch("event", &qg_event, "event/L");
            qa_g_tree->Branch("e", &qg_e, "e/D");
            qa_g_tree->Branch("theta_deg", &qg_theta_deg, "theta_deg/D");
            qa_g_tree->Branch("phi_deg", &qg_phi_deg, "phi_deg/D");
            qa_g_tree->Branch("gbt_score", &qg_gbt_score, "gbt_score/D");
            qa_g_tree->Branch("pcal_e", &qg_pcal_e, "pcal_e/D");
            qa_g_tree->Branch("pcal_lv", &qg_pcal_lv, "pcal_lv/D");
            qa_g_tree->Branch("pcal_lw", &qg_pcal_lw, "pcal_lw/D");
            qa_g_tree->Branch("beta", &qg_beta, "beta/D");
            qa_g_tree->Branch("prefiltered", &qg_prefiltered, "prefiltered/I");
            qa_g_tree->Branch("passed", &qg_passed, "passed/I");

            // --- qa_electron_stages: the failed_at -> cut mapping -------------
            // Written into the file rather than left for the Python to hard-code.
            // A plotting script carrying its own list of the six cuts is the
            // stale-label defect this program spends a section warning about,
            // re-created in a second language. This is built by iterating
            // kElectronStages, which is what the cutflow iterates, so an axis
            // drawn from it cannot name the cuts in an order the skim did not
            // apply.
            //
            // `label` is rendered by the same electron_cutflow_label() the cutflow
            // prints, so a QA figure's tick text and the cutflow's row are one
            // string from one source. `name` is the stable identifier: it carries
            // no threshold, so it cannot go stale when one moves.
            {
                int b_index = 0;
                std::array<char, 64> b_name{};
                std::array<char, 256> b_label{};
                auto* t = new TTree("qa_electron_stages",
                                    "The failed_at -> cut mapping for qa_electron, in the order the cuts are "
                                    "applied. Built from selection::kElectronStages, the array the cutflow "
                                    "iterates; `label` is rendered from the config by the same function");
                // /C branches: one NUL-terminated string per entry, which uproot
                // reads with no dictionary. snprintf bounds the copies, so a label
                // longer than the buffer truncates rather than overruns.
                t->Branch("index", &b_index, "index/I");
                t->Branch("name", b_name.data(), "name/C");
                t->Branch("label", b_label.data(), "label/C");

                for (const char* stage : pi0::selection::kElectronStages) {
                    b_index = pi0::selection::electron_stage_index(stage);
                    const std::string label = pi0::selection::electron_cutflow_label(stage, cuts);
                    b_name.fill('\0');
                    b_label.fill('\0');
                    std::snprintf(b_name.data(), b_name.size(), "%s", stage);
                    std::snprintf(b_label.data(), b_label.size(), "%s", label.c_str());
                    t->Fill();
                }
                t->Write();
            }
        }

        // ---- input ---------------------------------------------------------
        auto ds = std::make_unique<RHipoDS>(args.input, 0);
        ROOT::RDataFrame df(std::move(ds));

        // Fail loudly and in one go on a missing bank, rather than letting
        // RDataFrame throw "Unknown column" for whichever it hits first.
        const std::vector<std::string> needed = {
            "REC_Particle_pid",       "REC_Particle_px",         "REC_Particle_py",       "REC_Particle_pz",
            "REC_Particle_vz",        "REC_Particle_beta",       "REC_Particle_chi2pid",  "REC_Particle_status",
            "REC_Calorimeter_pindex", "REC_Calorimeter_layer",   "REC_Calorimeter_sector", "REC_Calorimeter_energy",
            "REC_Calorimeter_x",      "REC_Calorimeter_y",       "REC_Calorimeter_z",     "REC_Calorimeter_m2u",
            "REC_Calorimeter_m2v",    "REC_Calorimeter_lu",      "REC_Calorimeter_lv",    "REC_Calorimeter_lw",
            "REC_Traj_pindex",        "REC_Traj_detector",       "REC_Traj_layer",        "REC_Traj_edge",
            "REC_Event_helicity",     "RUN_config_run",          "RUN_config_event",
        };
        {
            const std::vector<std::string> have = df.GetColumnNames();
            std::vector<std::string> missing;
            for (const auto& c : needed) {
                if (std::find(have.begin(), have.end(), c) == have.end()) missing.push_back(c);
            }
            if (!missing.empty()) {
                std::ostringstream os;
                os << "the input is missing " << missing.size() << " required column(s):";
                for (const auto& m : missing) os << "\n  " << m;
                os << "\n(HIPO \"::\" and \".\" each become \"_\": REC::Particle.px -> REC_Particle_px)";
                throw std::runtime_error(os.str());
            }
        }

        // The neighbour window the classifier was trained on. Read from cuts.json
        // rather than hard-coded -- but see Features.hpp: editing these three does
        // NOT retune the classifier, it feeds it a distribution it has never seen.
        const pi0::photonid::FeatureCuts feat_cuts{cuts.photon.min_energy_gev, cuts.photon.theta_min_deg,
                                                   cuts.photon.theta_max_deg};

        // ---- entry count and progress -------------------------------------
        // The total the loop will visit: the file's tag-0 entry count, capped by
        // --max-events. Queried from a throwaway RHipoDS so the datasource feeding
        // `df` (and each worker's own) is untouched. It is the progress bar's
        // denominator and, when --threads > 1, the span the partitions tile.
        unsigned long file_entries = 0;
        {
            auto probe_ds = std::make_unique<RHipoDS>(args.input, 0);
            file_entries = probe_ds->GetEntries();
        }
        const std::size_t total_entries =
            args.max_events.has_value()
                ? std::min<std::size_t>(static_cast<std::size_t>(file_entries), *args.max_events)
                : static_cast<std::size_t>(file_entries);

        // The read-only context every event needs. All of it is safe to share
        // across worker threads: cuts and model are immutable, and
        // VzCorrector::correct() is const/thread-safe once the variant is set
        // (which happened above, before any worker exists).
        EventContext ctx;
        ctx.run = run;
        ctx.cuts = &cuts;
        ctx.feat_cuts = &feat_cuts;
        ctx.model = model;
        ctx.vz = &vz;
        ctx.vz_corrected = vz_corrected;
        ctx.vz_cuts = &vz_cuts;
        ctx.qa_enabled = args.qa_ntuple.has_value();

        // ---- the trees are filled HERE, single-threaded, in entry order ------
        // These take a finished row and Fill the corresponding tree. They are the
        // ONLY callers of TTree::Fill, and only the main thread ever runs them, so
        // the "Fill is not thread-safe" hazard never arises for any --threads.
        auto fill_slim = [&](const SlimRow& r) {
            b_run = r.b_run;
            b_event = r.b_event;
            b_helicity = r.b_helicity;
            b_q2 = r.b_q2;
            b_xb = r.b_xb;
            b_nu = r.b_nu;
            b_w = r.b_w;
            b_y = r.b_y;
            b_ex = r.b_ex;
            b_ey = r.b_ey;
            b_ez = r.b_ez;
            b_ee = r.b_ee;
            b_gpx = r.b_gpx;
            b_gpy = r.b_gpy;
            b_gpz = r.b_gpz;
            b_g_e_gamma_deg = r.b_g_e_gamma_deg;
            tree->Fill();
        };
        auto fill_qa_e = [&](const QaElectronRow& r) {
            qe_run = r.run;
            qe_event = r.event;
            qe_sector = r.sector;
            qe_p = r.p;
            qe_theta_deg = r.theta_deg;
            qe_phi_deg = r.phi_deg;
            qe_vz = r.vz;
            qe_vz_corrected = r.vz_corrected;
            qe_chi2pid = r.chi2pid;
            qe_sf = r.sf;
            qe_pcal_e = r.pcal_e;
            qe_ecin_e = r.ecin_e;
            qe_ecout_e = r.ecout_e;
            qe_pcal_lv = r.pcal_lv;
            qe_pcal_lw = r.pcal_lw;
            qe_dc_edge_r1 = r.dc_edge_r1;
            qe_dc_edge_r2 = r.dc_edge_r2;
            qe_dc_edge_r3 = r.dc_edge_r3;
            qe_failed_at = r.failed_at;
            qa_e_tree->Fill();
        };
        auto fill_qa_g = [&](const QaPhotonRow& r) {
            qg_run = r.run;
            qg_event = r.event;
            qg_e = r.e;
            qg_theta_deg = r.theta_deg;
            qg_phi_deg = r.phi_deg;
            qg_gbt_score = r.gbt_score;
            qg_pcal_e = r.pcal_e;
            qg_pcal_lv = r.pcal_lv;
            qg_pcal_lw = r.pcal_lw;
            qg_beta = r.beta;
            qg_prefiltered = r.prefiltered;
            qg_passed = r.passed;
            qa_g_tree->Fill();
        };
        // Drain one partition's buffered rows into the trees, in entry order.
        // Called in partition-index order, so the whole file is filled in entry
        // order -- which is the byte-for-byte identity the slim's downstream pool
        // depends on.
        auto drain = [&](PartitionBuffers& b) {
            for (const SlimRow& r : b.slim) fill_slim(r);
            if (qa_e_tree != nullptr)
                for (const QaElectronRow& r : b.qa_e) fill_qa_e(r);
            if (qa_g_tree != nullptr)
                for (const QaPhotonRow& r : b.qa_g) fill_qa_g(r);
        };

        // Progress on STDERR (stdout carries the cutflow and provenance). The bar
        // is TTY-aware and thread-safe, so every worker may call add() concurrently.
        pi0::Progress progress("stageA skim", static_cast<std::int64_t>(total_entries), std::cerr);

        const unsigned int n_threads = std::max(1u, args.threads);

        // The counters this run produced, filled by whichever path runs below.
        Counters counters;

        // ---- the event loop -------------------------------------------------
        // --max-events is a PREFIX of the file (the first N tag-0 entries), not a
        // sample; that is why n_all, the cutflow's first row and the provenance all
        // say so. The output is IDENTICAL across --threads: the two paths below emit
        // exactly the same rows in exactly entry order.
        if (n_threads == 1) {
            // ---- the pre-MT single-pass path (this is the farm path) ---------
            // One RHipoDS + RDataFrame + Foreach filling the trees directly, in
            // entry order, with O(1) buffering (the per-event scratch is drained and
            // cleared each event). Byte-for-byte what this program did before MT.
            // Range is single-thread only, which is exactly this path.
            ROOT::RDF::RNode node = df;
            if (args.max_events.has_value()) node = node.Range(0u, *args.max_events);

            PartitionBuffers scratch;  // cnt accumulates across events; rows are per-event
            node.Foreach(
                [&](const RVecI& p_pid, const RVecF& p_px, const RVecF& p_py, const RVecF& p_pz, const RVecF& p_vz,
                    const RVecF& p_beta, const RVecF& p_chi2pid, const RVecS& p_status,
                    const RVecS& c_pindex, const RVecS& c_layer, const RVecS& c_sector, const RVecF& c_energy,
                    const RVecF& c_x, const RVecF& c_y, const RVecF& c_z, const RVecF& c_m2u, const RVecF& c_m2v,
                    const RVecF& c_lu, const RVecF& c_lv, const RVecF& c_lw,
                    const RVecS& t_pindex, const RVecS& t_detector, const RVecS& t_layer, const RVecF& t_edge,
                    const RVecS& ev_helicity, const RVecI& rc_run, const RVecI& rc_event) {
                    process_event(p_pid, p_px, p_py, p_pz, p_vz, p_beta, p_chi2pid, p_status, c_pindex, c_layer,
                                  c_sector, c_energy, c_x, c_y, c_z, c_m2u, c_m2v, c_lu, c_lv, c_lw, t_pindex,
                                  t_detector, t_layer, t_edge, ev_helicity, rc_run, rc_event, ctx, scratch);
                    drain(scratch);
                    scratch.slim.clear();
                    scratch.qa_e.clear();
                    scratch.qa_g.clear();
                    progress.add();
                },
                needed);
            counters = std::move(scratch.cnt);
        } else {
            // ---- the opt-in multi-threaded path (bounded runs) ---------------
            // Read the events ONCE, sequentially, through the one shared RHipoDS --
            // the read cannot be partitioned across readers without paying RHipoDS's
            // whole-file record-index scan once PER reader (see the file header) --
            // buffering each event's raw columns. Then SCORE the buffer in parallel:
            // a FIXED number of contiguous index-partitions decided by the event
            // count, run in waves of `lanes` worker threads, each running
            // process_event over its slice into its own buffers. The partitions are
            // drained into the trees in index order, so the emission order is the
            // single-thread entry order.
            ROOT::EnableThreadSafety();

            // ---- phase 1: sequential read into the raw-event buffer ----------
            // The same Foreach the --threads 1 path runs, but instead of scoring it
            // COPIES each event's columns out of RHipoDS's per-event buffers so they
            // survive past the callback. progress tracks this read: on a HIPO-scan-
            // bound file this is where the wall time goes.
            std::vector<RawEvent> raw;
            raw.reserve(total_entries);
            {
                ROOT::RDF::RNode node = df;
                if (args.max_events.has_value()) node = node.Range(0u, *args.max_events);
                node.Foreach(
                    [&](const RVecI& p_pid, const RVecF& p_px, const RVecF& p_py, const RVecF& p_pz,
                        const RVecF& p_vz, const RVecF& p_beta, const RVecF& p_chi2pid, const RVecS& p_status,
                        const RVecS& c_pindex, const RVecS& c_layer, const RVecS& c_sector, const RVecF& c_energy,
                        const RVecF& c_x, const RVecF& c_y, const RVecF& c_z, const RVecF& c_m2u, const RVecF& c_m2v,
                        const RVecF& c_lu, const RVecF& c_lv, const RVecF& c_lw,
                        const RVecS& t_pindex, const RVecS& t_detector, const RVecS& t_layer, const RVecF& t_edge,
                        const RVecS& ev_helicity, const RVecI& rc_run, const RVecI& rc_event) {
                        raw.push_back(RawEvent{p_pid, p_px, p_py, p_pz, p_vz, p_beta, p_chi2pid, p_status, c_pindex,
                                               c_layer, c_sector, c_energy, c_x, c_y, c_z, c_m2u, c_m2v, c_lu, c_lv,
                                               c_lw, t_pindex, t_detector, t_layer, t_edge, ev_helicity, rc_run,
                                               rc_event});
                        progress.add();
                    },
                    needed);
            }

            // ---- phase 2 + 3: score in parallel, fill in index order ---------
            const std::size_t n_events = raw.size();
            const std::size_t n_part = pi0::partition_count(n_events, args.partition_target);
            const unsigned int lanes = static_cast<unsigned int>(std::min<std::size_t>(n_threads, n_part));

            std::cout << "  threads    : " << n_threads << "  (" << n_part << " data-decided partition(s) over "
                      << n_events << " buffered event(s); output is deterministic, byte-identical to --threads 1)\n";

            // Score ONE partition's contiguous slice of the raw buffer into `out`.
            // Pure computation over already-read events -- no ROOT, no RHipoDS, and
            // no shared mutable state except `out` -- so it parallelises cleanly.
            auto run_partition = [&](std::size_t part_idx, PartitionBuffers& out) {
                const pi0::PartitionRange pr = pi0::partition_range(part_idx, n_part, n_events);
                for (std::size_t i = pr.begin; i < pr.end; ++i) {
                    const RawEvent& e = raw[i];
                    process_event(e.p_pid, e.p_px, e.p_py, e.p_pz, e.p_vz, e.p_beta, e.p_chi2pid, e.p_status,
                                  e.c_pindex, e.c_layer, e.c_sector, e.c_energy, e.c_x, e.c_y, e.c_z, e.c_m2u,
                                  e.c_m2v, e.c_lu, e.c_lv, e.c_lw, e.t_pindex, e.t_detector, e.t_layer, e.t_edge,
                                  e.ev_helicity, e.rc_run, e.rc_event, ctx, out);
                }
            };

            for (std::size_t wave_start = 0; wave_start < n_part; wave_start += lanes) {
                const std::size_t wave_end = std::min<std::size_t>(wave_start + lanes, n_part);
                const std::size_t wcount = wave_end - wave_start;

                std::vector<PartitionBuffers> parts(wcount);
                // A worker's throw must not cross the thread boundary; capture it and
                // rethrow on this thread after the join.
                std::vector<std::exception_ptr> errs(wcount);

                if (wcount == 1) {
                    // Single lane -- the tail wave, or a --threads that resolved to
                    // one lane. Run inline; no std::thread.
                    try {
                        run_partition(wave_start, parts[0]);
                    } catch (...) {
                        errs[0] = std::current_exception();
                    }
                } else {
                    std::vector<std::thread> workers;
                    workers.reserve(wcount);
                    for (std::size_t k = 0; k < wcount; ++k) {
                        const std::size_t part_idx = wave_start + k;
                        workers.emplace_back([&, k, part_idx]() {
                            try {
                                run_partition(part_idx, parts[k]);
                            } catch (...) {
                                errs[k] = std::current_exception();
                            }
                        });
                    }
                    for (std::thread& t : workers) t.join();
                }

                for (std::size_t k = 0; k < wcount; ++k) {
                    if (errs[k]) std::rethrow_exception(errs[k]);
                }
                // INDEX ORDER: drain and merge partition wave_start, wave_start+1,
                // ... so the trees see rows in the single-thread entry order and the
                // counters (all order-independent) are summed.
                for (std::size_t k = 0; k < wcount; ++k) {
                    drain(parts[k]);
                    counters.merge(parts[k].cnt);
                }
            }
        }
        progress.finish();

        // ---- unpack the merged counters into the names the rest of main uses ---
        // The single-run check below reads runs_seen; the provenance and the cutflow
        // read the rest. electron_fail stays non-const: the cutflow indexes it with
        // operator[] (which inserts a 0 for a stage that never failed, as intended).
        const long long n_all = counters.n_all;
        const long long n_has_electron = counters.n_has_electron;
        const long long n_electron_pass = counters.n_electron_pass;
        const long long n_q2 = counters.n_q2;
        const long long n_w = counters.n_w;
        const long long n_y = counters.n_y;
        const long long n_has_photon = counters.n_has_photon;
        const long long n_written = counters.n_written;
        const long long n_run_zero = counters.n_run_zero;
        const long long n_qa_electron_passed = counters.n_qa_electron_passed;
        std::map<std::string, long long>& electron_fail = counters.electron_fail;
        const std::map<int, long long>& runs_seen = counters.runs_seen;

        // ---- the single-run assumption ---------------------------------------
        // Checked before ANYTHING is written. The model was selected once for the
        // whole file; if the file actually spans runs, some events were scored by
        // a model chosen for a different one, and the provenance header -- which
        // records a single run and a single model -- would be a lie. Refuse, and
        // take the half-written file with us: a file that exists is a file
        // somebody will eventually plot.
        if (runs_seen.size() > 1) {
            std::ostringstream os;
            os << "the input spans " << runs_seen.size() << " runs, but the GBT model is selected once per file:";
            for (const auto& [r, n] : runs_seen) os << "\n  run " << r << " : " << n << " events";
            os << "\nSplit the file by run and skim each part separately. (One input HIPO file per invocation is"
                  "\nthe contract; entry numbering restarts per file and the provenance records one run.)";
            out->Close();
            std::remove(args.output.c_str());
            // The QA file goes with it, and for the identical reason: a file that
            // exists is a file somebody will eventually plot. Its rows carry the
            // model chosen for one run applied to candidates from several, which
            // is exactly the confusion this refusal exists to prevent.
            if (qa_out != nullptr) {
                qa_out->Close();
                std::remove(args.qa_ntuple->c_str());
            }
            std::cerr << "error: " << os.str() << '\n';
            return 5;
        }
        if (!runs_seen.empty() && runs_seen.begin()->first != run) {
            // Only reachable with --run: the probe takes its answer from the same
            // column this map is built from. An override that contradicts the data
            // is worth saying out loud, but it is the operator's call, so it warns
            // rather than refuses -- overriding a wrong run number is exactly what
            // --run is for.
            std::cout << "\nWARNING: --run says " << run << " but the file's events report run "
                      << runs_seen.begin()->first << ".\n"
                      << "         The model was chosen for " << run << " and the provenance records the override.\n";
        }

        // ---- provenance ------------------------------------------------------
        Provenance prov;
        prov.add("stageA_skim.code_version", PI0_CODE_VERSION);
        prov.add("stageA_skim.created_utc", utc_now());
        prov.add("input.path", args.input);
        prov.add("config.path", args.config);
        prov.add("config.sha256", config_sha);
        prov.add("target", pi0::vertex::to_string(target));
        prov.add("polarity", cuts.beam.polarity);
        prov.add("beam.energy_gev", fmt(cuts.beam.energy_gev));
        prov.add("run", std::to_string(run));
        prov.add("run.source", run_source);
        prov.add("gbt.model", model_name);
        prov.add("gbt.pass", std::to_string(cuts.photon.gbt_pass));
        prov.add("gbt.threshold", fmt(cuts.photon.gbt_threshold));
        prov.add("gbt.allow_rga_fallback", cuts.photon.allow_rga_fallback ? "true" : "false");
        // The single most important line in this block. A file scored by a model
        // trained on other data must be identifiable as such after the fact.
        prov.add("gbt.fallback_used", fallback_used ? "TRUE -- photons scored by a model trained on OTHER data"
                                                    : "false");
        prov.add("vertex.correction_applied", vz_corrected ? "true" : "false (LD2 uses raw vz by construction)");
        prov.add("vertex.variant", vz_corrected ? vz.active_variant() : std::string("none"));
        prov.add("vz_params.path", PI0_VZ_PARAMS_FILE);
        prov.add("helicity.convention",
                 "REC::Event.helicity[0] -- HWP-CORRECTED. Per the CLAS12 bank definition, "
                 "REC::Event.helicity is 'online-delay-corrected helicity, with HWP-correction "
                 "(0=UDF)', as opposed to REC::Event.helicityRaw which is without. The half-wave-plate "
                 "state is applied upstream by the cooking, so no HWP table is needed here. "
                 "Caveat: the sign is therefore an INHERITED assumption -- no HWP-in/HWP-out closure "
                 "test exists in this chain. Settle it before publishing a signed asymmetry. "
                 "helicity == 0 means undefined and is dropped downstream.");
        prov.add("photon.min_multiplicity", "1 (deliberate: single-photon events feed the mixed-event donor pool)");
        // NOT cosmetic, and not omittable when the flag is absent. events.input_tag0
        // is the count this program READ, which is the file's total only if it
        // read all of it -- so the two keys must be read together, and a key
        // that is present only sometimes is a key nobody checks. "all" is
        // recorded explicitly so that "this was a full pass" is an assertion in
        // the file rather than an inference from a missing field.
        prov.add("events.max_events_requested",
                 args.max_events.has_value()
                     ? std::to_string(*args.max_events) +
                           " -- TRUNCATED RUN: --max-events was given, so events.input_tag0 below is "
                           "a PREFIX of the input, not the whole of it, and no yield or "
                           "normalisation from this file is complete"
                     : std::string("all (no --max-events; the whole input was read)"));
        prov.add("events.input_tag0", std::to_string(n_all));
        prov.add("events.written", std::to_string(n_written));
        // The parallelism used, and the guarantee that it did not touch the output.
        // Recorded so a file can be checked for the byte-identity claim after the
        // fact: same input + config + code_version at ANY --threads gives this same
        // slim and the same QA trees. n_threads == 1 is the farm's single-pass path.
        {
            const std::size_t recorded_parts =
                (n_threads == 1) ? 1 : pi0::partition_count(total_entries, args.partition_target);
            prov.add("compute.threads", std::to_string(n_threads));
            prov.add("compute.partitions", std::to_string(recorded_parts) +
                                               (n_threads == 1 ? " (single-pass path; --threads 1)"
                                                               : " (fixed by event count, not by thread count)"));
            prov.add("compute.output_determinism",
                     "byte-identical regardless of --threads: the entry range is split into a fixed number of "
                     "contiguous partitions decided by the event count, and each partition's rows are filled into "
                     "the trees single-threaded in partition-index (entry) order. The GBT scoring parallelises; the "
                     "emission order does not. Stage A has no floating-point reduction, so --threads 1 is EXACTLY "
                     "the pre-MT output, not merely equal-to-ulp.");
        }
        // Not cosmetic. These events had an unfilled RUN::config; the run branch
        // carries the file's run for them (see the fill site). Recorded so the
        // substitution is auditable and so a file where this count is large gets
        // looked at rather than trusted.
        prov.add("events.run_zero_in_bank",
                 std::to_string(n_run_zero) + " of " + std::to_string(n_all) +
                     " tag-0 events had RUN::config.run == 0 (bank not filled); the run branch carries the file's run");

        // In the SLIM's block, not only the QA file's. A reader holding the slim
        // must be able to tell that a per-candidate file was cut from this same
        // pass -- otherwise the QA file is an orphan whose numbers cannot be
        // traced back to the run that produced them, which is the hole this
        // project's provenance exists to close.
        prov.add("qa_ntuple.written",
                 args.qa_ntuple.has_value()
                     ? *args.qa_ntuple +
                           " -- DIAGNOSTIC ONLY: per-CANDIDATE pre-cut values, NOT per-event. No yield and no "
                           "normalisation may be taken from it. It changed nothing in this file."
                     : std::string("false (no --qa-ntuple)"));

        out->cd();
        prov.write(*out);
        prov.write_text(*out);
        tree->Write();
        out->Close();

        // ---- the QA ntuple's provenance --------------------------------------
        // THE SAME BLOCK, plus what is true only of this file.
        //
        // Copied rather than referenced: the QA file has to be readable on its
        // own. A diagnostic that records which cuts and which model produced it
        // can still be trusted a year later; one that says "see the slim" is one
        // nobody can trace once the two are in different directories. Same
        // argument as stageB_bin propagating Stage A's block forward rather than
        // dropping it.
        if (qa_out != nullptr) {
            // The row counts, checked rather than trusted. The -1s and the
            // independently-counted electron passes must agree: if they do not,
            // either electron_stage_index() mapped a rejection to -1 or the QA
            // fill drifted away from the verdict beside it -- and a QA file whose
            // "passed" rows are not the skim's passed rows is worse than no QA
            // file, because every plot drawn from it still looks right. Same
            // reasoning as the cutflow's own consistency check below.
            if (qa_e_tree->GetEntries() != n_has_electron) {
                throw std::runtime_error(
                    "the QA ntuple holds " + std::to_string(qa_e_tree->GetEntries()) + " electron rows but " +
                    std::to_string(n_has_electron) +
                    " trigger electrons were found. Every candidate must be recorded exactly once, before its "
                    "verdict is acted on.");
            }
            if (n_qa_electron_passed != n_electron_pass) {
                throw std::runtime_error(
                    "the QA ntuple marks " + std::to_string(n_qa_electron_passed) +
                    " electrons as failed_at == -1 but " + std::to_string(n_electron_pass) +
                    " passed the selection. failed_at does not agree with the verdict it was derived from.");
            }

            Provenance qa_prov = prov;
            qa_prov.add("qa_ntuple.path", *args.qa_ntuple);
            qa_prov.add("qa_ntuple.slim_path", args.output);
            // The loudest line in this block, and the reason the block is here. A
            // per-candidate file has an obvious number of rows and no obvious
            // reason not to divide by it.
            qa_prov.add("qa_ntuple.status",
                        "DIAGNOSTIC ONLY. These trees are PER-CANDIDATE, not per-event: qa_electron has one row "
                        "per trigger electron (pre-cut, rejected ones included) and qa_photon one row per pid==22 "
                        "candidate (pre-threshold). NO YIELD, EFFICIENCY OR NORMALISATION MAY BE TAKEN FROM THIS "
                        "FILE. The physics output of this pass is the slim at qa_ntuple.slim_path; this file "
                        "exists to plot the selection and what each cut removes.");
            qa_prov.add("qa_ntuple.qa_electron.rows",
                        std::to_string(qa_e_tree->GetEntries()) + " (one per trigger electron; failed_at == -1 on " +
                            std::to_string(n_qa_electron_passed) +
                            " of them, which is the cutflow's electron-pass count)");
            qa_prov.add("qa_ntuple.qa_photon.rows",
                        std::to_string(qa_g_tree->GetEntries()) +
                            " (one per pid==22 candidate IN EVENTS THAT ALREADY PASSED THE ELECTRON AND DIS CUTS "
                            "-- the photon loop is downstream of both, so this is not a photon count for the file)");
            qa_prov.add("qa_ntuple.gbt_score.absent",
                        "gbt_score is NaN wherever prefiltered == 1. The GBT pre-filter (energy, PCAL energy, "
                        "theta) runs BEFORE the classifier, so for those candidates the model was never evaluated "
                        "and no score exists. It is not zero: zero is a legal score meaning the classifier was "
                        "certain, and it was never asked.");
            qa_prov.add("qa_ntuple.failed_at",
                        "indexes the qa_electron_stages tree in this file, which is generated from "
                        "selection::kElectronStages -- the same array the cutflow's rows are built from. -1 means "
                        "the candidate passed all six cuts. The cuts short-circuit, so this is the FIRST failing "
                        "cut, not the only one: the tree is a SEQUENTIAL cutflow, exactly like the printed one.");

            qa_out->cd();
            qa_prov.write(*qa_out);
            qa_prov.write_text(*qa_out);
            qa_e_tree->Write();
            qa_g_tree->Write();
            qa_out->Close();
        }

        // ---- cutflow ---------------------------------------------------------
        std::ostringstream trigger_label;
        trigger_label << "has trigger electron (pid == " << pi0::selection::kPdgElectron << ", status < 0, |status| in ["
                      << cuts.electron.status_min << ", " << cuts.electron.status_max << "))";

        std::vector<CutflowRow> rows;
        // The first row's label is NOT the constant string it used to be. Under
        // --max-events "all events" is false, and this file's whole discipline
        // is that a cutflow label may not lie: the old analysis printed
        // "Momentum > 0.8 GeV" beside a p > 2.0 cut for years. A row that says
        // "all" while showing a truncated count is the same defect, and the
        // cutflow is the thing everyone quotes.
        std::ostringstream all_label;
        if (args.max_events.has_value()) {
            all_label << "events read (HIPO tag 0 = physics; TRUNCATED at --max-events "
                      << *args.max_events << ", NOT the whole file)";
        } else {
            all_label << "all events (HIPO tag 0 = physics)";
        }
        rows.push_back({all_label.str(), n_all});
        rows.push_back({trigger_label.str(), n_has_electron});

        long long surviving = n_has_electron;
        for (const char* stage : pi0::selection::kElectronStages) {
            surviving -= electron_fail[stage];
            rows.push_back({"electron: " + pi0::selection::electron_cutflow_label(stage, cuts), surviving});
        }

        // The electron rows are built by SUBTRACTING per-stage failures, but the
        // passes were counted independently. They must agree, and if they do not
        // then either a rejection was booked against a stage that is not in
        // kElectronStages or the two counters drifted -- in which case the table
        // is wrong and a wrong cutflow is worse than none, because it is the
        // thing everyone quotes. Cheap, so it is checked rather than assumed.
        if (surviving != n_electron_pass) {
            throw std::runtime_error("cutflow is inconsistent: subtracting the per-stage electron failures leaves " +
                                     std::to_string(surviving) + " but " + std::to_string(n_electron_pass) +
                                     " electrons were counted as passing. A stage id is missing from "
                                     "kElectronStages, or the counters have drifted.");
        }

        rows.push_back({"DIS: Q2 > " + fmt(cuts.dis.q2_min) + " GeV^2", n_q2});
        rows.push_back({"DIS: W > " + fmt(cuts.dis.w_min) + " GeV", n_w});
        rows.push_back({"DIS: y < " + fmt(cuts.dis.y_max), n_y});
        rows.push_back({">= 1 selected photon (GBT score > " + fmt(cuts.photon.gbt_threshold) + ")", n_has_photon});
        rows.push_back({"written to " + args.output, n_written});
        print_cutflow(rows);

        std::cout << "\nprovenance written to " << args.output << ":/provenance\n" << prov.as_text();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 4;
    }
}
