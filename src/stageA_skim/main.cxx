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
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TNamed.h>
#include <TObjString.h>
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
#include "photonid/Features.hpp"
#include "photonid/PhotonGBT.hpp"
#include "photonid/RunRangeModelMap.hpp"
#include "selection/ElectronSelection.hpp"
#include "selection/PhotonSelection.hpp"
#include "selection/SamplingFraction.hpp"
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
};

void print_usage(const char* argv0) {
    std::cerr
        << "usage: " << argv0 << " --input <file.hipo> --output <slim.root>\n"
        << "                        --config <cuts.json> --target <LD2|CxC|Cu|Sn>\n"
        << "                        [--run <N>] [--max-events <N>] [--qa-ntuple <qa.root>]\n"
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
        << "                slim. Rows are CANDIDATES, not events: take no yield from it.\n";
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

        // Counted independently of n_electron_pass so the two can be checked
        // against each other after the loop. See the check there.
        long long n_qa_electron_passed = 0;

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

        // ---- counters ------------------------------------------------------
        long long n_all = 0, n_has_electron = 0, n_electron_pass = 0;
        long long n_q2 = 0, n_w = 0, n_y = 0, n_has_photon = 0, n_written = 0;
        std::map<std::string, long long> electron_fail;  // stage id -> count

        // The single-run check the probe cannot do. The GBT model is chosen ONCE
        // per file, which is only meaningful if the file holds one run -- so that
        // assumption gets verified rather than trusted. Free here: the loop reads
        // RUN_config_run for every event anyway.
        std::map<int, long long> runs_seen;   // non-zero run -> events
        long long n_run_zero = 0;             // events whose RUN::config says 0

        // ---- the event loop -------------------------------------------------
        // --max-events is applied as an RDataFrame Range, i.e. the first N tag-0
        // entries -- a PREFIX of the file, not a sample of it. Fine for a smoke
        // test, which is what it is for; not a physics-representative subset,
        // which is why n_all, the cutflow's first row and the provenance all say
        // so rather than leaving the reader to notice.
        //
        // Range is single-thread only. This program never enables implicit MT,
        // so that costs nothing here -- but it is why this is a Range and not a
        // counter-and-break inside the lambda (which RDataFrame has no way to
        // honour) or an early return (which would still read the whole file).
        ROOT::RDF::RNode node = df;
        if (args.max_events.has_value()) node = node.Range(0u, *args.max_events);

        node.Foreach(
            [&](const RVecI& p_pid, const RVecF& p_px, const RVecF& p_py, const RVecF& p_pz, const RVecF& p_vz,
                const RVecF& p_beta, const RVecF& p_chi2pid, const RVecS& p_status,
                const RVecS& c_pindex, const RVecS& c_layer, const RVecS& c_sector, const RVecF& c_energy,
                const RVecF& c_x, const RVecF& c_y, const RVecF& c_z, const RVecF& c_m2u, const RVecF& c_m2v,
                const RVecF& c_lu, const RVecF& c_lv, const RVecF& c_lw,
                const RVecS& t_pindex, const RVecS& t_detector, const RVecS& t_layer, const RVecF& t_edge,
                const RVecS& ev_helicity, const RVecI& rc_run, const RVecI& rc_event) {
                ++n_all;

                // Single-row banks are RVecs too. Guard for empty, always.
                const int this_run = rc_run.empty() ? 0 : rc_run[0];
                const std::int64_t this_event = rc_event.empty() ? 0 : static_cast<std::int64_t>(rc_event[0]);
                if (this_run == 0) ++n_run_zero; else ++runs_seen[this_run];

                // HELICITY IS HWP-CORRECTED ALREADY -- nothing here needs to fix it.
                //
                // REC::Event carries TWO fields, and the CLAS12 bank definition
                // (data/bankdefs/hipo4/data.json) is explicit about the difference:
                //   helicity    -> "online-delay-corrected helicity, WITH HWP-correction (0=UDF)"
                //   helicityRaw -> "online-delay-corrected helicity (0=UDF)"
                // We read REC_Event_helicity, i.e. the CORRECTED one. The half-wave
                // plate state is per-run RCDB data applied UPSTREAM by the cooking
                // and already folded into this field, so an analysis reading
                // `helicity` (rather than `helicityRaw`) needs no HWP table of its
                // own. The superseded code did exactly the same thing, and its
                // "HWP-corrected" doc-comment was accurate, not stale.
                //
                // What IS true is that this makes the sign an INHERITED assumption:
                // it is only as good as the cooking's HWP bookkeeping, and no
                // HWP-in/HWP-out closure test exists anywhere in this chain. Worth
                // settling before publishing a SIGNED asymmetry -- but not a reason
                // to touch the value here.
                //
                // helicity == 0 means UNDEFINED; it is dropped downstream rather
                // than treated as a state.
                const int helicity = ev_helicity.empty() ? 0 : static_cast<int>(ev_helicity[0]);

                // ---- trigger electron -------------------------------------
                const auto e_row_opt = pi0::selection::find_trigger_electron(p_pid, p_status, p_px, p_py, p_pz, cuts);
                if (!e_row_opt.has_value()) return;
                ++n_has_electron;
                const std::size_t e = *e_row_opt;

                const pi0::photonid::CaloMap calo = pi0::photonid::CaloMap::build(
                    c_pindex, c_layer, c_sector, c_energy, c_x, c_y, c_z, c_m2u, c_m2v, c_lu, c_lv, c_lw);

                const double ex = p_px[e], ey = p_py[e], ez = p_pz[e];
                const double e_p = std::sqrt(ex * ex + ey * ey + ez * ez);
                const double e_theta_deg = theta_deg_of(ex, ey, ez);
                const double e_phi_deg = phi_deg_of(ex, ey);

                // Calorimeter quantities. A track with no ECAL row at all gets
                // zeros: sector 0 fails the sampling-fraction band (which returns
                // false, not throws, outside [1,6]) and lv/lw = 0 fails the PCAL
                // fiducial, so it is rejected either way.
                const pi0::photonid::CaloRowData* e_calo = calo.find(e);
                const double e_pcal_e = e_calo ? e_calo->pcal.e : 0.0;
                const double e_ecin_e = e_calo ? e_calo->ecin.e : 0.0;
                const double e_ecout_e = e_calo ? e_calo->ecout.e : 0.0;
                const int e_sector = e_calo ? e_calo->pcal.sector : 0;
                const double e_lv = e_calo ? e_calo->pcal.lv : 0.0;
                const double e_lw = e_calo ? e_calo->pcal.lw : 0.0;
                const double e_sf = e_p > 0.0 ? (e_pcal_e + e_ecin_e + e_ecout_e) / e_p : 0.0;

                // ---- vertex ------------------------------------------------
                // LD2 is tested on the RAW vz; the solid targets on the corrected
                // one. For Cu and Sn this is not merely a quality cut -- it is the
                // target assignment, since the two foils differ only in v_z.
                const double e_vz_raw = static_cast<double>(p_vz[e]);
                const double e_vz_used =
                    vz_corrected ? vz.correct(e_vz_raw, e_p, e_theta_deg, e_phi_deg, e_sector) : e_vz_raw;
                const bool vertex_passed = pi0::vertex::VzCorrector::pass_window(e_vz_used, vz_cuts);

                const double edge_r1 = traj_edge(t_pindex, t_detector, t_layer, t_edge, e, 6);
                const double edge_r2 = traj_edge(t_pindex, t_detector, t_layer, t_edge, e, 18);
                const double edge_r3 = traj_edge(t_pindex, t_detector, t_layer, t_edge, e, 36);

                const auto verdict = pi0::selection::pass_electron(static_cast<double>(p_chi2pid[e]), e_p, vertex_passed,
                                                                   e_sf, e_sector, e_lv, e_lw, edge_r1, edge_r2, edge_r3,
                                                                   cuts);

                // ---- QA row ------------------------------------------------
                // HERE, between the verdict and the return that acts on it: one
                // line later and the rejected candidates -- the entire point of
                // the file -- would be the ones missing from it.
                //
                // Every branch is a value pass_electron() was just handed, or the
                // verdict it returned. Nothing is recomputed for the plot, so a
                // figure drawn from this shows the cut that ran.
                if (qa_e_tree != nullptr) {
                    qe_run = run;
                    qe_event = this_event;
                    qe_sector = e_sector;
                    qe_p = e_p;
                    qe_theta_deg = e_theta_deg;
                    qe_phi_deg = e_phi_deg;
                    qe_vz = e_vz_raw;
                    // = vz when the target has no correction (LD2). The branch is
                    // still filled: a NaN here would make "corrected" mean two
                    // things, and the provenance already records which it was.
                    qe_vz_corrected = e_vz_used;
                    qe_chi2pid = static_cast<double>(p_chi2pid[e]);
                    qe_sf = e_sf;
                    qe_pcal_e = e_pcal_e;
                    qe_ecin_e = e_ecin_e;
                    qe_ecout_e = e_ecout_e;
                    qe_pcal_lv = e_lv;
                    qe_pcal_lw = e_lw;
                    qe_dc_edge_r1 = edge_r1;
                    qe_dc_edge_r2 = edge_r2;
                    qe_dc_edge_r3 = edge_r3;
                    // The SAME identifier the cutflow counts, mapped to the index
                    // of the row it belongs to. Not a second enumeration of the
                    // cuts -- see electron_stage_index().
                    qe_failed_at = pi0::selection::electron_stage_index(verdict.failed_at);
                    if (qe_failed_at < 0) ++n_qa_electron_passed;
                    qa_e_tree->Fill();
                }

                if (!verdict.passed) {
                    ++electron_fail[verdict.failed_at];
                    return;
                }
                ++n_electron_pass;

                // ---- DIS ---------------------------------------------------
                // E' = sqrt(p^2 + m_e^2). The electron mass is carried, not
                // dropped -- which is what the old code did too. The difference
                // from the massless E' = |p| is m_e^2/(2p) ~ 2.6e-8 GeV at
                // p = 5 GeV, i.e. a relative ~5e-9 on Q2: far below the
                // resolution these cuts are drawn against, and it changes no
                // event's verdict at q2_min = 1. It is here because it is free
                // and because "E'" should mean the electron's energy rather
                // than an approximation to it. tests/test_core.cpp pins the
                // size of the shift so the change is auditable rather than
                // invisible.
                //
                // NOTE what this does NOT make exact: compute_dis still forms
                // Q2 by the angle form 4*E*E'*sin^2(theta/2), which is itself
                // the massless-electron expression (cuts.json's /dis block says
                // so). Carrying m_e in E' does not turn it into -(k-k')^2.
                // Photons stay massless by convention: E_gamma = |p_gamma|.
                const double ee = pi0::energy_from_p(e_p, pi0::kElectronMassGeV);
                const pi0::DisKin dis = pi0::kin::compute_dis(ex, ey, ez, ee, cuts.beam.energy_gev);

                if (!(dis.q2 > cuts.dis.q2_min)) return;
                ++n_q2;
                if (!(dis.w > cuts.dis.w_min)) return;
                ++n_w;
                if (!(dis.y < cuts.dis.y_max)) return;
                ++n_y;

                // ---- photons ------------------------------------------------
                b_gpx.clear();
                b_gpy.clear();
                b_gpz.clear();
                b_g_e_gamma_deg.clear();

                const std::size_t n_rows = std::min({p_pid.size(), p_px.size(), p_py.size(), p_pz.size(), p_beta.size()});
                for (std::size_t r = 0; r < n_rows; ++r) {
                    if (p_pid[r] != pi0::selection::kPdgPhoton) continue;

                    const pi0::photonid::CaloRowData* g_calo = calo.find(r);
                    // No ECAL row -> PCAL energy 0 -> the pre-filter's E_PCAL > 0
                    // rejects it, and score_fn is never invoked. That matters:
                    // build_features() throws for a row with no calorimeter data,
                    // and this is what keeps it from ever seeing one.
                    const double g_pcal_e = g_calo ? g_calo->pcal.e : 0.0;
                    const double g_lv = g_calo ? g_calo->pcal.lv : 0.0;
                    const double g_lw = g_calo ? g_calo->pcal.lw : 0.0;

                    const double gx = p_px[r], gy = p_py[r], gz = p_pz[r];

                    // MOST CLUSTERS HAVE NO SCORE, AND THAT IS NOT A GAP TO FILL.
                    // The pre-filter runs before the callable, so for anything it
                    // rejects the GBT is never evaluated and no score exists. The
                    // QA row says so with prefiltered == 1 and a NaN, rather than
                    // a 0 -- which is a legal score, and would read as "the
                    // classifier was certain this was not a photon" instead of
                    // "the classifier was never asked". Forcing the GBT to run on
                    // everything just to fill the column would change what this
                    // program costs in order to describe what it does.
                    double g_score = std::numeric_limits<double>::quiet_NaN();
                    bool g_scored = false;

                    // The scored overload, not the eager one: the GBT is the
                    // expensive part of this program, and the pre-filter exists
                    // precisely so it is not evaluated on most clusters.
                    const bool ok = pi0::selection::pass_photon_scored(
                        p_pid[r], gx, gy, gz, g_pcal_e, g_lv, g_lw, static_cast<double>(p_beta[r]),
                        [&]() {
                            const std::vector<float> feats =
                                pi0::photonid::build_features(r, calo, p_pid, p_px, p_py, p_pz, feat_cuts);
                            // Captured, not returned twice: this is the number the
                            // threshold is about to be applied to, so recording it
                            // here is recording the cut rather than modelling it.
                            g_score = pi0::photonid::score(feats, model);
                            g_scored = true;
                            return g_score;
                        },
                        cuts);

                    // ---- QA row ---------------------------------------------
                    // Before the `continue`, so the rejected candidates are here.
                    if (qa_g_tree != nullptr) {
                        qg_run = run;
                        qg_event = this_event;
                        // The selection's OWN energy and angle, not a second sqrt
                        // and a second acos that agree today.
                        qg_e = pi0::selection::photon_energy_gev(gx, gy, gz);
                        qg_theta_deg = pi0::selection::photon_theta_deg(gx, gy, gz);
                        // phi is the one branch here no cut reads: the photon
                        // selection has no azimuthal term. It is carried because a
                        // fiducial map is drawn in (theta, phi) and this file
                        // exists to be plotted.
                        qg_phi_deg = phi_deg_of(gx, gy);
                        qg_gbt_score = g_score;
                        qg_pcal_e = g_pcal_e;
                        qg_pcal_lv = g_lv;
                        qg_pcal_lw = g_lw;
                        qg_beta = static_cast<double>(p_beta[r]);
                        qg_prefiltered = g_scored ? 0 : 1;
                        qg_passed = ok ? 1 : 0;
                        qa_g_tree->Fill();
                    }

                    if (!ok) continue;

                    b_gpx.push_back(gx);
                    b_gpy.push_back(gy);
                    b_gpz.push_back(gz);
                    // Precomputed HERE so that the pairing stage never needs the
                    // electron again: it is the only thing the e-gamma angle cut
                    // needs from it, and a slim file that carried the electron
                    // just for this would invite recomputing it inconsistently.
                    b_g_e_gamma_deg.push_back(pi0::kin::angle_between_deg(gx, gy, gz, ex, ey, ez));
                }

                // >= 1, NOT >= 2 -- single-photon events feed the mixed-event
                // donor pool downstream. See the header comment.
                if (b_gpx.empty()) return;
                ++n_has_photon;

                // The FILE's run, not this_run. They differ only when RUN::config
                // was not filled for this event, in which case this_run is 0 --
                // and a 0 in the run branch is a landmine for every downstream
                // per-run normalisation, QA plot and file_hash mixing seed. The
                // single-run check after the loop is what licenses this
                // substitution: once the file is known to hold exactly one run,
                // the file's run IS this event's run. The count of events whose
                // bank said 0 goes into the provenance, so the substitution is
                // recorded rather than hidden.
                b_run = run;
                b_event = this_event;
                b_helicity = helicity;
                b_q2 = dis.q2;
                b_xb = dis.xb;
                b_nu = dis.nu;
                b_w = dis.w;
                b_y = dis.y;
                b_ex = ex;
                b_ey = ey;
                b_ez = ez;
                b_ee = ee;
                tree->Fill();
                ++n_written;
            },
            needed);

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
