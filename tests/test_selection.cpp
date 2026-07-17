// Unit tests for the cut configuration and the electron/photon selection.
//
// These tests exist to make four specific historical failures impossible. All
// four shipped in the analysis this code replaces:
//
//   1. A CUTFLOW LABEL THAT LIED. The old cutflow printed "Momentum > 0.8 GeV"
//      while applying p > 2.0. The string was a literal that had drifted from
//      the constant beside it, and it was copied into a second algorithm and
//      into every log anyone ever quoted from.
//   2. CUT VALUES SCATTERED ACROSS FILES. Thresholds lived in C++ constants,
//      in six TOML files, and in dead config keys that were declared, set, and
//      never read. Three keys, two values, one effective cut.
//   3. SILENT DEFAULTS. A missing key produced a plausible number rather than
//      an error, so a config that did not say what you thought it said ran
//      anyway.
//   4. DEGREE/RADIAN CONFUSION AT A PUBLIC BOUNDARY, twice, independently. The
//      old photon pre-filter took `double theta` in RADIANS and converted
//      inside -- an unmarked radian at an API boundary.
//
// Method note: reference values here are stated as literals worked out from the
// SOURCES OF TRUTH (config/cuts.json for cuts, SamplingFractionService.hpp for
// the calibration), and never by calling the function under test. Where a
// polynomial is checked, its coefficients are re-transcribed from the source
// file INTO THIS TEST and the arithmetic is done locally -- so a transcription
// typo in SamplingFraction.cpp has to be made twice, identically, to pass.
//
// The cuts.json values asserted below are transcribed from that file, which is
// the single source of truth. If they ever disagree with it, IT WINS and this
// file is what changes -- but a disagreement means someone moved a cut, and
// that should be a deliberate act with a physics reason, which is exactly why
// these assertions are here to force the conversation.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <nlohmann/json.hpp>

#include <ROOT/RVec.hxx>

#include "config/Cuts.hpp"
#include "core/Constants.hpp"
#include "selection/ElectronSelection.hpp"
#include "selection/PhotonSelection.hpp"
#include "selection/SamplingFraction.hpp"
#include "vertex/VzCorrector.hpp"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using ROOT::VecOps::RVec;

namespace sel = pi0::selection;

namespace {

// PI0_CUTS_JSON is the absolute path to config/cuts.json, injected by
// tests/meson.build. The tests read the REAL file: a test that parses a
// hand-written copy of the config proves only that the copy parses.
const char* cuts_json_path() { return PI0_CUTS_JSON; }

const pi0::Cuts& shipped_cuts() {
    static const pi0::Cuts c = pi0::Cuts::load(cuts_json_path());
    return c;
}

/// Write `doc` to a scratch file and return the path. Used to build a config
/// that is the real one MINUS one key, so the loader is tested against a
/// document that is realistic in every respect except the defect under test.
std::string write_temp_json(const nlohmann::json& doc, const std::string& tag) {
    const std::string path = std::string(PI0_TEST_TMPDIR) + "/cuts_" + tag + ".json";
    std::ofstream out(path);
    REQUIRE(out.good());
    out << doc.dump(2);
    out.close();
    return path;
}

nlohmann::json parsed_cuts_json() {
    std::ifstream in(cuts_json_path());
    REQUIRE(in.good());
    nlohmann::json doc;
    in >> doc;
    return doc;
}

/// Remove a dotted key path from `doc`. Fails the test if it was not there --
/// otherwise a renamed key would make the "throws on missing key" test pass
/// vacuously, by deleting nothing.
void erase_path(nlohmann::json& doc, const std::string& dotted) {
    nlohmann::json* cur = &doc;
    std::string part;
    std::size_t start = 0;
    for (;;) {
        const std::size_t dot = dotted.find('.', start);
        part = dotted.substr(start, dot == std::string::npos ? std::string::npos : dot - start);
        if (dot == std::string::npos) break;
        REQUIRE(cur->contains(part));
        cur = &(*cur)[part];
        start = dot + 1;
    }
    REQUIRE(cur->contains(part));
    cur->erase(part);
}

/// A set of electron arguments that passes every cut, so that each test can
/// spoil exactly one of them and watch that one stage report itself.
struct GoodElectron {
    double chi2pid = 0.0;
    double p_gev = 4.0;
    bool vertex_passed = true;
    double sampling_fraction = 0.0;  // filled in by make_good_electron()
    int pcal_sector = 3;
    double pcal_lv_cm = 20.0;
    double pcal_lw_cm = 20.0;
    double dc_edge_r1_cm = 5.0;
    double dc_edge_r2_cm = 5.0;
    double dc_edge_r3_cm = 15.0;

    [[nodiscard]] sel::ElectronCutResult apply(const pi0::Cuts& cuts) const {
        return sel::pass_electron(chi2pid, p_gev, vertex_passed, sampling_fraction, pcal_sector,
                                  pcal_lv_cm, pcal_lw_cm, dc_edge_r1_cm, dc_edge_r2_cm,
                                  dc_edge_r3_cm, cuts);
    }
};

/// The polarity the SHIPPED config selects. The electron helpers below use this
/// rather than naming Outbending, so that they keep testing pass_electron's
/// actual behaviour if cuts.json ever switches polarity. (The sampling-fraction
/// tests further up do name a polarity explicitly -- those are about the tables
/// themselves, not about what the config selects.)
sel::Polarity shipped_polarity() {
    return sel::polarity_from_string(shipped_cuts().beam.polarity);
}

GoodElectron make_good_electron(const pi0::Cuts& cuts) {
    GoodElectron e;
    // Sit exactly on the band centre, so the SF stage passes for any n_sigma > 0.
    e.sampling_fraction =
        sel::sf_mu(e.p_gev, e.pcal_sector, sel::polarity_from_string(cuts.beam.polarity));
    return e;
}

}  // namespace

// ===========================================================================
// Cuts::load -- the values
// ===========================================================================

TEST_CASE("Cuts::load parses config/cuts.json and yields the documented values", "[cuts]") {
    const pi0::Cuts& c = shipped_cuts();

    SECTION("beam") {
        CHECK_THAT(c.beam.energy_gev, WithinAbs(10.53, 1e-12));
        CHECK(c.beam.polarity == "outbending");
    }

    SECTION("the momentum cut is 2.0 GeV, and has always been 2.0") {
        // THE HEADLINE TRAP. The old cutflow's label said 0.8 while the code
        // applied 2.0 (ElectronCutsService::MIN_ELECTRON_MOMENTUM). Any cutflow
        // table copied from the old logs inherits the lie. 0.8 is not a value
        // this analysis has ever used for anything.
        CHECK_THAT(c.electron.min_momentum_gev, WithinAbs(2.0, 1e-12));
        CHECK(c.electron.min_momentum_gev != 0.8);
    }

    SECTION("electron") {
        CHECK_THAT(c.electron.chi2pid_min, WithinAbs(-5.0, 1e-12));
        CHECK_THAT(c.electron.chi2pid_max, WithinAbs(5.0, 1e-12));
        CHECK_THAT(c.electron.sf_n_sigma, WithinAbs(3.5, 1e-12));
        CHECK_THAT(c.electron.dc_edge_r1_cm, WithinAbs(1.68, 1e-12));
        CHECK_THAT(c.electron.dc_edge_r2_cm, WithinAbs(2.0, 1e-12));
        CHECK_THAT(c.electron.dc_edge_r3_cm, WithinAbs(8.75, 1e-12));
        // |status| in [2000, 4000): inclusive lower, EXCLUSIVE upper.
        CHECK(c.electron.status_min == 2000);
        CHECK(c.electron.status_max == 4000);
    }

    SECTION("dis") {
        CHECK_THAT(c.dis.q2_min, WithinAbs(1.0, 1e-12));
        CHECK_THAT(c.dis.w_min, WithinAbs(2.0, 1e-12));
        CHECK_THAT(c.dis.y_max, WithinAbs(0.85, 1e-12));
    }

    SECTION("photon") {
        CHECK_THAT(c.photon.gbt_threshold, WithinAbs(0.78, 1e-12));
        CHECK_THAT(c.photon.min_energy_gev, WithinAbs(0.2, 1e-12));
        CHECK_THAT(c.photon.theta_min_deg, WithinAbs(5.0, 1e-12));
        CHECK_THAT(c.photon.theta_max_deg, WithinAbs(35.0, 1e-12));
        CHECK_THAT(c.photon.beta_min, WithinAbs(0.9, 1e-12));
        CHECK_THAT(c.photon.beta_max, WithinAbs(1.1, 1e-12));
        CHECK(c.photon.gbt_pass == 1);
    }

    SECTION("allow_rga_fallback is TRUE -- a decision, not a default") {
        // CHANGED 2026-07-17. This was false, and this assertion read
        // `== false`. No GBT model exists for RG-D -- the lookup table tops out
        // at run 16772 and RG-D runs start at 18305 -- so with it false
        // stageA_skim exits 3 on every RG-D run and there is no skim at all.
        // Running the production at all means taking the RG-A INBENDING PASS-1
        // model, and taking it makes the production A STUDY, NOT A MEASUREMENT.
        //
        // This assertion is a TRIPWIRE IN BOTH DIRECTIONS, which is why it is
        // updated rather than deleted: flip the key and this test fails, and
        // whoever flipped it reads this. Setting it back to false is the right
        // move the day an RG-D-trained model exists -- PROVENANCE.md calls that
        // the single highest-value fix available to the photon selection.
        //
        // WHAT THIS DOES NOT WEAKEN: the code's refusal is tested against a
        // literal, not against the shipped config --
        //   test_photonid.cpp, "select_model REFUSES an RG-D run when the
        //   fallback is not allowed"
        // still requires select_model(18500, 1, false) to throw. The safety
        // property lives in the code; this key is only the policy, and every
        // output made under it carries gbt.fallback_used = TRUE.
        CHECK(c.photon.allow_rga_fallback == true);
    }

    SECTION("photon PCAL fiducial is 14.0 and the electron's is 9.0 -- NOT unified") {
        // cuts.json says the difference is intentional, in both blocks. This
        // assertion exists so that a future "cleanup" that unifies them fails
        // here rather than silently changing the photon acceptance.
        CHECK_THAT(c.photon.pcal_lv_min_cm, WithinAbs(14.0, 1e-12));
        CHECK_THAT(c.photon.pcal_lw_min_cm, WithinAbs(14.0, 1e-12));
        CHECK_THAT(c.electron.pcal_lv_min_cm, WithinAbs(9.0, 1e-12));
        CHECK_THAT(c.electron.pcal_lw_min_cm, WithinAbs(9.0, 1e-12));
        CHECK(c.photon.pcal_lv_min_cm > c.electron.pcal_lv_min_cm);
    }

    SECTION("pairing") {
        CHECK_THAT(c.pairing.mass_window_gev, WithinAbs(0.2, 1e-12));
        CHECK_THAT(c.pairing.min_mgg_gev, WithinAbs(0.001, 1e-12));
        CHECK_THAT(c.pairing.e_gamma_min_angle_deg, WithinAbs(8.0, 1e-12));
        CHECK_THAT(c.pairing.open_a_deg, WithinAbs(17.561, 1e-12));
        CHECK_THAT(c.pairing.open_b_inv_gev, WithinAbs(0.756, 1e-12));
        CHECK_THAT(c.pairing.open_offset_deg, WithinAbs(1.0, 1e-12));
    }

    SECTION("vertex -- LD2 is a raw flat window, uncorrected") {
        // The cryotarget is one long cell: no foil peaks to separate, so no
        // correction and no peak rule. The bounds are tab:vz verbatim.
        CHECK(c.vertex.ld2.correction_enabled == false);
        CHECK(c.vertex.ld2.rule == pi0::vertex::VzRule::RawWindow);
        CHECK_THAT(c.vertex.ld2.vz_min_cm, WithinAbs(-15.0, 1e-12));
        CHECK_THAT(c.vertex.ld2.vz_max_cm, WithinAbs(5.0, 1e-12));
    }

    SECTION("vertex -- Cu is the -7.9 foil and Sn is the -2.9 foil") {
        // ***** THE NAMING TRAP, PINNED AGAINST THE SHIPPED FILE *****
        //
        // The fit that produced these called its two parameter columns 'lo' and
        // 'hi', and 'hi' is the MORE NEGATIVE (upstream) foil -- the old
        // service's doc-comment said the opposite. cuts.json drops lo/hi and
        // names targets, so the trap is unreachable from here; what remains is
        // to pin WHICH NUMBER each target carries, because for Cu and Sn the
        // vertex window is not a quality cut, it IS the target assignment.
        // Swapping them still compiles, still runs, and silently analyses the
        // wrong nucleus.
        //
        // tests/test_vertex.cpp pins the acceptance semantics against the same
        // numbers, transcribed; this pins the file itself.
        REQUIRE(c.vertex.cu.peaks.size() == 1u);
        CHECK_THAT(c.vertex.cu.peaks[0].mu_cm, WithinAbs(-7.861, 1e-12));
        CHECK_THAT(c.vertex.cu.peaks[0].sigma_cm, WithinAbs(0.415, 1e-12));

        REQUIRE(c.vertex.sn.peaks.size() == 1u);
        CHECK_THAT(c.vertex.sn.peaks[0].mu_cm, WithinAbs(-2.916, 1e-12));
        CHECK_THAT(c.vertex.sn.peaks[0].sigma_cm, WithinAbs(0.370, 1e-12));

        // Cu is upstream of Sn. The one-line statement of the whole trap.
        CHECK(c.vertex.cu.peaks[0].mu_cm < c.vertex.sn.peaks[0].mu_cm);

        // Both are corrected, and each selects exactly ONE peak -- Cu does not
        // accept Sn's foil, which is what "one peak each" buys.
        CHECK(c.vertex.cu.correction_enabled);
        CHECK(c.vertex.sn.correction_enabled);
        CHECK(c.vertex.cu.rule == pi0::vertex::VzRule::CorrectedPeaks);
        CHECK(c.vertex.sn.rule == pi0::vertex::VzRule::CorrectedPeaks);
        CHECK_THAT(c.vertex.cu.n_sigma, WithinAbs(3.0, 1e-12));
        CHECK_THAT(c.vertex.sn.n_sigma, WithinAbs(3.0, 1e-12));
    }

    SECTION("vertex -- the resulting Cu/Sn windows are the note's tab:vz") {
        // cuts.json states the arithmetic it expects: "Resulting window -9.106
        // to -6.616 cm (= -7.861 -+ 3*0.415), which matches the note's tab:vz --
        // if your code produces a different window from these inputs, the code
        // is wrong." So check the code against the note, not just the inputs
        // against themselves.
        const auto& cu = c.vertex.cu.peaks[0];
        CHECK_THAT(cu.mu_cm - c.vertex.cu.n_sigma * cu.sigma_cm, WithinAbs(-9.106, 1e-9));
        CHECK_THAT(cu.mu_cm + c.vertex.cu.n_sigma * cu.sigma_cm, WithinAbs(-6.616, 1e-9));

        const auto& sn = c.vertex.sn.peaks[0];
        CHECK_THAT(sn.mu_cm - c.vertex.sn.n_sigma * sn.sigma_cm, WithinAbs(-4.026, 1e-9));
        CHECK_THAT(sn.mu_cm + c.vertex.sn.n_sigma * sn.sigma_cm, WithinAbs(-1.806, 1e-9));

        // The two windows do not overlap -- if they did, the Cu/Sn assignment
        // would be ambiguous for some events rather than exclusive.
        CHECK(cu.mu_cm + c.vertex.cu.n_sigma * cu.sigma_cm <
              sn.mu_cm - c.vertex.sn.n_sigma * sn.sigma_cm);
    }

    SECTION("vertex -- CxC takes EITHER foil, with its own peak parameters") {
        // Two carbon foils in the same two z slots, so either peak is signal.
        REQUIRE(c.vertex.cxc.peaks.size() == 2u);
        CHECK(c.vertex.cxc.correction_enabled);
        CHECK(c.vertex.cxc.rule == pi0::vertex::VzRule::CorrectedPeaks);
        CHECK_THAT(c.vertex.cxc.peaks[0].mu_cm, WithinAbs(-7.887, 1e-12));
        CHECK_THAT(c.vertex.cxc.peaks[0].sigma_cm, WithinAbs(0.395, 1e-12));
        CHECK_THAT(c.vertex.cxc.peaks[1].mu_cm, WithinAbs(-2.906, 1e-12));
        CHECK_THAT(c.vertex.cxc.peaks[1].sigma_cm, WithinAbs(0.373, 1e-12));

        // NOT the CuSn numbers reused: cuts.json says they are separate
        // measurements and "must not be deduplicated". A well-meaning tidy-up
        // that shares one peak table between CuSn and CxC fails here.
        CHECK(c.vertex.cxc.peaks[0].mu_cm != c.vertex.cu.peaks[0].mu_cm);
        CHECK(c.vertex.cxc.peaks[1].mu_cm != c.vertex.sn.peaks[0].mu_cm);
    }

    SECTION("vertex -- the sigma that OVERRIDES the params file is the one here") {
        // The correction polynomials come from data/Vz/vz_corrector_params.txt;
        // the cut window comes from this file, and for cusn_outbending's
        // upstream peak they disagree: the params file records sigma = 0.3851,
        // this says 0.415 (+8%). 0.415 is what produced the published
        // selection. Open question Q2 in the note -- recorded, not resolved.
        // tests/test_vertex.cpp pins the params file's side of the same fact.
        CHECK_THAT(c.vertex.cu.peaks[0].sigma_cm, WithinAbs(0.415, 1e-12));
        CHECK(std::fabs(c.vertex.cu.peaks[0].sigma_cm - 0.3851194323) > 0.02);
    }

    SECTION("the mass window is wide, and is NOT a +-30 MeV pi0 selection") {
        // Documentation claiming +-30 MeV is wrong: the old C++ default was
        // 0.03 but every shipped config overrode it to 0.2, so the docstring
        // described a cut that never ran. |m_gg - 0.135| < 0.2 reduces to
        // m_gg < 0.335 with NO lower bound, which is by design -- the sidebands
        // are what the background fit needs.
        CHECK(c.pairing.mass_window_gev > 0.1);
    }
}

// ===========================================================================
// Cuts::load -- failing loudly
// ===========================================================================

TEST_CASE("Cuts::load THROWS on a missing key, and names it", "[cuts][loud]") {
    // Every key below is deleted from an OTHERWISE COMPLETE copy of the real
    // config, so each case isolates exactly one missing key.
    const auto key = GENERATE(std::string("beam.energy_gev"),
                              std::string("beam.polarity"),
                              std::string("electron.min_momentum_gev"),
                              std::string("electron.chi2pid_min"),
                              std::string("electron.sampling_fraction.n_sigma"),
                              std::string("electron.trigger.status_abs_max"),
                              std::string("dis.q2_min"),
                              std::string("photon.gbt_threshold"),
                              std::string("photon.allow_rga_fallback"),
                              std::string("photon.beta_max"),
                              std::string("pairing.mass_window_gev"),
                              std::string("pairing.opening_angle.a_deg"),
                              // The /vertex block. Until these were read, this
                              // block sat in cuts.json unread while a hard-coded
                              // table in vertex/VzCorrector.cpp supplied the
                              // windows -- the "declared, set, never read"
                              // defect this file exists to end. A missing key
                              // here must be as loud as any other.
                              std::string("vertex.targets.ld2.vz_min_cm"),
                              std::string("vertex.targets.ld2.rule"),
                              std::string("vertex.targets.ld2.correction_enabled"),
                              std::string("vertex.targets.cu.n_sigma"),
                              std::string("vertex.targets.cu.peaks"),
                              std::string("vertex.targets.sn.peaks"),
                              std::string("vertex.targets.cxc.peaks"),
                              std::string("vertex.targets.cxc.correction_enabled"));

    INFO("deleted key: " << key);
    nlohmann::json doc = parsed_cuts_json();
    erase_path(doc, key);
    const std::string path = write_temp_json(doc, "missing");

    // Not merely "throws": the message must name the key, or the failure is a
    // puzzle rather than a report.
    REQUIRE_THROWS_WITH(pi0::Cuts::load(path), ContainsSubstring(key));
    std::remove(path.c_str());
}

TEST_CASE("Cuts::load THROWS rather than defaulting", "[cuts][loud]") {
    SECTION("a nonexistent file") {
        REQUIRE_THROWS_AS(pi0::Cuts::load("/nonexistent/definitely/not/cuts.json"),
                          std::runtime_error);
    }

    SECTION("a file that is not JSON") {
        const std::string path = std::string(PI0_TEST_TMPDIR) + "/not_json.json";
        {
            std::ofstream out(path);
            out << "this is not json";
        }
        REQUIRE_THROWS_AS(pi0::Cuts::load(path), std::runtime_error);
        std::remove(path.c_str());
    }

    SECTION("an empty object -- i.e. every key missing at once") {
        const std::string path = write_temp_json(nlohmann::json::object(), "empty");
        REQUIRE_THROWS_AS(pi0::Cuts::load(path), std::runtime_error);
        std::remove(path.c_str());
    }

    SECTION("a key of the wrong type") {
        nlohmann::json doc = parsed_cuts_json();
        doc["electron"]["min_momentum_gev"] = "2.0";  // a string, not a number
        const std::string path = write_temp_json(doc, "wrongtype");
        REQUIRE_THROWS_WITH(pi0::Cuts::load(path), ContainsSubstring("min_momentum_gev"));
        std::remove(path.c_str());
    }

    SECTION("an unknown polarity") {
        nlohmann::json doc = parsed_cuts_json();
        doc["beam"]["polarity"] = "sideways";
        const std::string path = write_temp_json(doc, "polarity");
        REQUIRE_THROWS_WITH(pi0::Cuts::load(path), ContainsSubstring("polarity"));
        std::remove(path.c_str());
    }
}

TEST_CASE("Cuts::load refuses an inbending config", "[cuts][loud][vertex]") {
    // NOT a limitation that snuck in -- a refusal that replaces a silent wrong
    // answer. cuts.json's /vertex block is the OUTBENDING set, as its beam block
    // says of every polarity-dependent number in the file. The windows used to
    // come from a hard-coded table that ALSO carried inbending values the file
    // never claimed to have; that table is deleted, so there is now nothing
    // behind an inbending request but the outbending numbers.
    //
    // Loading them for inbending data would be exactly the silent
    // wrong-number failure that deleting the table was meant to prevent, so the
    // loader refuses. Adding inbending means adding a per-polarity /vertex
    // block, which is what cuts.json's own _polarity_comment prescribes.
    nlohmann::json doc = parsed_cuts_json();
    doc["beam"]["polarity"] = "inbending";
    const std::string path = write_temp_json(doc, "inbending");

    // The message must say WHY, not merely that something is wrong.
    REQUIRE_THROWS_WITH(pi0::Cuts::load(path), ContainsSubstring("OUTBENDING"));
    std::remove(path.c_str());

    // ...and the shipped file is outbending, so this refusal is not firing on
    // the real config.
    CHECK(shipped_cuts().beam.polarity == "outbending");
}

TEST_CASE("Cuts::load refuses a vertex block that disagrees with itself", "[cuts][loud][vertex]") {
    SECTION("an unknown rule") {
        nlohmann::json doc = parsed_cuts_json();
        doc["vertex"]["targets"]["cu"]["rule"] = "voigt_window";
        const std::string path = write_temp_json(doc, "vzrule");
        REQUIRE_THROWS_WITH(pi0::Cuts::load(path), ContainsSubstring("vertex.targets.cu.rule"));
        std::remove(path.c_str());
    }

    SECTION("correction_enabled contradicting the rule") {
        // 'corrected_peaks' peaks are fitted on the CORRECTED v_z; a raw_window
        // bounds the RAW one. A config that says corrected_peaks with the
        // correction off would have the skim cut a raw vertex against corrected
        // peaks -- a window silently drawn around the wrong quantity.
        nlohmann::json doc = parsed_cuts_json();
        doc["vertex"]["targets"]["cu"]["correction_enabled"] = false;
        const std::string path = write_temp_json(doc, "vzcorr");
        REQUIRE_THROWS_WITH(pi0::Cuts::load(path),
                            ContainsSubstring("vertex.targets.cu.correction_enabled"));
        std::remove(path.c_str());
    }

    SECTION("an empty peak list") {
        // Accepts nothing at all -- a selection nobody meant to write.
        nlohmann::json doc = parsed_cuts_json();
        doc["vertex"]["targets"]["sn"]["peaks"] = nlohmann::json::array();
        const std::string path = write_temp_json(doc, "vzpeaks");
        REQUIRE_THROWS_WITH(pi0::Cuts::load(path), ContainsSubstring("vertex.targets.sn.peaks"));
        std::remove(path.c_str());
    }

    SECTION("a peak missing its sigma, named by its index") {
        // The message must locate the defect inside the array, not just say
        // "peaks is wrong".
        nlohmann::json doc = parsed_cuts_json();
        doc["vertex"]["targets"]["cxc"]["peaks"][1].erase("sigma_cm");
        const std::string path = write_temp_json(doc, "vzpeaksig");
        REQUIRE_THROWS_WITH(pi0::Cuts::load(path),
                            ContainsSubstring("vertex.targets.cxc.peaks[1].sigma_cm"));
        std::remove(path.c_str());
    }
}

TEST_CASE("Cuts::load refuses a config that contradicts core/Constants.hpp", "[cuts][loud]") {
    // Constants.hpp and cuts.json each spell out the beam energy, the proton
    // mass and the pi0 mass. Constants.hpp's own header warns that a duplicated
    // constant which drifts out of sync is "the classic way an analysis
    // silently changes physics", and cuts.json exists partly because the old
    // tree carried THREE beam energies that agreed only by coincidence.
    //
    // The duplication is nonetheless still real -- so the loader pins them
    // together. These tests are what keeps that pin honest.
    const auto [key_a, key_b, bad_value] =
        GENERATE(std::make_tuple("beam", "energy_gev", 10.6),
                 std::make_tuple("beam", "proton_mass_gev", 0.9383),
                 std::make_tuple("pairing", "pi0_mass_gev", 0.135));

    INFO("drifted key: " << key_a << "." << key_b);
    nlohmann::json doc = parsed_cuts_json();
    doc[key_a][key_b] = bad_value;
    const std::string path = write_temp_json(doc, "drift");

    REQUIRE_THROWS_WITH(pi0::Cuts::load(path), ContainsSubstring("Constants.hpp"));
    std::remove(path.c_str());
}

TEST_CASE("Cuts::load refuses a config whose pid keys contradict the code", "[cuts][loud]") {
    // The Cuts struct carries no pid field: a PDG code is an identity, not a
    // threshold. cuts.json still records both, so the loader verifies them --
    // otherwise the file could quietly describe a selection the code does not
    // implement, which is what the old tree's dead TOML keys were.
    nlohmann::json doc = parsed_cuts_json();
    doc["photon"]["pid"] = 2112;  // a neutron, say
    const std::string path = write_temp_json(doc, "pid");
    REQUIRE_THROWS_WITH(pi0::Cuts::load(path), ContainsSubstring("photon.pid"));
    std::remove(path.c_str());
}

// ===========================================================================
// Sampling fraction
// ===========================================================================

TEST_CASE("sampling fraction mu/sigma reproduce a hand-evaluated polynomial", "[sf]") {
    // Coefficients RE-TRANSCRIBED here, independently, from the source of truth
    //   clas-analysis-1/clas12/include/clas12/services/SamplingFractionService.hpp
    //   -> initialize_data() -> m_data["OB"], third SectorSFParameters entry.
    // Typed in again on purpose: if SamplingFraction.cpp has a transcription
    // typo, this test only passes if the identical typo was made twice.
    SECTION("outbending sector 3") {
        constexpr double mu_a = 2.258053E-01, mu_b = 9.073840E-03, mu_c = -1.240745E-03,
                         mu_d = 3.579321E-05;
        constexpr double sg_a = 2.399538E-02, sg_b = -4.208668E-03, sg_c = 4.503779E-04,
                         sg_d = -1.783519E-05;

        const auto p = GENERATE(0.0, 1.0, 2.5, 4.0, 7.25, 10.53);
        INFO("p = " << p);

        const double mu_expected = mu_a + mu_b * p + mu_c * p * p + mu_d * p * p * p;
        const double sg_expected = sg_a + sg_b * p + sg_c * p * p + sg_d * p * p * p;

        CHECK_THAT(sel::sf_mu(p, 3, sel::Polarity::Outbending), WithinRel(mu_expected, 1e-12));
        CHECK_THAT(sel::sf_sigma(p, 3, sel::Polarity::Outbending), WithinRel(sg_expected, 1e-12));
    }

    SECTION("inbending sector 1 -- the other table is not a copy of the first") {
        constexpr double mu_a = 2.311088E-01, mu_b = 4.309933E-03, mu_c = -5.451627E-04,
                         mu_d = 4.965823E-06;
        const double p = 4.0;
        const double mu_expected = mu_a + mu_b * p + mu_c * p * p + mu_d * p * p * p;
        CHECK_THAT(sel::sf_mu(p, 1, sel::Polarity::Inbending), WithinRel(mu_expected, 1e-12));

        // And the two polarities really are different data at the same sector.
        CHECK(sel::sf_mu(p, 1, sel::Polarity::Inbending) !=
              sel::sf_mu(p, 1, sel::Polarity::Outbending));
    }

    SECTION("a fully worked value, computed offline") {
        // Sector 3 OB at p = 4.0:
        //   mu = 0.2258053 + 9.073840e-3*4 - 1.240745e-3*16 + 3.579321e-5*64
        //      = 0.2258053 + 0.03629536 - 0.01985192 + 0.0022907654
        //      = 0.2445395054
        CHECK_THAT(sel::sf_mu(4.0, 3, sel::Polarity::Outbending),
                   WithinAbs(0.2445395054, 1e-9));
    }
}

TEST_CASE("sampling fraction pass() is symmetric about mu and strict at the edge", "[sf]") {
    const int sector = 2;
    const double p = 5.0;
    const double n_sigma = shipped_cuts().electron.sf_n_sigma;  // 3.5, from config
    const auto pol = sel::Polarity::Outbending;

    const double mu = sel::sf_mu(p, sector, pol);
    const double sigma = sel::sf_sigma(p, sector, pol);
    const double half = n_sigma * sigma;
    REQUIRE(sigma > 0.0);  // otherwise the rest of this test says nothing

    SECTION("the centre passes") { CHECK(sel::pass(mu, p, sector, pol, n_sigma)); }

    SECTION("symmetric about mu") {
        const auto frac = GENERATE(0.0, 0.25, 0.5, 0.9, 0.999);
        INFO("offset = " << frac << " * n_sigma * sigma");
        CHECK(sel::pass(mu + frac * half, p, sector, pol, n_sigma) ==
              sel::pass(mu - frac * half, p, sector, pol, n_sigma));
        CHECK(sel::pass(mu + frac * half, p, sector, pol, n_sigma));
    }

    SECTION("rejects just outside the band, both sides") {
        const double eps = 1e-9;
        CHECK_FALSE(sel::pass(mu + half + eps, p, sector, pol, n_sigma));
        CHECK_FALSE(sel::pass(mu - half - eps, p, sector, pol, n_sigma));
    }

    SECTION("the boundary itself is EXCLUDED -- the comparison is strict") {
        CHECK_FALSE(sel::pass(mu + half, p, sector, pol, n_sigma));
        CHECK_FALSE(sel::pass(mu - half, p, sector, pol, n_sigma));
    }

    SECTION("just inside passes") {
        CHECK(sel::pass(std::nextafter(mu + half, mu), p, sector, pol, n_sigma));
        CHECK(sel::pass(std::nextafter(mu - half, mu), p, sector, pol, n_sigma));
    }
}

TEST_CASE("sampling fraction pass() rejects an invalid sector without throwing", "[sf]") {
    const auto pol = sel::Polarity::Outbending;
    const double p = 4.0;
    const double mu_valid = sel::sf_mu(p, 1, pol);  // a value that would pass in sector 1

    const auto sector = GENERATE(-1, 0, 7, 99);
    INFO("sector = " << sector);
    CHECK_FALSE(sel::pass(mu_valid, p, sector, pol, 3.5));
    // sf_params(), by contrast, throws: it has no "reject" answer to give.
    CHECK_THROWS_AS(sel::sf_params(sector, pol), std::out_of_range);
}

TEST_CASE("all six sectors are populated and distinct", "[sf]") {
    // A table half-filled by a copy-paste slip would still pass a single-sector
    // test. Every sector must be real data, and no two sectors are identical.
    const auto pol = GENERATE(sel::Polarity::Outbending, sel::Polarity::Inbending);
    for (int s = 1; s <= 6; ++s) {
        INFO("sector " << s);
        // Sampling fraction of the CLAS12 ECAL is ~0.25 and the width ~1e-2.
        CHECK(sel::sf_mu(4.0, s, pol) > 0.15);
        CHECK(sel::sf_mu(4.0, s, pol) < 0.35);
        CHECK(sel::sf_sigma(4.0, s, pol) > 0.0);
        for (int t = s + 1; t <= 6; ++t) {
            CHECK(sel::sf_mu(4.0, s, pol) != sel::sf_mu(4.0, t, pol));
        }
    }
}

TEST_CASE("polarity_from_string round-trips and rejects nonsense", "[sf]") {
    CHECK(sel::polarity_from_string("outbending") == sel::Polarity::Outbending);
    CHECK(sel::polarity_from_string("inbending") == sel::Polarity::Inbending);
    CHECK(std::string(sel::to_string(sel::Polarity::Outbending)) == "outbending");
    CHECK(std::string(sel::to_string(sel::Polarity::Inbending)) == "inbending");

    // The shipped config must name a polarity this code understands.
    CHECK_NOTHROW(sel::polarity_from_string(shipped_cuts().beam.polarity));

    CHECK_THROWS_AS(sel::polarity_from_string("OB"), std::runtime_error);
    CHECK_THROWS_AS(sel::polarity_from_string("Outbending"), std::runtime_error);
    CHECK_THROWS_AS(sel::polarity_from_string(""), std::runtime_error);
}

// ===========================================================================
// find_trigger_electron
// ===========================================================================

TEST_CASE("find_trigger_electron picks the first FD-status electron", "[electron]") {
    const pi0::Cuts& c = shipped_cuts();

    SECTION("a positive-status electron is ignored, the negative one is taken") {
        // Row 0: an electron with POSITIVE status -- not trigger-flagged.
        // Row 1: an electron with status -2100 -- the trigger electron.
        const RVec<int> pid{11, 11};
        const RVec<short> status{2100, -2100};
        const RVec<float> px{1.0f, 0.5f}, py{0.0f, 0.0f}, pz{3.0f, 2.0f};

        const auto found = sel::find_trigger_electron(pid, status, px, py, pz, c);
        REQUIRE(found.has_value());
        CHECK(*found == 1u);
    }

    SECTION("FIRST match wins, not the highest momentum") {
        // This ordering dependence is inherited from the old code and preserved
        // deliberately. Row 0 is softer than row 1; row 0 still wins.
        const RVec<int> pid{11, 11};
        const RVec<short> status{-2100, -2100};
        const RVec<float> px{0.1f, 5.0f}, py{0.0f, 0.0f}, pz{2.5f, 9.0f};

        const auto found = sel::find_trigger_electron(pid, status, px, py, pz, c);
        REQUIRE(found.has_value());
        CHECK(*found == 0u);
    }

    SECTION("non-electrons are skipped whatever their status") {
        const RVec<int> pid{22, 2212, -11, 11};
        const RVec<short> status{-2100, -2100, -2100, -2100};
        const RVec<float> px{1.0f, 1.0f, 1.0f, 1.0f}, py{0.0f, 0.0f, 0.0f, 0.0f},
            pz{3.0f, 3.0f, 3.0f, 3.0f};

        const auto found = sel::find_trigger_electron(pid, status, px, py, pz, c);
        REQUIRE(found.has_value());
        CHECK(*found == 3u);  // the positron at row 2 is NOT an electron
    }

    SECTION("a zero momentum vector disqualifies a row") {
        const RVec<int> pid{11, 11};
        const RVec<short> status{-2100, -2100};
        const RVec<float> px{0.0f, 0.5f}, py{0.0f, 0.0f}, pz{0.0f, 2.0f};

        const auto found = sel::find_trigger_electron(pid, status, px, py, pz, c);
        REQUIRE(found.has_value());
        CHECK(*found == 1u);
    }

    SECTION("the |status| window is [2000, 4000): inclusive low, EXCLUSIVE high") {
        // -1999 is below the band; -4000 is the first value ABOVE it, and must
        // be rejected -- the upper bound is exclusive. -2000 and -3999 are in.
        const auto probe = [&](short st) {
            const RVec<int> pid{11};
            const RVec<short> status{st};
            const RVec<float> px{1.0f}, py{0.0f}, pz{3.0f};
            return sel::find_trigger_electron(pid, status, px, py, pz, c).has_value();
        };

        CHECK_FALSE(probe(-1999));
        CHECK(probe(-2000));   // inclusive lower bound
        CHECK(probe(-3999));
        CHECK_FALSE(probe(-4000));  // EXCLUSIVE upper bound
        CHECK_FALSE(probe(-4001));
    }

    SECTION("an event with no qualifying electron yields nullopt") {
        const RVec<int> pid{2212, 22};
        const RVec<short> status{-2100, -2100};
        const RVec<float> px{1.0f, 1.0f}, py{0.0f, 0.0f}, pz{3.0f, 3.0f};
        CHECK_FALSE(sel::find_trigger_electron(pid, status, px, py, pz, c).has_value());
    }

    SECTION("an empty bank yields nullopt") {
        CHECK_FALSE(sel::find_trigger_electron({}, {}, {}, {}, {}, c).has_value());
    }

    SECTION("a truncated column does not cause an out-of-bounds read") {
        // A malformed bank must yield "no electron", not a crash.
        const RVec<int> pid{11, 11};
        const RVec<short> status{-2100};  // short by one row
        const RVec<float> px{0.0f, 1.0f}, py{0.0f, 0.0f}, pz{0.0f, 3.0f};
        // Row 0 has zero momentum and row 1 is beyond the shortest column.
        CHECK_FALSE(sel::find_trigger_electron(pid, status, px, py, pz, c).has_value());
    }
}

// ===========================================================================
// pass_electron
// ===========================================================================

TEST_CASE("pass_electron accepts a good electron", "[electron]") {
    const pi0::Cuts& c = shipped_cuts();
    const auto result = make_good_electron(c).apply(c);
    INFO("failed_at = " << (result.failed_at ? result.failed_at : "(none)"));
    CHECK(result.passed);
    CHECK(result.failed_at == nullptr);
}

TEST_CASE("each electron cut can be made to fail individually, and reports itself",
          "[electron]") {
    const pi0::Cuts& c = shipped_cuts();

    const auto expect_failure = [&](const GoodElectron& e, const char* stage) {
        const auto result = e.apply(c);
        INFO("expected stage: " << stage);
        INFO("actual stage:   " << (result.failed_at ? result.failed_at : "(passed)"));
        CHECK_FALSE(result.passed);
        REQUIRE(result.failed_at != nullptr);
        CHECK(std::string(result.failed_at) == stage);
    };

    SECTION("chi2pid, below the floor") {
        auto e = make_good_electron(c);
        e.chi2pid = c.electron.chi2pid_min - 0.1;
        expect_failure(e, sel::electron_stage::kChi2Pid);
    }

    SECTION("chi2pid, above the ceiling") {
        auto e = make_good_electron(c);
        e.chi2pid = c.electron.chi2pid_max + 0.1;
        expect_failure(e, sel::electron_stage::kChi2Pid);
    }

    SECTION("chi2pid, exactly on the boundary -- the cut is STRICT") {
        auto e = make_good_electron(c);
        e.chi2pid = c.electron.chi2pid_max;
        expect_failure(e, sel::electron_stage::kChi2Pid);
    }

    SECTION("momentum") {
        auto e = make_good_electron(c);
        e.p_gev = c.electron.min_momentum_gev - 0.1;  // 1.9 GeV
        // Re-centre the SF so this can only fail on momentum.
        e.sampling_fraction = sel::sf_mu(e.p_gev, e.pcal_sector, shipped_polarity());
        expect_failure(e, sel::electron_stage::kMomentum);
    }

    SECTION("momentum, exactly at the threshold -- STRICT") {
        auto e = make_good_electron(c);
        e.p_gev = c.electron.min_momentum_gev;  // exactly 2.0
        e.sampling_fraction = sel::sf_mu(e.p_gev, e.pcal_sector, shipped_polarity());
        expect_failure(e, sel::electron_stage::kMomentum);
    }

    SECTION("momentum: 1.0 GeV is REJECTED -- proof the cut is 2.0, not 0.8") {
        // If anyone ever 'restores' the old label's 0.8, this fails. A 1.0 GeV
        // electron passes a 0.8 cut and fails a 2.0 cut; the applied cut is 2.0.
        auto e = make_good_electron(c);
        e.p_gev = 1.0;
        e.sampling_fraction = sel::sf_mu(e.p_gev, e.pcal_sector, shipped_polarity());
        expect_failure(e, sel::electron_stage::kMomentum);
    }

    SECTION("vertex -- the caller's verdict, passed straight through") {
        auto e = make_good_electron(c);
        e.vertex_passed = false;
        expect_failure(e, sel::electron_stage::kVertex);
    }

    SECTION("sampling fraction, too high") {
        auto e = make_good_electron(c);
        const double sigma = sel::sf_sigma(e.p_gev, e.pcal_sector, shipped_polarity());
        e.sampling_fraction += (c.electron.sf_n_sigma + 0.01) * sigma;
        expect_failure(e, sel::electron_stage::kSamplingFraction);
    }

    SECTION("sampling fraction, too low") {
        auto e = make_good_electron(c);
        const double sigma = sel::sf_sigma(e.p_gev, e.pcal_sector, shipped_polarity());
        e.sampling_fraction -= (c.electron.sf_n_sigma + 0.01) * sigma;
        expect_failure(e, sel::electron_stage::kSamplingFraction);
    }

    SECTION("an unplaceable track fails AT the sampling-fraction stage") {
        auto e = make_good_electron(c);
        e.pcal_sector = 0;  // no valid PCAL sector
        expect_failure(e, sel::electron_stage::kSamplingFraction);
    }

    SECTION("PCAL fiducial, lv") {
        auto e = make_good_electron(c);
        e.pcal_lv_cm = c.electron.pcal_lv_min_cm - 0.1;
        expect_failure(e, sel::electron_stage::kPcalFiducial);
    }

    SECTION("PCAL fiducial, lw") {
        auto e = make_good_electron(c);
        e.pcal_lw_cm = c.electron.pcal_lw_min_cm - 0.1;
        expect_failure(e, sel::electron_stage::kPcalFiducial);
    }

    SECTION("PCAL fiducial, exactly on the threshold -- STRICT") {
        auto e = make_good_electron(c);
        e.pcal_lv_cm = c.electron.pcal_lv_min_cm;  // exactly 9.0
        expect_failure(e, sel::electron_stage::kPcalFiducial);
    }

    SECTION("PCAL fiducial: an electron at 10 cm PASSES the 9.0 cut") {
        // The mirror of the photon test below: at lv = 10 the electron passes
        // and a photon would fail. The two thresholds are different on purpose.
        auto e = make_good_electron(c);
        e.pcal_lv_cm = 10.0;
        e.pcal_lw_cm = 10.0;
        CHECK(e.apply(c).passed);
        CHECK(10.0 < c.photon.pcal_lv_min_cm);  // ... and 10 would fail as a photon
    }

    SECTION("DC edge R1") {
        auto e = make_good_electron(c);
        e.dc_edge_r1_cm = c.electron.dc_edge_r1_cm - 0.01;
        expect_failure(e, sel::electron_stage::kDcEdge);
    }

    SECTION("DC edge R2") {
        auto e = make_good_electron(c);
        e.dc_edge_r2_cm = c.electron.dc_edge_r2_cm - 0.01;
        expect_failure(e, sel::electron_stage::kDcEdge);
    }

    SECTION("DC edge R3") {
        auto e = make_good_electron(c);
        e.dc_edge_r3_cm = c.electron.dc_edge_r3_cm - 0.01;
        expect_failure(e, sel::electron_stage::kDcEdge);
    }

    SECTION("DC edge R3, exactly on the threshold -- STRICT") {
        auto e = make_good_electron(c);
        e.dc_edge_r3_cm = c.electron.dc_edge_r3_cm;  // exactly 8.75
        expect_failure(e, sel::electron_stage::kDcEdge);
    }
}

TEST_CASE("pass_electron short-circuits: the FIRST failing stage is reported", "[electron]") {
    const pi0::Cuts& c = shipped_cuts();

    // Spoil chi2pid AND the DC edge. chi2pid is applied first, so that is what
    // must be reported -- this is what makes the result usable as a sequential
    // cutflow rather than an arbitrary pick among the failures.
    auto e = make_good_electron(c);
    e.chi2pid = 99.0;
    e.dc_edge_r1_cm = 0.0;

    const auto result = e.apply(c);
    CHECK_FALSE(result.passed);
    REQUIRE(result.failed_at != nullptr);
    CHECK(std::string(result.failed_at) == sel::electron_stage::kChi2Pid);
}

TEST_CASE("cutflow labels are derived from the config, and cannot go stale", "[electron][label]") {
    const pi0::Cuts& c = shipped_cuts();

    SECTION("the momentum label says 2, and does NOT say 0.8") {
        // THE REGRESSION TEST FOR THE HEADLINE BUG. The old label was the
        // literal string "Momentum > 0.8 GeV" beside code applying 2.0.
        const std::string label = sel::electron_cutflow_label(sel::electron_stage::kMomentum, c);
        INFO("label = " << label);
        CHECK_THAT(label, ContainsSubstring("2"));
        CHECK_THAT(label, !ContainsSubstring("0.8"));
    }

    SECTION("the label tracks the config rather than the code") {
        // Move the threshold; the label must move with it. This is the property
        // that makes the old bug unreproducible: nothing renders a number that
        // is not the number being applied.
        pi0::Cuts moved = c;
        moved.electron.min_momentum_gev = 3.25;
        CHECK_THAT(sel::electron_cutflow_label(sel::electron_stage::kMomentum, moved),
                   ContainsSubstring("3.25"));
    }

    SECTION("every stage has a label, and none is empty") {
        for (const char* stage : sel::kElectronStages) {
            INFO("stage = " << stage);
            CHECK_FALSE(sel::electron_cutflow_label(stage, c).empty());
        }
    }

    SECTION("the PCAL label says 9, not the photon's 14") {
        CHECK_THAT(sel::electron_cutflow_label(sel::electron_stage::kPcalFiducial, c),
                   ContainsSubstring("9"));
        CHECK_THAT(sel::electron_cutflow_label(sel::electron_stage::kPcalFiducial, c),
                   !ContainsSubstring("14"));
    }

    SECTION("an unknown stage throws rather than returning a plausible string") {
        CHECK_THROWS_AS(sel::electron_cutflow_label("not_a_stage", c), std::invalid_argument);
    }

    SECTION("every failed_at value is a legal label input") {
        // Guards against a stage constant existing that kElectronStages forgot.
        auto e = make_good_electron(c);
        e.chi2pid = 99.0;
        const auto result = e.apply(c);
        REQUIRE(result.failed_at != nullptr);
        CHECK_NOTHROW(sel::electron_cutflow_label(result.failed_at, c));
    }
}

// ===========================================================================
// pass_photon
// ===========================================================================

namespace {

/// A photon at theta degrees in the x-z plane with energy `e_gev`. Massless, so
/// |p| = E.
struct PhotonArgs {
    double px{}, py{}, pz{};
    double pcal_energy_gev = 0.5;
    double pcal_lv_cm = 20.0;
    double pcal_lw_cm = 20.0;
    double beta = 1.0;
    double gbt_score = 0.99;
    int pid = sel::kPdgPhoton;

    [[nodiscard]] bool apply(const pi0::Cuts& cuts) const {
        return sel::pass_photon(pid, px, py, pz, pcal_energy_gev, pcal_lv_cm, pcal_lw_cm, beta,
                                gbt_score, cuts);
    }
};

PhotonArgs photon_at(double theta_deg, double e_gev) {
    PhotonArgs g;
    const double theta_rad = theta_deg * pi0::kDegToRad;
    g.px = e_gev * std::sin(theta_rad);
    g.py = 0.0;
    g.pz = e_gev * std::cos(theta_rad);
    return g;
}

}  // namespace

TEST_CASE("pass_photon accepts a good photon", "[photon]") {
    CHECK(photon_at(20.0, 2.0).apply(shipped_cuts()));
}

TEST_CASE("each photon cut can be made to fail individually", "[photon]") {
    const pi0::Cuts& c = shipped_cuts();

    SECTION("pid must be 22") {
        auto g = photon_at(20.0, 2.0);
        g.pid = 2112;  // a neutron: the very background the GBT exists to reject
        CHECK_FALSE(g.apply(c));
    }

    SECTION("energy floor rejects below and accepts above") {
        // The INCLUSIVITY of this bound (>=, unlike most cuts here) is pinned
        // exactly further down, on passes_gbt_prefilter -- photon_at() builds
        // |p| out of E*sin/E*cos, so it cannot place the energy exactly on the
        // floor and this section does not pretend to.
        CHECK(photon_at(20.0, c.photon.min_energy_gev * 1.01).apply(c));
        CHECK_FALSE(photon_at(20.0, c.photon.min_energy_gev * 0.5).apply(c));
    }

    SECTION("PCAL energy must be positive -- the shower starts in the preshower") {
        auto g = photon_at(20.0, 2.0);
        g.pcal_energy_gev = 0.0;
        CHECK_FALSE(g.apply(c));
    }

    SECTION("theta window rejects outside, accepts inside") {
        CHECK(photon_at(c.photon.theta_min_deg + 0.5, 2.0).apply(c));
        CHECK(photon_at(c.photon.theta_max_deg - 0.5, 2.0).apply(c));
        CHECK_FALSE(photon_at(c.photon.theta_min_deg - 0.5, 2.0).apply(c));
        CHECK_FALSE(photon_at(c.photon.theta_max_deg + 0.5, 2.0).apply(c));
    }

    SECTION("theta window is INCLUSIVE at both ends") {
        // Tested on passes_gbt_prefilter, which takes theta_deg directly, and
        // NOT by building a photon "at exactly 5 degrees" with photon_at().
        //
        // That was the first attempt and it failed, correctly: the
        // deg -> rad -> (px, pz) -> acos -> deg round-trip lands about 8e-15
        // BELOW the angle asked for (4.999999999999992 for 5.0), so such a
        // photon is genuinely outside an inclusive >= 5.0 window and the code
        // was right to reject it. A test cannot construct an exact boundary
        // through trig, so it must not claim to. Here the boundary is exact
        // because it is passed in as a number.
        //
        // Nothing physical rides on the boundary itself -- a reconstructed
        // photon landing on it to 1e-15 is measure-zero -- but which side the
        // bound falls on is documented in cuts.json, so it is pinned.
        const double e_ok = 2.0;
        const double epcal_ok = 0.5;

        CHECK(sel::passes_gbt_prefilter(e_ok, epcal_ok, c.photon.theta_min_deg, c));  // 5: IN
        CHECK(sel::passes_gbt_prefilter(e_ok, epcal_ok, c.photon.theta_max_deg, c));  // 35: IN
        CHECK_FALSE(sel::passes_gbt_prefilter(
            e_ok, epcal_ok, std::nextafter(c.photon.theta_min_deg, 0.0), c));
        CHECK_FALSE(sel::passes_gbt_prefilter(
            e_ok, epcal_ok, std::nextafter(c.photon.theta_max_deg, 90.0), c));
    }

    SECTION("the energy floor is INCLUSIVE, checked exactly") {
        // Same reasoning: photon_at() cannot place |p| exactly on the floor
        // either (it is built from E*sin/E*cos and re-normed), so the exact
        // boundary is checked on the pre-filter, which takes the energy
        // directly.
        CHECK(sel::passes_gbt_prefilter(c.photon.min_energy_gev, 0.5, 20.0, c));  // 0.2: IN
        CHECK_FALSE(sel::passes_gbt_prefilter(std::nextafter(c.photon.min_energy_gev, 0.0), 0.5,
                                              20.0, c));
    }

    SECTION("GBT score, strictly above the threshold") {
        auto g = photon_at(20.0, 2.0);
        g.gbt_score = c.photon.gbt_threshold + 0.01;
        CHECK(g.apply(c));
        g.gbt_score = c.photon.gbt_threshold;  // exactly 0.78 -- STRICT, so out
        CHECK_FALSE(g.apply(c));
        g.gbt_score = c.photon.gbt_threshold - 0.01;
        CHECK_FALSE(g.apply(c));
    }

    SECTION("PCAL fiducial is the photon's 14.0, NOT the electron's 9.0") {
        auto g = photon_at(20.0, 2.0);
        // 10 cm clears the electron's 9.0 but must NOT clear the photon's 14.0.
        // If someone ever unifies the two thresholds, this is what fails.
        g.pcal_lv_cm = 10.0;
        g.pcal_lw_cm = 10.0;
        CHECK_FALSE(g.apply(c));
        CHECK(10.0 > c.electron.pcal_lv_min_cm);  // ... and it would pass as an electron

        g.pcal_lv_cm = c.photon.pcal_lv_min_cm;  // exactly 14.0 -- STRICT, so out
        g.pcal_lw_cm = 20.0;
        CHECK_FALSE(g.apply(c));

        g.pcal_lv_cm = 14.1;
        CHECK(g.apply(c));
    }

    SECTION("beta window is STRICT at both ends") {
        auto g = photon_at(20.0, 2.0);
        g.beta = c.photon.beta_min;  // exactly 0.9
        CHECK_FALSE(g.apply(c));
        g.beta = c.photon.beta_max;  // exactly 1.1
        CHECK_FALSE(g.apply(c));
        g.beta = 0.95;
        CHECK(g.apply(c));
        g.beta = 0.5;
        CHECK_FALSE(g.apply(c));
    }
}

TEST_CASE("the photon theta cut is in DEGREES", "[photon][units]") {
    // THE DEGREE/RADIAN GUARD. The old pre-filter took `double theta` in
    // radians and converted inside -- an unmarked radian at a public boundary,
    // and one of the two independent degree/radian bugs in the old analysis.
    //
    // A photon at 20 degrees is 0.349 radians. If this API were secretly taking
    // radians, 0.349 would fall below the 5 (deg) floor and the photon would be
    // rejected; if a caller passed radians to this degree API, likewise. So a
    // 20-degree photon passing, and a 0.349-"degree" photon failing, pins the
    // unit down.
    const pi0::Cuts& c = shipped_cuts();

    CHECK(photon_at(20.0, 2.0).apply(c));

    // 20 degrees expressed in radians, fed to a degree API: must be rejected.
    CHECK_FALSE(photon_at(20.0 * pi0::kDegToRad, 2.0).apply(c));

    // And photon_theta_deg really returns degrees.
    const auto g = photon_at(20.0, 2.0);
    CHECK_THAT(sel::photon_theta_deg(g.px, g.py, g.pz), WithinAbs(20.0, 1e-9));
    CHECK_THAT(sel::photon_theta_deg(0.0, 0.0, 1.0), WithinAbs(0.0, 1e-9));
    CHECK_THAT(sel::photon_theta_deg(1.0, 0.0, 0.0), WithinAbs(90.0, 1e-9));
}

TEST_CASE("photon energy is |p|: massless", "[photon]") {
    // E = |p| with no calorimeter correction anywhere in this chain.
    const pi0::Cuts& c = shipped_cuts();

    // A photon whose components give |p| just under the floor is rejected, and
    // just over is accepted -- i.e. the floor is applied to |p|, not to pz or
    // to any single component.
    CHECK_FALSE(photon_at(20.0, 0.199).apply(c));
    CHECK(photon_at(20.0, 0.201).apply(c));

    // A zero momentum vector has no direction; it must be rejected, not NaN its
    // way through.
    PhotonArgs zero;
    zero.px = zero.py = zero.pz = 0.0;
    CHECK_FALSE(zero.apply(c));
    CHECK(std::isnan(sel::photon_theta_deg(0.0, 0.0, 0.0)));
}

TEST_CASE("pass_photon_scored evaluates the model only when the pre-filter passes",
          "[photon]") {
    // The pre-filter exists so the GBT is not run on most clusters. This pins
    // that: the callable must not be invoked for a photon the cheap cuts
    // already rejected.
    const pi0::Cuts& c = shipped_cuts();

    int calls = 0;
    const auto score = [&calls] {
        ++calls;
        return 0.99;
    };

    SECTION("not called for a non-photon") {
        auto g = photon_at(20.0, 2.0);
        CHECK_FALSE(sel::pass_photon_scored(2112, g.px, g.py, g.pz, g.pcal_energy_gev,
                                            g.pcal_lv_cm, g.pcal_lw_cm, g.beta, score, c));
        CHECK(calls == 0);
    }

    SECTION("not called for a photon below the energy floor") {
        auto g = photon_at(20.0, 0.05);
        CHECK_FALSE(sel::pass_photon_scored(g.pid, g.px, g.py, g.pz, g.pcal_energy_gev,
                                            g.pcal_lv_cm, g.pcal_lw_cm, g.beta, score, c));
        CHECK(calls == 0);
    }

    SECTION("not called for a photon outside the theta window") {
        auto g = photon_at(45.0, 2.0);
        CHECK_FALSE(sel::pass_photon_scored(g.pid, g.px, g.py, g.pz, g.pcal_energy_gev,
                                            g.pcal_lv_cm, g.pcal_lw_cm, g.beta, score, c));
        CHECK(calls == 0);
    }

    SECTION("called exactly once for a photon that reaches the classifier") {
        auto g = photon_at(20.0, 2.0);
        CHECK(sel::pass_photon_scored(g.pid, g.px, g.py, g.pz, g.pcal_energy_gev, g.pcal_lv_cm,
                                      g.pcal_lw_cm, g.beta, score, c));
        CHECK(calls == 1);
    }

    SECTION("agrees with the eager overload") {
        const auto g = photon_at(20.0, 2.0);
        const auto probe = GENERATE(0.5, 0.78, 0.9);
        INFO("score = " << probe);
        CHECK(sel::pass_photon_scored(g.pid, g.px, g.py, g.pz, g.pcal_energy_gev, g.pcal_lv_cm,
                                      g.pcal_lw_cm, g.beta, [probe] { return probe; }, c) ==
              sel::pass_photon(g.pid, g.px, g.py, g.pz, g.pcal_energy_gev, g.pcal_lv_cm,
                               g.pcal_lw_cm, g.beta, probe, c));
    }
}
