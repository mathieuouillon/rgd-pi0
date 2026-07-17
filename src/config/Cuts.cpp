#include "config/Cuts.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "core/Constants.hpp"

// JSON reader: nlohmann/json.
//
// Chosen over a hand-rolled parser because the meson wrap is genuinely trivial
// (external/nlohmann_json.wrap, `meson wrap install nlohmann_json`, with the
// system pkg-config package used when present). It is header-only and PRIVATE
// to this translation unit: Cuts.hpp does not include it, so nothing else in
// the project acquires the dependency.
//
// nlohmann is used ONLY as a parser. Every key access below goes through the
// Reader helper, never through operator[] -- nlohmann's operator[] on a const
// object is UB for a missing key, and on a non-const object it default-inserts
// null, which is precisely the silent-default failure this file exists to
// prevent. `at()` throws, but throws a message naming only the leaf key; Reader
// names the full path and the file.

namespace pi0 {
namespace {

using nlohmann::json;

/// Navigates dotted key paths in the parsed document, and throws with the full
/// path and the file name on anything unexpected.
class Reader {
   public:
    Reader(const json& root, std::string file) : m_root(root), m_file(std::move(file)) {}

    /// \param path dotted key path, e.g. "electron.sampling_fraction.n_sigma".
    [[nodiscard]] const json& node(const std::string& path) const {
        const json* cur = &m_root;
        std::string walked;
        std::istringstream parts(path);
        std::string part;
        while (std::getline(parts, part, '.')) {
            if (!cur->is_object()) {
                fail(path, "is unreachable: '" + walked + "' is not an object");
            }
            const auto it = cur->find(part);
            walked += (walked.empty() ? "" : ".") + part;
            if (it == cur->end()) {
                fail(path, "is MISSING (no '" + walked + "')");
            }
            cur = &(*it);
        }
        return *cur;
    }

    [[nodiscard]] double num(const std::string& path) const {
        const json& n = node(path);
        if (!n.is_number()) fail(path, wrong_type("a number", n));
        return n.get<double>();
    }

    [[nodiscard]] int integer(const std::string& path) const {
        const json& n = node(path);
        if (!n.is_number_integer()) fail(path, wrong_type("an integer", n));
        return n.get<int>();
    }

    [[nodiscard]] bool boolean(const std::string& path) const {
        const json& n = node(path);
        if (!n.is_boolean()) fail(path, wrong_type("a boolean", n));
        return n.get<bool>();
    }

    [[nodiscard]] std::string str(const std::string& path) const {
        const json& n = node(path);
        if (!n.is_string()) fail(path, wrong_type("a string", n));
        return n.get<std::string>();
    }

    [[nodiscard]] const json& array(const std::string& path) const {
        const json& n = node(path);
        if (!n.is_array()) fail(path, wrong_type("an array", n));
        return n;
    }

    /// A number read from an arbitrary node rather than from a path off the
    /// root -- for array elements, which node() cannot address.
    ///
    /// `display_path` is what the error message names, so it must spell the
    /// full path the reader would have walked (e.g.
    /// "vertex.targets.cxc.peaks[1].mu_cm"). Same discipline as node(): the
    /// message is useless if it names only the leaf.
    [[nodiscard]] double num_in(const json& parent, const std::string& key,
                                const std::string& display_path) const {
        if (!parent.is_object()) fail(display_path, "is unreachable: its parent is not an object");
        const auto it = parent.find(key);
        if (it == parent.end()) fail(display_path, "is MISSING");
        if (!it->is_number()) fail(display_path, wrong_type("a number", *it));
        return it->get<double>();
    }

    [[noreturn]] void fail(const std::string& path, const std::string& why) const {
        throw std::runtime_error("Cuts::load: " + m_file + ": key '" + path + "' " + why);
    }

   private:
    [[nodiscard]] static std::string wrong_type(const std::string& wanted, const json& got) {
        return "is not " + wanted + " (found " + std::string(got.type_name()) + ")";
    }

    const json& m_root;
    std::string m_file;
};

/// True when two doubles agree to within a double round-trip of each other.
///
/// This is NOT a physics tolerance. It exists only so that "the same decimal
/// literal, parsed by two different parsers" compares equal; any real drift
/// between a constant and its config twin is many orders of magnitude larger.
bool same_number(double a, double b) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= 1e-12 * scale;
}

/// core/Constants.hpp and config/cuts.json both spell out the beam energy, the
/// proton mass and the pi0 mass. Constants.hpp's own header warns that "a
/// duplicated constant that drifts out of sync is the classic way an analysis
/// silently changes physics", and cuts.json was written to end exactly that
/// (the old tree carried THREE beam energies that agreed by coincidence).
///
/// The duplication is nonetheless real: the physics code reads the constant,
/// the config documents the value, and nothing forces them to be the same
/// number. So force it here. If they ever disagree the program refuses to
/// start, which is the only outcome that cannot silently change physics.
void require_agrees_with_constant(const Reader& r, const std::string& path, double constant,
                                  const std::string& constant_name) {
    const double from_json = r.num(path);
    if (!same_number(from_json, constant)) {
        std::ostringstream os;
        os.precision(17);
        os << "disagrees with core/Constants.hpp: cuts.json says " << from_json << " but "
           << constant_name << " is " << constant
           << ". These are the same physical quantity in two places. Fix the one that is "
              "wrong -- do not make the code prefer one silently.";
        r.fail(path, os.str());
    }
}

/// Parse one /vertex/targets/<target> block.
///
/// `base` is that block's dotted path, e.g. "vertex.targets.cu". Every key is
/// required: there is no default window, and a target whose block is incomplete
/// is a target whose selection nobody can state.
///
/// THE NAMING TRAP IS UNREACHABLE FROM HERE, BY DESIGN. cuts.json lists each
/// target's peaks explicitly instead of naming the fit's 'lo'/'hi' columns, so
/// "which foil is Cu" is a number in the file (-7.861, the more negative,
/// upstream one), not a routing decision in code that could be inverted while
/// still compiling and still filling histograms. See vertex/VzCorrector.hpp.
pi0::vertex::VzTargetCuts read_vertex_target(const Reader& r, const std::string& base) {
    using pi0::vertex::VzPeak;
    using pi0::vertex::VzRule;
    using pi0::vertex::VzTargetCuts;

    VzTargetCuts t;
    t.correction_enabled = r.boolean(base + ".correction_enabled");

    const std::string rule = r.str(base + ".rule");
    if (rule == "raw_window") {
        t.rule = VzRule::RawWindow;
        t.vz_min_cm = r.num(base + ".vz_min_cm");
        t.vz_max_cm = r.num(base + ".vz_max_cm");
    } else if (rule == "corrected_peaks") {
        t.rule = VzRule::CorrectedPeaks;
        t.n_sigma = r.num(base + ".n_sigma");

        const json& peaks = r.array(base + ".peaks");
        if (peaks.empty()) {
            r.fail(base + ".peaks",
                   "is an empty array. A corrected_peaks target with no peaks accepts NOTHING, "
                   "which is a selection nobody meant to write.");
        }
        for (std::size_t i = 0; i < peaks.size(); ++i) {
            const std::string at = base + ".peaks[" + std::to_string(i) + "]";
            VzPeak pk;
            pk.mu_cm = r.num_in(peaks[i], "mu_cm", at + ".mu_cm");
            pk.sigma_cm = r.num_in(peaks[i], "sigma_cm", at + ".sigma_cm");
            t.peaks.push_back(pk);
        }
    } else {
        r.fail(base + ".rule",
               "is '" + rule + "', which is neither 'raw_window' nor 'corrected_peaks'");
    }

    // The rule says WHICH v_z is being cut, and correction_enabled says whether
    // that v_z was corrected. They are two keys describing one decision, so a
    // config that disagrees with itself would have the skim cut a raw vertex
    // against peaks fitted on the corrected distribution -- a window silently
    // drawn around the wrong quantity. Refuse instead.
    const bool wants_correction = (t.rule == VzRule::CorrectedPeaks);
    if (t.correction_enabled != wants_correction) {
        r.fail(base + ".correction_enabled",
               std::string("is ") + (t.correction_enabled ? "true" : "false") + " but rule is '" +
                   rule + "'. The peaks of 'corrected_peaks' are fitted on the CORRECTED v_z and "
                          "'raw_window' bounds the RAW v_z, so these two keys must agree.");
    }
    return t;
}

/// Read one edge array of /mixing/pool_grid into a Grid1D.
///
/// These are the CONFIGURATION PRODUCT EDGES and they are hand-authored, so
/// unlike Grid A's they genuinely do belong in this file. `Grid1D` is reused as
/// the axis type (it is the project's only one) but the EDGES are this block's
/// own -- see stageB_bin/PoolGrid.hpp on why the pool grid and Grid A must not be
/// unified despite currently holding identical numbers.
///
/// Monotonicity and finiteness are NOT checked here; PoolGridCuts::validate()
/// owns that and load() calls it below, so the rule lives in one place and also
/// covers a grid assembled in code that never came through this loader.
Grid1D read_pool_axis(const Reader& r, const std::string& path, const std::string& name) {
    Grid1D g;
    g.name = name;
    const json& arr = r.array(path);
    g.edges.reserve(arr.size());
    for (std::size_t i = 0; i < arr.size(); ++i) {
        const std::string at = path + "[" + std::to_string(i) + "]";
        if (!arr[i].is_number()) {
            r.fail(at, "is not a number (found " + std::string(arr[i].type_name()) + ")");
        }
        g.edges.push_back(arr[i].get<double>());
    }
    return g;
}

/// Read /mixing/pool_grid/n_photons_classes.
///
/// `max` is REQUIRED but may be null, and null means OPEN-ENDED -- the >=4 class
/// absorbs every higher multiplicity. Same discipline as extraction.fit_bounds's
/// null amp_max: a null is a real state that the type represents (std::optional),
/// never a sentinel and never a large finite stand-in. An ABSENT max is a
/// different thing from a null one and is rejected: "no upper bound" is a claim
/// the config has to make out loud, not one a reader infers from silence.
std::vector<MultClass> read_mult_classes(const Reader& r, const std::string& path) {
    std::vector<MultClass> out;
    const json& arr = r.array(path);
    for (std::size_t i = 0; i < arr.size(); ++i) {
        const std::string at = path + "[" + std::to_string(i) + "]";
        const json& e = arr[i];
        if (!e.is_object()) {
            r.fail(at, "is not an object (found " + std::string(e.type_name()) + ")");
        }

        MultClass c;

        const auto lit = e.find("label");
        if (lit == e.end()) r.fail(at + ".label", "is MISSING");
        if (!lit->is_string()) {
            r.fail(at + ".label", "is not a string (found " + std::string(lit->type_name()) + ")");
        }
        c.label = lit->get<std::string>();

        const auto mit = e.find("min");
        if (mit == e.end()) r.fail(at + ".min", "is MISSING");
        if (!mit->is_number_integer()) {
            r.fail(at + ".min",
                   "is not an integer (found " + std::string(mit->type_name()) + ")");
        }
        c.min = mit->get<int>();

        const auto xit = e.find("max");
        if (xit == e.end()) {
            r.fail(at + ".max",
                   "is MISSING. An open-ended class must say so with an explicit null, because a "
                   "missing key and an unbounded class are different claims and only one of them "
                   "is deliberate.");
        }
        if (xit->is_null()) {
            c.max.reset();  // open-ended
        } else if (xit->is_number_integer()) {
            c.max = xit->get<int>();
        } else {
            r.fail(at + ".max", "is neither an integer nor null (found " +
                                    std::string(xit->type_name()) +
                                    "). null means OPEN-ENDED; there is no other spelling of it.");
        }

        out.push_back(std::move(c));
    }
    return out;
}

}  // namespace

Cuts Cuts::load(const std::string& json_path) {
    std::ifstream in(json_path);
    if (!in) {
        throw std::runtime_error("Cuts::load: cannot open '" + json_path + "'");
    }

    json root;
    try {
        in >> root;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("Cuts::load: '" + json_path + "' is not valid JSON: " + e.what());
    }
    if (!root.is_object()) {
        throw std::runtime_error("Cuts::load: '" + json_path + "' is valid JSON but not an object");
    }

    const Reader r(root, json_path);
    Cuts c;

    // -- beam ---------------------------------------------------------------
    c.beam.energy_gev = r.num("beam.energy_gev");
    c.beam.polarity = r.str("beam.polarity");
    if (c.beam.polarity != "inbending" && c.beam.polarity != "outbending") {
        r.fail("beam.polarity",
               "is '" + c.beam.polarity + "', which is neither 'inbending' nor 'outbending'");
    }

    // -- electron -----------------------------------------------------------
    c.electron.chi2pid_min = r.num("electron.chi2pid_min");
    c.electron.chi2pid_max = r.num("electron.chi2pid_max");
    c.electron.min_momentum_gev = r.num("electron.min_momentum_gev");
    c.electron.sf_n_sigma = r.num("electron.sampling_fraction.n_sigma");
    c.electron.pcal_lv_min_cm = r.num("electron.pcal_lv_min_cm");
    c.electron.pcal_lw_min_cm = r.num("electron.pcal_lw_min_cm");
    c.electron.dc_edge_r1_cm = r.num("electron.dc_edge_r1_cm");
    c.electron.dc_edge_r2_cm = r.num("electron.dc_edge_r2_cm");
    c.electron.dc_edge_r3_cm = r.num("electron.dc_edge_r3_cm");
    c.electron.status_min = r.integer("electron.trigger.status_abs_min");
    c.electron.status_max = r.integer("electron.trigger.status_abs_max");

    // -- dis ----------------------------------------------------------------
    // q2_max, w_max and y_min are deliberately absent from cuts.json (none is
    // reachable at 10.53 GeV, so none ever rejected an event). Do not add reads
    // for them here expecting to find something.
    c.dis.q2_min = r.num("dis.q2_min");
    c.dis.w_min = r.num("dis.w_min");
    c.dis.y_max = r.num("dis.y_max");

    // -- photon -------------------------------------------------------------
    c.photon.gbt_threshold = r.num("photon.gbt_threshold");
    c.photon.min_energy_gev = r.num("photon.min_energy_gev");
    c.photon.theta_min_deg = r.num("photon.theta_min_deg");
    c.photon.theta_max_deg = r.num("photon.theta_max_deg");
    c.photon.pcal_lv_min_cm = r.num("photon.pcal_lv_min_cm");
    c.photon.pcal_lw_min_cm = r.num("photon.pcal_lw_min_cm");
    c.photon.beta_min = r.num("photon.beta_min");
    c.photon.beta_max = r.num("photon.beta_max");
    c.photon.gbt_pass = r.integer("photon.gbt_pass");
    c.photon.allow_rga_fallback = r.boolean("photon.allow_rga_fallback");

    // -- vertex -------------------------------------------------------------
    // The windows the vz cut is drawn against. Until these were read here they
    // were a hard-coded table in vertex/VzCorrector.cpp while this block sat in
    // cuts.json unread -- the exact "declared, set, and never read" defect this
    // file exists to end, with the added twist that the table also carried
    // INBENDING windows the file never claimed to have.
    //
    // That table is gone, so there is nothing to fall back to, so the polarity
    // must be checked rather than assumed: cuts.json's beam block says every
    // polarity-dependent number in the file is the outbending set, and applying
    // an outbending window to inbending data is precisely the silent
    // wrong-number failure the table's deletion was meant to prevent. Refuse.
    if (c.beam.polarity != "outbending") {
        r.fail("beam.polarity",
               "is 'inbending', but /vertex carries the OUTBENDING windows only -- as beam's own "
               "_polarity_comment says of every polarity-dependent number in this file. There is "
               "no inbending window anywhere in the project to fall back to (the hard-coded table "
               "that once held one is deleted, deliberately). Adding inbending requires a "
               "per-polarity /vertex block, not an edit of these numbers.");
    }
    c.vertex.ld2 = read_vertex_target(r, "vertex.targets.ld2");
    c.vertex.cxc = read_vertex_target(r, "vertex.targets.cxc");
    c.vertex.cu = read_vertex_target(r, "vertex.targets.cu");
    c.vertex.sn = read_vertex_target(r, "vertex.targets.sn");

    // -- pairing ------------------------------------------------------------
    c.pairing.mass_window_gev = r.num("pairing.mass_window_gev");
    c.pairing.min_mgg_gev = r.num("pairing.min_mgg_gev");
    c.pairing.e_gamma_min_angle_deg = r.num("pairing.e_gamma_min_angle_deg");
    c.pairing.open_a_deg = r.num("pairing.opening_angle.a_deg");
    c.pairing.open_b_inv_gev = r.num("pairing.opening_angle.b_inv_gev");
    c.pairing.open_offset_deg = r.num("pairing.opening_angle.offset_deg");

    // -- sidis ---------------------------------------------------------------
    //
    // Under /pairing in the JSON, in Cuts::sidis in the struct: z = E_pi0 / nu
    // needs the event's nu, so it is not something find_gg_pairs could apply,
    // and putting it in PairingCuts would imply it does. Same
    // key-path-is-not-field-name mapping as electron.sf_n_sigma, and this
    // function is the only place it appears.
    //
    // These two keys were DECLARED AND NEVER READ until make_grid needed them --
    // this file's own defect, of exactly the kind its header warns about. See
    // the note on Cuts::sidis.
    c.sidis.z_min = r.num("pairing.z_min");
    c.sidis.z_max = r.num("pairing.z_max");
    if (!(c.sidis.z_min < c.sidis.z_max)) {
        r.fail("pairing.z_max", "is " + std::to_string(c.sidis.z_max) + ", which is not above pairing.z_min (" +
                                    std::to_string(c.sidis.z_min) + "). The window z_min < z < z_max is empty, so the "
                                    "analysis would select no pi0 at all.");
    }

    // -- binning -------------------------------------------------------------
    //
    // The grid SHAPE only. The EDGES are computed from data by
    // src/tools/make_grid and live in config/binning/*.json with their own
    // provenance block -- see the note on Cuts::binning for why they are not,
    // and must not be, in this file.
    const auto read_n_bins = [&](const std::string& path) {
        const int n = r.integer(path);
        // Zero bins is not a degenerate configuration to tolerate, it is a typo.
        // One is legal and occasionally wanted (integrating an axis away), so
        // the floor is 1 rather than 2.
        if (n < 1) r.fail(path, "is " + std::to_string(n) + ", but a grid axis needs at least one bin.");
        return n;
    };
    c.binning.n_q2 = read_n_bins("binning.grid_a.n_q2");
    c.binning.n_xb = read_n_bins("binning.grid_a.n_xb");
    c.binning.n_z = read_n_bins("binning.grid_b.n_z");
    c.binning.n_pt2 = read_n_bins("binning.grid_b.n_pt2");

    // -- mgg_histogram -------------------------------------------------------
    //
    // The m_gg axis the same-event spectrum, the mixed-event spectrum and the
    // per-(4D bin, m_gg bin) kinematic sums all share. DECLARED AND NEVER READ
    // until stageB_bin consumed it -- see the note on Cuts::mgg_histogram.
    c.mgg_histogram.min_gev = r.num("mgg_histogram.min_gev");
    c.mgg_histogram.max_gev = r.num("mgg_histogram.max_gev");
    c.mgg_histogram.bins = r.integer("mgg_histogram.bins");
    if (!(c.mgg_histogram.min_gev < c.mgg_histogram.max_gev)) {
        r.fail("mgg_histogram.max_gev",
               "is " + std::to_string(c.mgg_histogram.max_gev) + ", which is not above mgg_histogram.min_gev (" +
                   std::to_string(c.mgg_histogram.min_gev) +
                   "). An empty or inverted mass axis does not throw downstream -- it silently histograms nothing.");
    }
    if (c.mgg_histogram.bins < 1) {
        r.fail("mgg_histogram.bins",
               "is " + std::to_string(c.mgg_histogram.bins) + ", but the mass axis needs at least one bin.");
    }

    // -- bsa -----------------------------------------------------------------
    //
    // The phi_h axis. Its range is NOT a free choice: SidisKin::phi_h_deg is
    // atan2-derived and lands in [-180, 180], so an axis over [0, 360] would put
    // every negative phi_h off the grid, and Grid1D::find() returns -1 there --
    // i.e. half the data would be dropped BY AN AXIS, silently, with no error
    // raised anywhere. Checked against the producer's range rather than trusted.
    c.bsa.n_phi_bins = r.integer("bsa.n_phi_bins");
    if (c.bsa.n_phi_bins < 1) {
        r.fail("bsa.n_phi_bins",
               "is " + std::to_string(c.bsa.n_phi_bins) + ", but the phi_h axis needs at least one bin.");
    }
    c.bsa.phi_min_deg = r.num("bsa.phi_min_deg");
    c.bsa.phi_max_deg = r.num("bsa.phi_max_deg");
    if (!(c.bsa.phi_min_deg < c.bsa.phi_max_deg)) {
        r.fail("bsa.phi_max_deg", "is " + std::to_string(c.bsa.phi_max_deg) + ", which is not above bsa.phi_min_deg (" +
                                      std::to_string(c.bsa.phi_min_deg) + ").");
    }
    if (c.bsa.phi_min_deg > -180.0 || c.bsa.phi_max_deg < 180.0) {
        r.fail("bsa", "spans [" + std::to_string(c.bsa.phi_min_deg) + ", " + std::to_string(c.bsa.phi_max_deg) +
                          "] degrees, which does not cover the full [-180, 180] that pi0::kin::compute_sidis "
                          "produces (phi_h is atan2-derived, Trento convention). A phi_h outside the axis lands off "
                          "the grid and is DISCARDED WITHOUT AN ERROR, so a narrow axis does not fail -- it quietly "
                          "deletes data. Widen the axis; do not run with a lossy one.");
    }

    // -- mixing --------------------------------------------------------------
    //
    // The frozen donor pool's grid and knobs (stageB_bin/DonorPool.hpp). This
    // whole block was DECLARED AND NEVER READ until DonorPool consumed it -- see
    // the note on Cuts::mixing.
    //
    // NOTE the pool grid's Q^2/x_B edges are NOT Cuts::binning's Grid A, however
    // identical the two look today. Nothing below cross-checks them against Grid
    // A's shape or its edges, deliberately: they are independent grids that
    // currently coincide only because Grid A's file is a placeholder, and a check
    // asserting the coincidence would pass every test until make_grid first ran
    // on real data and fail forever after.
    const int donors = r.integer("mixing.donors_per_bin");
    if (donors < 1) {
        r.fail("mixing.donors_per_bin",
               "is " + std::to_string(donors) +
                   ", but a pool of that depth holds nothing and the mixed background would be "
                   "empty. It is a knob, but it is not that kind of knob.");
    }
    c.mixing.donors_per_bin = static_cast<std::size_t>(donors);

    // "file_hash" is the only mode anything implements. Refusing the rest is the
    // same rule as beam.polarity: a config may not name a scheme that does not
    // exist and then run anyway on whatever the code happens to do instead.
    c.mixing.seed_mode = r.str("mixing.seed_mode");
    if (c.mixing.seed_mode != "file_hash") {
        r.fail("mixing.seed_mode",
               "is '" + c.mixing.seed_mode +
                   "', which nothing implements. The only supported mode is 'file_hash': the seed "
                   "is derived from the input file's identity, so (file, seed, donors_per_bin) "
                   "fully determines the pool and therefore every mixed pair.");
    }

    c.mixing.pool_grid.q2 = read_pool_axis(r, "mixing.pool_grid.q2_edges_gev2", "q2");
    c.mixing.pool_grid.xb = read_pool_axis(r, "mixing.pool_grid.xb_edges", "xb");
    c.mixing.pool_grid.n_photons_classes =
        read_mult_classes(r, "mixing.pool_grid.n_photons_classes");

    // PoolGridCuts owns what a usable pool grid is; this loader owns only the
    // mapping from JSON to it. Re-thrown through Reader::fail so that a bad grid
    // reports like every other bad key -- naming the file and the path -- and so
    // that Cuts::load()'s documented "throws std::runtime_error" holds.
    try {
        c.mixing.pool_grid.validate();
    } catch (const std::invalid_argument& e) {
        r.fail("mixing.pool_grid", std::string("is not a usable pool grid: ") + e.what());
    }

    // -- cross-checks against core/Constants.hpp ----------------------------
    require_agrees_with_constant(r, "beam.energy_gev", kBeamEnergyGeV, "kBeamEnergyGeV");
    require_agrees_with_constant(r, "beam.proton_mass_gev", kProtonMassGeV, "kProtonMassGeV");
    require_agrees_with_constant(r, "pairing.pi0_mass_gev", kPi0MassGeV, "kPi0MassGeV");

    // -- cross-checks on keys the struct deliberately does NOT carry --------
    //
    // A PDG code is an identity, not a threshold, so `Cuts` has no pid field
    // and the selection code names 11 and 22 as constants. Likewise the two
    // require_* flags describe cuts that this analysis always applies. Reading
    // them here and refusing to run on a mismatch keeps cuts.json honest: it
    // means no key in that file can quietly describe a selection the code does
    // not implement, which is the failure mode the old tree's dead TOML keys
    // are the monument to.
    const auto require_int_is = [&](const std::string& path, int expected, const std::string& why) {
        const int got = r.integer(path);
        if (got != expected) {
            r.fail(path, "is " + std::to_string(got) + ", but this analysis is built around " +
                             std::to_string(expected) + ". " + why);
        }
    };
    const auto require_true = [&](const std::string& path, const std::string& why) {
        if (!r.boolean(path)) {
            r.fail(path, "is false, but this analysis always applies that cut. " + why);
        }
    };

    require_int_is("electron.trigger.pid", 11,
                   "The trigger electron's PDG code is not a tunable threshold; it is not carried "
                   "in the Cuts struct.");
    require_int_is("photon.pid", 22,
                   "The photon's PDG code is not a tunable threshold; it is not carried in the "
                   "Cuts struct.");
    require_true("electron.trigger.require_negative_status",
                 "find_trigger_electron() requires status < 0 unconditionally.");
    require_true("photon.require_pcal_energy",
                 "The GBT pre-filter requires E_PCAL > 0 unconditionally.");

    return c;
}

}  // namespace pi0
