#pragma once

/// \file Cuts.hpp
/// \brief The runtime image of config/cuts.json: every cut value the analysis
///        applies, loaded once, from one file.
///
/// THE RULE THIS TYPE EXISTS TO ENFORCE: no cut value may be hard-coded. A
/// threshold spelled as a literal in a translation unit is a bug, because the
/// analysis this replaces scattered its cuts across C++ constants, six TOML
/// files and a set of dead config keys -- and the TOML keys were declared, set,
/// and never read, so anyone who "tuned" them changed nothing. Everything a cut
/// depends on comes through this struct or it does not exist.
///
/// Cuts::load() FAILS LOUDLY. There is no default for any key: a missing key
/// throws rather than silently selecting a plausible number, because a silent
/// default is exactly how a cutflow starts describing a selection that never ran.
///
/// Units: energies and momenta in GeV, angles in DEGREES, distances in cm.

#include <string>

#include "core/Pairing.hpp"
#include "stageB_bin/PoolGrid.hpp"
#include "vertex/VzCorrector.hpp"

namespace pi0 {

/// Every cut value in the analysis, as loaded from config/cuts.json.
///
/// The nested anonymous structs mirror the JSON's top-level blocks one-for-one.
/// Field names are NOT always the JSON key path (e.g. `electron.sf_n_sigma`
/// comes from `electron.sampling_fraction.n_sigma`); Cuts::load() owns that
/// mapping and is the only place it appears.
///
/// The members are value-initialised so that a default-constructed `Cuts` is
/// deterministically zero rather than indeterminate. That is a guard against
/// undefined behaviour, NOT a set of defaults: load() writes every field or
/// throws, so a zero can never reach the physics through the supported path.
struct Cuts {
    struct {
        double energy_gev{};    ///< beam energy, GeV. Pinned to pi0::kBeamEnergyGeV at load.
        std::string polarity{}; ///< "inbending" or "outbending"; anything else throws.
    } beam;

    struct {
        double chi2pid_min{}, chi2pid_max{};  ///< open interval, strict at both ends
        double min_momentum_gev{};            ///< p > this. It is 2.0. It has ALWAYS been 2.0.
        double sf_n_sigma{};                  ///< from electron.sampling_fraction.n_sigma
        double pcal_lv_min_cm{}, pcal_lw_min_cm{};  ///< 9.0 -- NOT the photon's 14.0
        double dc_edge_r1_cm{}, dc_edge_r2_cm{}, dc_edge_r3_cm{};  ///< DC layers 6 / 18 / 36
        int status_min{}, status_max{};       ///< |status| window, [min, max): inclusive lower, EXCLUSIVE upper
    } electron;

    struct {
        double q2_min{};  ///< GeV^2, strict
        double w_min{};   ///< GeV, strict
        double y_max{};   ///< strict. No q2_max/w_max: neither is reachable at this beam energy.
    } dis;

    struct {
        double gbt_threshold{};   ///< sigmoid(model output) > this
        double min_energy_gev{};  ///< E_gamma >= this. INCLUSIVE, unlike most cuts here.
        double theta_min_deg{}, theta_max_deg{};   ///< INCLUSIVE at both ends, DEGREES
        double pcal_lv_min_cm{}, pcal_lw_min_cm{}; ///< 14.0 -- NOT the electron's 9.0. Do not unify.
        double beta_min{}, beta_max{};             ///< open interval, strict at both ends
        int gbt_pass{};                            ///< cooking pass, keys the model lookup with the run number
        bool allow_rga_fallback{};                 ///< DEFAULT FALSE. See cuts.json before touching.
    } photon;

    /// The gamma-gamma pairing block. Defined in core/Pairing.hpp so that the
    /// pairing code needs no dependency on this header.
    PairingCuts pairing{};

    /// SIDIS-level cuts on a reconstructed pi0, from cuts.json's /pairing block.
    ///
    /// NOT part of PairingCuts, on purpose: find_gg_pairs() cannot apply these
    /// and must not appear to. z = E_pi0 / nu needs the event's nu, which the
    /// pairing code has never seen and has no business seeing. The JSON keeps
    /// the keys under /pairing for historical reasons and load() owns the
    /// mapping, exactly as it does for electron.sf_n_sigma.
    ///
    /// PROVENANCE: pairing.z_min and pairing.z_max sat in cuts.json DECLARED AND
    /// UNREAD until 2026-07-16 -- the same "declared, set, never read" defect
    /// this file's preamble holds up as the thing to avoid, reproduced inside
    /// the file meant to end it. Nothing had needed them yet; make_grid did, and
    /// that is what found it. The lesson is that the rule is enforced by
    /// consumers, not by intent.
    struct {
        double z_min{};  ///< z > this. STRICT (cuts.json: "EXCLUSIVE at both ends")
        double z_max{};  ///< z < this. STRICT.
    } sidis;

    /// The SHAPE of the two factorized equal-statistics grids -- bin counts per
    /// axis, from cuts.json's /binning block.
    ///
    /// THE EDGES ARE NOT HERE AND WILL NEVER BE HERE. They are a measurement
    /// over a dataset, computed by src/tools/make_grid and written to
    /// config/binning/grid_A_q2_xb.json and grid_B_z_pt2.json with their own
    /// provenance block. Stage B loads those files and hashes them into its
    /// output. Two reasons for the split, and both are the point of the rewrite:
    ///
    ///   * an edge array retyped into this file by hand is an edge array that
    ///     drifts from the data that defined it, with nothing to catch it;
    ///   * the superseded analysis's kd-tree edges were never archived anywhere
    ///     except /work, so the binning of its production is unrecoverable if
    ///     that disk is lost. Edges belong in version control, next to the
    ///     record of which files and which config produced them.
    ///
    /// So: the shape is a human's design choice and lives here; the edges are
    /// data and live in their own version-controlled file.
    ///
    /// "Equal statistics" is a statement about each axis's MARGINAL -- see
    /// cuts.json's /binning/_equal_statistics_is_marginal_comment. The 56 cells
    /// of Grid A do NOT hold equal counts, and nothing in this project pretends
    /// they do.
    struct {
        int n_q2{};   ///< Grid A, Q^2 axis. Default 8.
        int n_xb{};   ///< Grid A, x_B axis. Default 7.
        int n_z{};    ///< Grid B, z axis. Default 5.
        int n_pt2{};  ///< Grid B, p_T^2 axis. Default 5.
    } binning;

    /// The m_gg histogram axis, from cuts.json's /mgg_histogram block: 200
    /// uniform bins over [0.0, 0.3] GeV as shipped.
    ///
    /// A CUT VALUE, NOT A PLOTTING PREFERENCE, and that block's own _comment
    /// says so: the sideband normalisation alpha and the +-3 sigma yield window
    /// are both sums over these bins, so changing the binning changes the yield.
    /// It is shared by the same-event and mixed-event spectra -- they are
    /// subtracted from each other bin by bin, which is only meaningful on one
    /// axis -- and by the per-(4D bin, m_gg bin) kinematic sums, so that the
    /// abscissa can be sideband-subtracted exactly like the yield.
    ///
    /// PROVENANCE: this block sat in cuts.json DECLARED AND NEVER READ until
    /// stageB_bin needed it. That is the FOURTH time this file's founding defect
    /// has been found inside the file written to end it (see Cuts::sidis for the
    /// second and Cuts::mixing for the third). Nothing had needed it yet; the
    /// rule is enforced by consumers, not by intent.
    struct {
        double min_gev{};  ///< low edge of the first bin
        double max_gev{};  ///< top edge of the last bin
        int bins{};        ///< uniform bin count over [min_gev, max_gev]
    } mgg_histogram;

    /// The phi_h axis of the beam-spin-asymmetry histogram, from /bsa.
    ///
    /// DEGREES, Trento convention, matching SidisKin::phi_h_deg -- which is
    /// atan2-derived and therefore lands in [-180, 180]. An axis over [0, 360]
    /// would put every negative phi_h out of range and silently discard half the
    /// data, so the range is not a free choice: it is fixed by the producer.
    ///
    /// There is no helicity key here. helicity == 0 means UNDEFINED (see
    /// stageA_skim's helicity.convention provenance) and is not a state to
    /// configure around; stageB_bin fills the BSA for helicity == +-1 only, and
    /// helicity == 0 events still enter every other output. See stageB_bin's
    /// header for why that asymmetry is correct rather than an omission.
    struct {
        int n_phi_bins{};       ///< uniform bins over [phi_min_deg, phi_max_deg]. 12 as shipped.
        double phi_min_deg{};   ///< -180
        double phi_max_deg{};   ///< +180
    } bsa;

    /// The event-mixing block, from /mixing. Defined in stageB_bin/PoolGrid.hpp
    /// for the same reason PairingCuts lives in core/Pairing.hpp and VertexCuts
    /// in vertex/VzCorrector.hpp: the mixing code must not depend on this header,
    /// and this way `Cuts` composes it rather than either owning the other.
    ///
    /// PROVENANCE: the whole /mixing block -- donors_per_bin, seed_mode, and all
    /// 224 pool bins' worth of grid -- sat in cuts.json DECLARED AND NEVER READ
    /// until DonorPool needed it. That is the third time this file's own founding
    /// defect has been found inside the file written to end it (see Cuts::sidis
    /// for the second). The rule is enforced by CONSUMERS, not by intent: a key
    /// nothing reads is dead however good the reason it was written.
    ///
    /// NOTE mixing.pool_grid's Q^2/x_B axes are NOT Cuts::binning's Grid A. They
    /// are the hand-authored configuration product edges, shared with the N_DIS
    /// normalisation denominator; Grid A is machine-fitted equal-statistics edges
    /// in config/binning/grid_A_q2_xb.json. They hold the same numbers today only
    /// because Grid A's file is still a placeholder. Do not unify them, and do not
    /// let one stand in for the other -- see PoolGrid.hpp and that file's own
    /// _coincidence_warning.
    MixingCuts mixing{};

    /// The per-target vertex windows, from /vertex/targets. Defined in
    /// vertex/VzCorrector.hpp for the same reason PairingCuts lives in
    /// core/Pairing.hpp: the vertex code must not depend on this header, and
    /// this way `Cuts` composes it rather than either owning the other.
    ///
    /// For Cu and Sn this block is not merely a quality cut -- it IS the target
    /// assignment, since the two foils are distinguished only by v_z.
    ///
    /// OUTBENDING ONLY. cuts.json's /vertex block is the outbending set and
    /// there is no hard-coded table behind it, so load() refuses an inbending
    /// config rather than silently apply these numbers to it.
    pi0::vertex::VertexCuts vertex{};

    /// Parse `json_path` and return a fully populated Cuts.
    ///
    /// \throws std::runtime_error if the file cannot be opened, is not valid
    ///         JSON, is missing ANY key this loader reads, has a key of the
    ///         wrong JSON type, carries an unknown polarity, or disagrees with
    ///         core/Constants.hpp about a physics constant. There is no partial
    ///         success and no defaulting: you get every value or an exception.
    [[nodiscard]] static Cuts load(const std::string& json_path);
};

}  // namespace pi0
