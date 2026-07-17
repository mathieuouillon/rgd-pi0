#pragma once

/// \file VzCorrector.hpp
/// \brief Vertex-z correction for the RG-D solid targets, plus the per-target
///        acceptance window. Port of clas-analysis-1's VzCorrector.hpp +
///        VzCorrectorService.hpp. No ROOT, no HIPO.
///
/// Units: lengths in cm, momenta in GeV/c, angles in DEGREES.
///
/// ---------------------------------------------------------------------------
/// WHAT THIS IS
/// ---------------------------------------------------------------------------
/// The RG-D solid-target runs put two foils in the beam. Reconstructed v_z
/// smears with (sector, p, phi, theta) in a way that blurs the two foil peaks
/// into each other. The correction, fitted offline in Python and exported to
/// data/Vz/vz_corrector_params.txt, standardises each event so that foil events
/// land on a common reference peak with a common width, at which point a flat
/// n-sigma window separates the foils.
///
/// For each event we look up cubic-in-theta coefficients from the cell
/// (sector, p-bin, phi-bin) of the active variant, evaluate (mu, sigma) for
/// BOTH peaks, form the two z-scores, standardise each to its reference
/// (mu0, sigma0), then sigmoid-blend the two hypotheses by the posterior
/// P(lo | vz) under a two-Gaussian mixture with equal priors.
///
/// ---------------------------------------------------------------------------
/// *** NAMING TRAP -- READ THIS BEFORE TOUCHING THE lo/hi PARAMETERS ***
/// ---------------------------------------------------------------------------
/// 'lo' and 'hi' name the PARAMETER SET, not the z ordering. They are just the
/// two columns of the fit, and the fit happened to call the DOWNSTREAM foil
/// 'lo'. Concretely, for outbending:
///
///     lo peak: mu = -2.916 cm   <- the LESS negative, DOWNSTREAM foil  -> Sn
///     hi peak: mu = -7.861 cm   <- the MORE negative, UPSTREAM   foil  -> Cu
///
/// So mu_hi is the more negative number. 'hi' does NOT mean "higher z".
///
/// The old VzCorrectorService.hpp doc-comment got this backwards -- it said
///     "Cu : corrected vz inside the hi-peak (downstream) window"
///     "Sn : corrected vz inside the lo-peak (upstream) window"
/// The ROUTING there was right (Cu->hi, Sn->lo); the parentheticals were
/// exactly inverted. The hi peak at -7.861 is UPSTREAM, the lo peak at -2.916
/// is DOWNSTREAM. This port keeps the routing and fixes the words.
///
/// lo/hi survive HERE, in `Params`, because the params file's columns are
/// literally named that way and `correct()` must read them. They do NOT survive
/// into the CUT: config/cuts.json names targets explicitly and lists each
/// target's peaks, precisely so the trap is unreachable from the selection side:
///     /vertex/targets/cu/peaks = [{mu_cm: -7.861, sigma_cm: 0.415}]   (the hi peak)
///     /vertex/targets/sn/peaks = [{mu_cm: -2.916, sigma_cm: 0.370}]   (the lo peak)
/// Two independent sources -- the old service's routing and a config written
/// from the analysis note that never saw the old header -- agree that Cu is the
/// -7.9 foil. tests/test_vertex.cpp pins the acceptance semantics and
/// tests/test_selection.cpp pins the shipped numbers, so neither can regress
/// silently.
///
/// ---------------------------------------------------------------------------
/// DEGREES, EXPLICITLY
/// ---------------------------------------------------------------------------
/// phi_bin_index() TAKES DEGREES. Every angular parameter in this header is
/// named *_deg and is in degrees; there is no radian quantity anywhere in this
/// file, so there is nothing named *_rad. The sector centres below are degrees.
/// This is called out because the superseded analysis shipped a sector helper
/// documented as radians and fed degrees. The unit is in the parameter name,
/// and tests/test_vertex.cpp asserts the degree convention directly by feeding
/// sector centres in degrees and demanding the centre bin back -- a radian
/// implementation would fail it.

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace pi0::vertex {

/// RG-D target. LD2 is the liquid-deuterium cryotarget; the rest are the solid
/// foils. CxC is the two-carbon-foil configuration (both foils are carbon, so
/// either peak is signal); Cu and Sn are the two foils of the CuSn
/// configuration, resolved one at a time by the peak window.
enum class Target { LD2, CxC, Cu, Sn };

/// Torus polarity. This project is OUTBENDING-only for now, but the params file
/// carries the inbending variants and they stay loadable.
enum class Polarity { Inbending, Outbending };

/// Human-readable tag, for error messages and logs.
[[nodiscard]] std::string to_string(Target t);
[[nodiscard]] std::string to_string(Polarity p);

/// One reference peak of the CORRECTED v_z distribution, cm.
///
/// These are the Voigt-fit Gaussian-equivalent peaks measured on the
/// already-corrected distribution -- which is what the cut must be drawn
/// against, and is NOT the same thing as the (mu0, sigma0) reference peaks in
/// `Params` that the correction standardises TO. The two disagree for
/// cusn_outbending's hi peak: the params file says sigma = 0.3851, the code that
/// produced the published numbers used 0.415 (+8%), and config/cuts.json
/// independently records 0.415. See `VertexCuts` below and open question Q2 in
/// the note.
struct VzPeak {
    double mu_cm{};     ///< Peak centre, cm.
    double sigma_cm{};  ///< Peak width, cm.
};

/// How a target's vertex cut is drawn. Mirrors config/cuts.json's
/// /vertex/targets/<target>/rule.
enum class VzRule {
    RawWindow,       ///< Flat (min, max) on the RAW v_z, exclusive. LD2 only.
    CorrectedPeaks,  ///< |v_z_corr - mu| < n_sigma * sigma for ANY listed peak.
};

/// The vertex cut for ONE target, as carried by config/cuts.json's
/// /vertex/targets/<target> block. There is no default: pi0::Cuts::load() writes
/// every field or throws.
///
/// *** WHICH PEAKS A TARGET SELECTS IS DATA, NOT CODE. ***
/// The old code routed by a switch on the target (Cu -> hi, Sn -> lo, CxC ->
/// either), which is where the lo/hi naming trap could bite. Here the routing IS
/// the `peaks` list: Cu's list holds only the -7.9 peak, Sn's only the -2.9 one,
/// and CxC's holds both because both of its foils are carbon. Acceptance is
/// "inside ANY listed peak", so the same rule expresses all three.
struct VzTargetCuts {
    /// Is the vz correction applied before the cut? False for LD2 (one long
    /// cryotarget cell, no foil peaks to separate), true for the solid foils.
    /// Tied to `rule`: RawWindow cuts the RAW v_z, CorrectedPeaks the corrected
    /// one, and pi0::Cuts::load() refuses a config that mixes them.
    bool correction_enabled{};

    VzRule rule{VzRule::RawWindow};

    double vz_min_cm{};  ///< RawWindow only. Exclusive.
    double vz_max_cm{};  ///< RawWindow only. Exclusive.

    double n_sigma{};           ///< CorrectedPeaks only. Window half-width, in sigmas.
    std::vector<VzPeak> peaks;  ///< CorrectedPeaks only. Cu/Sn carry one, CxC two.
};

/// config/cuts.json's whole /vertex block: the acceptance window for every RG-D
/// target.
///
/// ---------------------------------------------------------------------------
/// THESE OVERRIDE THE PARAMS FILE'S OWN REFERENCE PEAKS. Q2 in the note.
/// ---------------------------------------------------------------------------
/// `Params`' (mu0, sigma0) are the reference peaks the correction standardises
/// TO. The cut is drawn against the Voigt-fit peaks measured on the corrected
/// distribution, which mostly agree with them -- but not everywhere:
///
///     cusn_outbending sigma_hi:  params file says 0.3851, cuts.json says 0.415  (+8%)
///
/// cuts.json's value is the one that produced the published numbers. The split
/// is deliberate and must be preserved: the CORRECTION POLYNOMIALS come from
/// data/Vz/vz_corrector_params.txt, the CUT WINDOW comes from config/cuts.json.
/// The discrepancy itself is unexplained and is recorded, not resolved.
///
/// ---------------------------------------------------------------------------
/// OUTBENDING ONLY
/// ---------------------------------------------------------------------------
/// cuts.json's /vertex block is the outbending set, as its beam block says of
/// every polarity-dependent number in the file. There is no inbending window
/// here and no hard-coded table to fall back to -- pi0::Cuts::load() refuses an
/// inbending config outright rather than apply these numbers to it. Adding
/// inbending means adding a per-polarity block to cuts.json, not editing these.
struct VertexCuts {
    VzTargetCuts ld2{};
    VzTargetCuts cxc{};
    VzTargetCuts cu{};
    VzTargetCuts sn{};

    /// The block for one target.
    [[nodiscard]] const VzTargetCuts& for_target(Target t) const;
};

/// Multi-variant vz correction + per-target acceptance window.
///
/// Load once with `load()`, select a variant with `set_variant()`, then call
/// `correct()` per event. `pass_window()` is static and needs no loaded file.
///
/// Thread safety: `load()` and `set_variant()` mutate; `correct()` and
/// `pass_window()` are const/static and safe to call concurrently once the
/// variant is set.
class VzCorrector {
   public:
    // ---- geometry / schema constants ------------------------------------
    static constexpr int kNSectors = 6;  ///< CLAS12 sectors.
    static constexpr int kNPhi     = 6;  ///< Phi bins per sector.
    static constexpr int kNPoly    = 4;  ///< Cubic: 4 coefficients, lowest order first.

    /// Sector phi centres in DEGREES, indexed by `sector - 1`.
    static constexpr std::array<double, kNSectors> kSectorPhiCenterDeg = {
        0.0, 60.0, 120.0, 180.0, -120.0, -60.0,
    };

    static constexpr double kPhiHalfWidthDeg = 30.0;                              ///< Sector half-width, deg.
    static constexpr double kPhiBinWidthDeg  = 2.0 * kPhiHalfWidthDeg / kNPhi;    ///< 10 deg.

    /// Fallback theta clip domain for cells whose fitted domain is NaN (the
    /// sentinel the exporter writes for empty cells). Mirrors TH_LO / TH_HI in
    /// vz_correction.py. NOT a cut: it is the domain the cubic was fitted over,
    /// so it belongs to the correction, not to the selection.
    static constexpr double kThetaLoFallbackDeg = 6.0;
    static constexpr double kThetaHiFallbackDeg = 30.0;

    /// Floor on the per-cell predicted sigma. A numerical guard against a cubic
    /// that dips to <= 0 inside its domain and would divide by ~0. Not a cut.
    static constexpr double kSigmaFloor = 0.01;

    // NOTE: there is deliberately NO window table here -- no kNSigma, no LD2
    // range, no per-target peaks. Every one of those is a CUT VALUE and lives in
    // config/cuts.json, which this class reads through `VertexCuts`. A cut
    // spelled as a constant in this header is a bug: see config/README.md.

    /// Per-variant parameter set: 6 sectors x n_p p-bins x 6 phi-bins cells of
    /// cubic-in-theta coefficients, plus the reference peaks and the p-binning.
    struct Params {
        double mu0_lo{};        ///< Reference lo-peak centre post-correction, cm.
        double sigma0_lo{1.0};  ///< Reference lo-peak width  post-correction, cm.
        double mu0_hi{};        ///< Reference hi-peak centre post-correction, cm.
        double sigma0_hi{1.0};  ///< Reference hi-peak width  post-correction, cm.
        double p_lo{};          ///< Lower edge of the first p-bin, GeV/c.
        double p_binw{1.0};     ///< p-bin width, GeV/c.
        int    n_p{};           ///< Number of p-bins.

        /// Cubic coefficients per cell, lowest order first.
        std::vector<std::array<double, kNPoly>> mu_lo;
        std::vector<std::array<double, kNPoly>> mu_hi;
        std::vector<std::array<double, kNPoly>> sig_lo;
        std::vector<std::array<double, kNPoly>> sig_hi;
        /// Per-cell [theta_min, theta_max] clip domain in DEGREES. NaN entries
        /// fall back to [kThetaLoFallbackDeg, kThetaHiFallbackDeg].
        std::vector<std::array<double, 2>> theta_dom;

        /// Flat cell index: (sector-1) * n_p * kNPhi + ip * kNPhi + iphi.
        /// Matches the row order the exporter writes.
        [[nodiscard]] std::size_t cell_index(int sector, int ip, int iphi) const {
            return static_cast<std::size_t>((sector - 1) * n_p * kNPhi + ip * kNPhi + iphi);
        }

        /// Number of cells: kNSectors * n_p * kNPhi. 576 for the shipped n_p=16.
        [[nodiscard]] std::size_t n_cells() const {
            return static_cast<std::size_t>(kNSectors) * static_cast<std::size_t>(n_p) *
                   static_cast<std::size_t>(kNPhi);
        }
    };

    VzCorrector() = default;

    /// Load every variant from data/Vz/vz_corrector_params.txt.
    ///
    /// Wire format (whitespace-separated; '#' lines ignored):
    ///     n_variants n_phi n_sectors n_poly
    ///     variant_key mu0_lo sigma0_lo mu0_hi sigma0_hi p_lo p_binw n_p
    ///     sector ip iphi  c_mu_lo[4] c_mu_hi[4] c_sig_lo[4] c_sig_hi[4]  th_min th_max
    ///     ... (n_sectors * n_p * n_phi rows per variant)
    ///
    /// FAILS LOUDLY (std::runtime_error) on: unopenable file, missing/malformed
    /// global header, schema mismatch against kNPhi/kNSectors/kNPoly, n_p <= 0,
    /// duplicate variant key, truncated block, malformed row, an out-of-range
    /// (sector, ip, iphi) triplet, or any cell left unwritten / written twice.
    ///
    /// The last two are STRICTER than the old loader, deliberately. The old one
    /// took (sector, ip, iphi) from the row and indexed with it WITHOUT any
    /// range check -- a typo'd sector in the file was an out-of-bounds vector
    /// write, i.e. silent memory corruption rather than an error. For a
    /// well-formed file the behaviour is identical.
    [[nodiscard]] static VzCorrector load(const std::string& path);

    /// Select the active variant. Cu/Sn share the "cusn" fit; CxC has its own.
    /// \throws std::runtime_error for LD2 (no variant exists -- LD2 is
    ///         uncorrected by construction) or if the variant is not loaded.
    void set_variant(Target target, Polarity polarity);

    /// Variant key for a target/polarity, e.g. "cusn_outbending".
    /// `std::nullopt` for LD2, which has no correction.
    [[nodiscard]] static std::optional<std::string> variant_key(Target target, Polarity polarity);

    /// Map track phi to a phi-bin index within its sector.
    ///
    /// *** TAKES DEGREES. *** Wraps (phi_deg - sector_centre_deg) into
    /// [-180, 180), clips to [-30, +30) deg, then bins into kNPhi bins of
    /// kPhiBinWidthDeg. A track at its sector's centre lands in bin 3.
    ///
    /// \param phi_deg Track phi in the lab frame, DEGREES. Any wrapping.
    /// \param sector  CLAS12 sector, [1, 6]. Out-of-range values are clamped.
    /// \return        Bin index in [0, kNPhi).
    [[nodiscard]] static int phi_bin_index(double phi_deg, int sector);

    /// Apply the correction to one event. Requires `set_variant()` first.
    ///
    /// Out-of-range inputs are clamped, not rejected: the p-bin index to
    /// [0, n_p-1] (with the shipped p_lo=2, p_binw=0.5, n_p=16 this is p
    /// clamped to [2, 10] GeV/c), theta to the cell's fitted domain, and sector
    /// to [1, 6].
    ///
    /// \param vz        Raw vertex z, cm.
    /// \param p         Track momentum, GeV/c.
    /// \param theta_deg Track polar angle, DEGREES.
    /// \param phi_deg   Track azimuthal angle in the lab frame, DEGREES.
    /// \param sector    CLAS12 sector, [1, 6].
    /// \return          Corrected vz, cm.
    /// \throws std::runtime_error if no variant is active.
    [[nodiscard]] double correct(double vz, double p, double theta_deg, double phi_deg, int sector) const;

    /// Does this event pass the vertex cut described by `tc`?
    ///
    /// The two rules, both straight from config/cuts.json:
    ///   * RawWindow      -> `vz_cm` is the RAW vz, tested against the flat
    ///                       range (vz_min_cm, vz_max_cm), EXCLUSIVE at both
    ///                       ends (matching the old service's
    ///                       `vz > min && vz < max`). This is LD2.
    ///   * CorrectedPeaks -> `vz_cm` is the CORRECTED vz from `correct()`,
    ///                       accepted if |vz_cm - mu| < n_sigma * sigma for ANY
    ///                       peak in `tc.peaks`.
    ///
    /// For the outbending config that means Cu accepts only the -7.9 (upstream)
    /// foil, Sn only the -2.9 (downstream) one, and CxC either -- but that
    /// routing is `tc.peaks`, not a switch here. For Cu and Sn this is not
    /// merely a quality cut: it IS the target assignment, since the two foils
    /// differ only in v_z.
    ///
    /// \param vz_cm  Corrected vz for CorrectedPeaks, RAW vz for RawWindow.
    [[nodiscard]] static bool pass_window(double vz_cm, const VzTargetCuts& tc);

    /// `pass_window` for one target of a loaded /vertex block. Equivalent to
    /// `pass_window(vz_cm, v.for_target(target))`.
    [[nodiscard]] static bool pass_window(double vz_cm, Target target, const VertexCuts& v);

    /// Read-only access to a loaded variant's parameters. For tests and
    /// diagnostics.
    /// \throws std::runtime_error if that variant is not loaded.
    [[nodiscard]] const Params& params_for(Target target, Polarity polarity) const;

    /// Is a variant present in the loaded file?
    [[nodiscard]] bool has_variant(Target target, Polarity polarity) const;

    /// Every loaded variant key, sorted. `unordered_map` order is not stable,
    /// so this sorts to keep callers and tests deterministic.
    [[nodiscard]] std::vector<std::string> variant_keys() const;

    /// Key of the active variant, empty if none is selected.
    [[nodiscard]] const std::string& active_variant() const { return m_active_key; }

   private:
    [[nodiscard]] const Params& active() const;

    std::unordered_map<std::string, Params> m_variants;
    const Params* m_active = nullptr;
    std::string   m_active_key;
};

}  // namespace pi0::vertex
