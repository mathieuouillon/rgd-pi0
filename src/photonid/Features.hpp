#pragma once

/// \file Features.hpp
/// \brief The 45-float input vector for the CLAS12 GBT photon classifier.
///
/// This is a PORT, not a redesign. It reproduces
///   /Users/mathieuouillon/Documents/tmp/clas-analysis-1/clas12/include/clas12/services/PhotonCutsService.hpp
/// (`filter_photon`, `get_calo_map`, `get_particle_calo_vector`,
/// `pass_pid_purity_filter`) feature-for-feature and bit-for-bit, reading
/// RDataFrame RVec columns instead of hipo::bank.
///
/// The five CatBoost models in photonid/models/ were trained on exactly these
/// 45 numbers in exactly this order. A change to the value, the order, or the
/// arithmetic precision of any one of them silently changes the classifier's
/// output -- it does not fail, it just gets the wrong answer. Every "why is
/// this written so awkwardly" in Features.cpp has that as its answer, and each
/// such spot carries a BIT-IDENTICAL comment.
///
/// Units: momenta and energies in GeV. Angles are in DEGREES at every boundary
/// of this header, per the project convention -- EXCEPT the feature vector
/// itself, whose elements 2 (gTheta) and the R_* neighbour angles are in
/// RADIANS because that is what the models were trained on. Those live inside
/// an opaque `std::vector<float>` and are never read back as an angle by this
/// project, so they cannot leak into a degree/radian mix-up at an API boundary.

#include <cstddef>
#include <map>
#include <vector>

#include <ROOT/RVec.hxx>

namespace pi0::photonid {

/// Length of the GBT feature vector. Structural, not a tunable: all five
/// headers in photonid/models/ declare `FloatFeatureCount = 45`.
inline constexpr std::size_t kNumFeatures = 45;

/// ECAL layer ids as they appear in REC::Calorimeter.layer.
inline constexpr int kLayerPcal = 1;
inline constexpr int kLayerEcin = 4;
inline constexpr int kLayerEcout = 7;

/// Calorimeter quantities for one particle in one ECAL layer.
///
/// Every field defaults to ZERO, and that zero is load-bearing: the classifier
/// distinguishes "this particle has no PCAL cluster" from "it has one" by
/// testing `pcal.x == 0` / `pcal.e > 0`, not by testing for the layer's
/// presence. See `Features.cpp::calo_position`.
struct CaloLayerData {
    double e{0.0};                  ///< deposited energy, GeV
    double x{0.0}, y{0.0}, z{0.0};  ///< cluster position, cm
    double m2u{0.0}, m2v{0.0};      ///< second moments of the shower profile
    double lu{0.0}, lv{0.0}, lw{0.0};  ///< fiducial distances along u/v/w, cm
    int sector{0};                  ///< 1..6, 0 if unset
};

/// The three ECAL layers for one particle.
///
/// NOTE: only `pcal` ever carries `e`, `m2u`, `m2v` in the feature vector --
/// the classifier was trained that way. `ecin`/`ecout` DO have those fields
/// populated here (the old code left them at zero), because the electron and
/// photon fiducial cuts downstream want them. That is safe: the feature builder
/// never reads `ecin.e` / `ecout.e` / their moments, so populating them cannot
/// move the score.
struct CaloRowData {
    CaloLayerData pcal{};
    CaloLayerData ecin{};
    CaloLayerData ecout{};
};

/// pindex -> per-layer calorimeter data, built once per event.
///
/// Replaces the old `std::map<int, CaloRowData> get_calo_map(const hipo::bank&)`.
class CaloMap {
   public:
    /// Build from the REC::Calorimeter columns of one event.
    ///
    /// All twelve RVecs must have the same length (one entry per calorimeter
    /// row); a mismatch throws std::invalid_argument rather than reading past
    /// the end of the short one.
    ///
    /// BIT-IDENTICAL: when two rows share a pindex AND a layer, the LAST row
    /// wins, exactly as the old code's unconditional `calo_map[pindex].pcal_x = x`
    /// assignment did. Layers other than PCAL/ECIN/ECOUT are ignored.
    ///
    /// BIT-IDENTICAL: the `detector` column is deliberately NOT a parameter and
    /// NOT filtered on. REC::Calorimeter only ever contains ECAL (detector 7),
    /// and the old code applied no detector cut; adding one here could only
    /// change behaviour, never fix it.
    [[nodiscard]] static CaloMap build(
        const ROOT::VecOps::RVec<short>& pindex,
        const ROOT::VecOps::RVec<short>& layer,
        const ROOT::VecOps::RVec<short>& sector,
        const ROOT::VecOps::RVec<float>& energy,
        const ROOT::VecOps::RVec<float>& x,
        const ROOT::VecOps::RVec<float>& y,
        const ROOT::VecOps::RVec<float>& z,
        const ROOT::VecOps::RVec<float>& m2u,
        const ROOT::VecOps::RVec<float>& m2v,
        const ROOT::VecOps::RVec<float>& lu,
        const ROOT::VecOps::RVec<float>& lv,
        const ROOT::VecOps::RVec<float>& lw);

    /// \return the calorimeter data for `pindex`, or nullptr if that particle
    ///         has NO calorimeter row at all.
    ///
    /// The nullptr case is not the same as "has rows but no PCAL": the old code
    /// rejected the first outright (`calo_map.find(row) == end() -> return false`)
    /// and let the second through to the PID-purity filter. Keep them distinct.
    [[nodiscard]] const CaloRowData* find(std::size_t pindex) const;

    /// \return true iff `pindex` has at least one calorimeter row.
    [[nodiscard]] bool contains(std::size_t pindex) const { return find(pindex) != nullptr; }

    /// \return the number of distinct particles with calorimeter data.
    [[nodiscard]] std::size_t size() const { return m_rows.size(); }

   private:
    std::map<int, CaloRowData> m_rows;
};

/// The cut values the feature builder needs. Populate from `Cuts::photon`:
///
///     pi0::photonid::FeatureCuts{cuts.photon.min_energy_gev,
///                                cuts.photon.theta_min_deg,
///                                cuts.photon.theta_max_deg}
///
/// Passed in rather than read from a Cuts here so that photonid/ does not
/// depend on config/, and so that these are testable without a cuts.json.
///
/// *** THESE THREE ARE NOT FREE PARAMETERS. *** They are baked into the trained
/// models: the GBTs learned the feature distributions produced by a 0.2 GeV /
/// 5 deg / 35 deg neighbour window. Retuning them in config/cuts.json will NOT
/// retune the classifier, it will feed the classifier a distribution it has
/// never seen and quietly degrade it. They live in cuts.json (and are read from
/// there rather than hard-coded) because the project forbids hard-coded cut
/// values, and because the same three numbers really are selection cuts when
/// applied to the candidate photon. Treat editing them as retraining work.
struct FeatureCuts {
    double min_energy_gev{};  ///< photon.min_energy_gev, GeV
    double theta_min_deg{};   ///< photon.theta_min_deg, DEGREES
    double theta_max_deg{};   ///< photon.theta_max_deg, DEGREES
};

/// The classifier's PID-purity pre-filter for the CANDIDATE photon.
///
/// Ports `filter_photon`'s early-outs: the candidate must have calorimeter
/// data at all, and must then satisfy
///     E >= min_energy_gev  &&  E_pcal > 0  &&  theta in [theta_min, theta_max]
/// where E = |p| (photons are massless) and theta comes from the MOMENTUM
/// vector, not the calorimeter position.
///
/// The caller must have already established that `cand_row` is a photon
/// (pid == 22); this function does not check, matching the old code.
///
/// Call this BEFORE build_features: the old `filter_photon` bailed out here and
/// never built a vector, so a candidate that fails must never be scored.
///
/// \param cand_row  row index into the REC::Particle columns.
[[nodiscard]] bool passes_pid_purity(
    std::size_t cand_row,
    const CaloMap& calo,
    const ROOT::VecOps::RVec<float>& px,
    const ROOT::VecOps::RVec<float>& py,
    const ROOT::VecOps::RVec<float>& pz,
    const FeatureCuts& cuts);

/// Build the 45-float GBT feature vector for the photon at `cand_row`.
///
/// Layout (indices are into the returned vector), reproducing
/// PhotonCutsService.hpp lines 440-465 exactly. Note the neighbour blocks are
/// GROUPED BY VARIABLE, not by neighbour -- all three R_gamma, then all three
/// dE_gamma, and so on:
///
///     [ 0] gE                 candidate energy = |p|, GeV
///     [ 1] gEpcal             candidate PCAL energy, GeV
///     [ 2] gTheta             candidate polar angle, RADIANS (see file docs)
///     [ 3] gm2u               candidate PCAL m2u
///     [ 4] gm2v               candidate PCAL m2v
///     [ 5] R_e                angle to the nearest electron, RADIANS
///     [ 6] dE_e               gE - E of that electron, GeV
///     [ 7.. 9] R_gamma    [0..2]   3 nearest photons, RADIANS
///     [10..12] dE_gamma   [0..2]
///     [13..15] Epcal_gamma[0..2]
///     [16..18] m2u_gamma  [0..2]
///     [19..21] m2v_gamma  [0..2]
///     [22..23] R_ch       [0..1]   2 nearest charged hadrons, RADIANS
///     [24..25] dE_ch      [0..1]
///     [26..27] Epcal_ch   [0..1]
///     [28..29] m2u_ch     [0..1]
///     [30..31] m2v_ch     [0..1]
///     [32..33] R_nh       [0..1]   2 nearest neutral hadrons, RADIANS
///     [34..35] dE_nh      [0..1]
///     [36..37] Epcal_nh   [0..1]
///     [38..39] m2u_nh     [0..1]
///     [40..41] m2v_nh     [0..1]
///     [42] num_photons_0_1     neighbour photons with R < 0.1  rad
///     [43] num_photons_0_2     neighbour photons with R < 0.2  rad
///     [44] num_photons_0_35    neighbour photons with R < 0.35 rad
///
/// Unfilled neighbour slots stay 0 -- and 0 is also how the insertion sort
/// spells "empty" (`R[i] == 0`), a quirk of the original preserved verbatim.
///
/// R is the 3D angle between CALORIMETER HIT POSITIONS, not between momenta.
/// dE is `gE - E_neighbour`.
///
/// \param cand_row  row index into the REC::Particle columns. MUST be a photon
///                  that already passed passes_pid_purity().
/// \param calo      this event's CaloMap.
/// \param pid,px,py,pz  the REC::Particle columns for this event.
/// \param cuts      the neighbour window (see FeatureCuts).
///
/// \throws std::invalid_argument if the columns disagree in length, if
///         `cand_row` is out of range, or if `cand_row` has no calorimeter data
///         (the old code returned false there; here it is a caller bug, because
///         passes_pid_purity() already reports that case).
///
/// \return exactly kNumFeatures floats.
///
/// *** DELIBERATE DEVIATION FROM THE REQUESTED SIGNATURE: no electron index. ***
/// The contract sketched `build_features(cand_row, calo, <REC_Particle cols>,
/// <electron index>)`. The old code takes NO electron index: features [5] and
/// [6] are the nearest of ANY row with pid == 11 in the event (see the
/// `else if (pid == 11)` branch), and REC::Particle routinely holds several.
/// Restricting the search to a caller-chosen scattered electron would change
/// the feature -- and hence the score -- on every multi-electron event, which
/// the "score must be bit-identical" requirement forbids. Faithfulness wins;
/// the parameter is therefore absent rather than accepted and ignored.
[[nodiscard]] std::vector<float> build_features(
    std::size_t cand_row,
    const CaloMap& calo,
    const ROOT::VecOps::RVec<int>& pid,
    const ROOT::VecOps::RVec<float>& px,
    const ROOT::VecOps::RVec<float>& py,
    const ROOT::VecOps::RVec<float>& pz,
    const FeatureCuts& cuts);

}  // namespace pi0::photonid
