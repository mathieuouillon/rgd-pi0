#pragma once

/// \file ElectronSelection.hpp
/// \brief Finding the scattered electron and deciding whether it is a good one.
///
/// Units: momenta and energies in GeV, distances in cm, angles in DEGREES.

#include <array>
#include <cstddef>
#include <optional>
#include <string>

#include <ROOT/RVec.hxx>

#include "config/Cuts.hpp"

namespace pi0::selection {

/// PDG code of the electron.
///
/// This is an identity, not a threshold, which is why `Cuts` carries no pid
/// field: there is no configuration in which the scattered electron is not an
/// electron. cuts.json still records it (`electron.trigger.pid`) as
/// documentation, and Cuts::load() refuses to run if that key ever stops saying
/// 11 -- so this constant cannot drift away from the config unnoticed.
inline constexpr int kPdgElectron = 11;

/// Index of the trigger electron in REC::Particle, or nullopt if the event has
/// none.
///
/// The trigger electron is the FIRST row satisfying all of:
///   * pid == 11,
///   * status < 0,
///   * cuts.electron.status_min <= |status| < cuts.electron.status_max
///     (inclusive lower bound, EXCLUSIVE upper -- the Forward Detector band),
///   * a non-zero momentum vector (px, py, pz not all exactly zero).
///
/// "FIRST matching row", not "highest momentum". That is an ordering dependence
/// inherited from the old code and preserved deliberately: changing it would
/// change the selection, so it must be changed on purpose or not at all.
///
/// NOTE ON THE NAME: there is NO hardware trigger-bit cut in this analysis.
/// "Trigger electron" is a REC::Particle status definition -- the negative sign
/// and the 2000-4000 band -- and nothing is read from RUN::config::trigger. The
/// name is inherited; do not go looking for the trigger-bit test.
///
/// \param pid,status,px,py,pz  REC::Particle columns. The types are the ones
///        RHipoDS actually produces with n_inspect = 0 (every column is an
///        RVec, including single-row banks); pid is RVec<int>, status is
///        RVec<short>, and the momenta are RVec<float>.
/// \return the row index, usable to index any REC::Particle column.
///
/// Rows beyond the shortest supplied column are not examined: a truncated bank
/// yields no electron rather than an out-of-bounds read.
[[nodiscard]] std::optional<std::size_t> find_trigger_electron(
    const ROOT::VecOps::RVec<int>& pid,
    const ROOT::VecOps::RVec<short>& status,
    const ROOT::VecOps::RVec<float>& px,
    const ROOT::VecOps::RVec<float>& py,
    const ROOT::VecOps::RVec<float>& pz,
    const Cuts& cuts);

/// Stable identifiers for the electron cuts, in application order.
///
/// These are the values `ElectronCutResult::failed_at` takes. They are stage
/// IDENTIFIERS, not display labels: they never mention a number, so they cannot
/// go stale when a threshold moves. For something to show a human, call
/// electron_cutflow_label(), which renders the threshold from `Cuts`.
namespace electron_stage {
inline constexpr const char* kChi2Pid = "chi2pid";
inline constexpr const char* kMomentum = "momentum";
inline constexpr const char* kVertex = "vertex";
inline constexpr const char* kSamplingFraction = "sampling_fraction";
inline constexpr const char* kPcalFiducial = "pcal_fiducial";
inline constexpr const char* kDcEdge = "dc_edge";
}  // namespace electron_stage

/// The six stages, in the order pass_electron() applies them. A cutflow can be
/// built by iterating this and counting `failed_at`.
inline constexpr std::array<const char*, 6> kElectronStages = {
    electron_stage::kChi2Pid,          electron_stage::kMomentum,
    electron_stage::kVertex,           electron_stage::kSamplingFraction,
    electron_stage::kPcalFiducial,     electron_stage::kDcEdge,
};

/// The verdict on one electron.
struct ElectronCutResult {
    bool passed{};

    /// The stage that rejected it: one of the electron_stage constants, or
    /// nullptr iff `passed`. Because the cuts short-circuit, this is the FIRST
    /// failing stage, not the only one -- a cutflow built from it is a
    /// sequential cutflow, which is what it has always been.
    const char* failed_at{};
};

/// A human-readable cutflow row label for `stage`, with every threshold
/// rendered FROM `cuts`.
///
/// USE THIS RATHER THAN WRITING A LABEL BY HAND. The old cutflow printed the
/// row "Momentum > 0.8 GeV" while applying p > 2.0 -- the label was a stale
/// string literal that had drifted from the constant beside it, and it was
/// copied into a second algorithm and into every log anyone ever quoted. A
/// label that is computed from the value it describes cannot lie. If you find
/// yourself typing a number into a label, that is the bug reappearing.
///
/// \throws std::invalid_argument if `stage` is not one of kElectronStages.
[[nodiscard]] std::string electron_cutflow_label(const char* stage, const Cuts& cuts);

/// Apply the electron cuts, in order, short-circuiting on the first failure.
///
/// The order is inherited from the old ElectronSelectionAlg and is preserved:
///   1. chi2pid    strictly inside (chi2pid_min, chi2pid_max)
///   2. momentum   p > min_momentum_gev
///   3. vertex     the caller's verdict (see `vertex_passed`)
///   4. sampling fraction, per sector and polarity
///   5. PCAL fiducial   lv > 9.0 AND lw > 9.0
///   6. DC edge         R1 > 1.68 AND R2 > 2.0 AND R3 > 8.75
///
/// ALL comparisons are STRICT. Every threshold comes from `cuts`.
///
/// \param chi2pid            REC::Particle.chi2pid for this row.
/// \param p_gev              |p| of the electron, GeV.
/// \param vertex_passed      THE CALLER'S VERTEX VERDICT. The vertex cut is not
///        implemented here and must not be: it is target-dependent, it needs the
///        per-cell correction parameterisation, and for the solid foils it is
///        what DEFINES which nucleus the event belongs to. It lives in
///        pi0::vertex. Duplicating even its window here would recreate the
///        two-tables-that-disagree problem the old tree had between its pi0 path
///        and its electron_analysis binary. Pass the answer in.
/// \param sampling_fraction  (E_PCAL + E_ECIN + E_ECOUT) / p, dimensionless.
/// \param pcal_sector        PCAL sector, 1-6. Outside that range the sampling
///        fraction stage fails (see selection::pass), so the electron is
///        rejected AT that stage -- which is the honest answer for a track the
///        calorimeter cannot place.
/// \param pcal_lv_cm,pcal_lw_cm  REC::Calorimeter lv/lw for the PCAL cluster, cm.
/// \param dc_edge_r1_cm,dc_edge_r2_cm,dc_edge_r3_cm  REC::Traj `edge` at DC
///        layers 6, 18 and 36 respectively, cm.
/// \param cuts               the loaded configuration. Not optional, by design.
[[nodiscard]] ElectronCutResult pass_electron(double chi2pid,
                                              double p_gev,
                                              bool vertex_passed,
                                              double sampling_fraction,
                                              int pcal_sector,
                                              double pcal_lv_cm,
                                              double pcal_lw_cm,
                                              double dc_edge_r1_cm,
                                              double dc_edge_r2_cm,
                                              double dc_edge_r3_cm,
                                              const Cuts& cuts);

}  // namespace pi0::selection
