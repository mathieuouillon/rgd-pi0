#pragma once

/// \file PhotonSelection.hpp
/// \brief The photon selection, applied once at skim time.
///
/// Photons are identified by a Catboost GBT classifier, not by cut-based PID.
/// This file owns the CUTS around that classifier; it does not own the
/// classifier. The score is computed in pi0::photonid, from the verbatim
/// Catboost model headers, and arrives here as a number.
///
/// WARNING, inherited from the old design and recorded in cuts.json: this
/// selection is FROZEN at skim time. The skim carries no PID, status, beta or
/// vertex column for photons, so nothing downstream can re-cut or vary any of
/// it. Varying a value in the photon block as a systematic means re-running the
/// skim over the full dataset.
///
/// Units: momenta and energies in GeV, distances in cm, angles in DEGREES.

#include <cmath>

#include "config/Cuts.hpp"
#include "core/Constants.hpp"

namespace pi0::selection {

/// PDG code of the photon. An identity, not a threshold -- see kPdgElectron in
/// ElectronSelection.hpp for why this is a constant and not a `Cuts` field.
/// Cuts::load() verifies cuts.json's `photon.pid` still agrees.
inline constexpr int kPdgPhoton = 22;

/// Polar angle of a photon, DEGREES.
///
/// \return NaN for an identically-zero momentum vector, which has no direction.
///         Every caller here tests the energy first, and E = |p| = 0 fails the
///         energy floor, so the NaN cannot reach a decision -- and if it did, it
///         would compare false against both theta bounds and reject the photon,
///         which is the safe direction.
[[nodiscard]] inline double photon_theta_deg(double px, double py, double pz) {
    const double p = std::sqrt(px * px + py * py + pz * pz);
    if (p == 0.0) return std::nan("");
    return std::acos(pz / p) * kRadToDeg;
}

/// The classifier's purity pre-filter: the cheap gate that a photon must pass
/// before the GBT is asked about it.
///
///     E_gamma >= min_energy_gev     INCLUSIVE  (unlike most cuts here)
///     E_PCAL  >  0                  the shower must start in the preshower
///     theta_min_deg <= theta <= theta_max_deg   INCLUSIVE at both ends
///
/// Exposed separately because it is worth applying before evaluating the model:
/// the GBT is the expensive part of the skim, and this filter is three
/// comparisons.
///
/// \param theta_deg  polar angle in DEGREES. The old code's equivalent took the
///        parameter as `double theta`, in radians, and converted inside -- an
///        unmarked radian at a public boundary, which is one of the two
///        degree/radian bugs this project exists to not repeat. If you have
///        radians, convert at the call site or use photon_theta_deg().
[[nodiscard]] bool passes_gbt_prefilter(double e_gamma_gev,
                                        double pcal_energy_gev,
                                        double theta_deg,
                                        const Cuts& cuts);

/// The full photon selection, with the classifier score supplied.
///
/// Order, inherited from the old PhotonSelectionAlg cutflow:
///   1. pid == 22
///   2. the GBT pre-filter (energy, PCAL energy, theta)
///   3. gbt_score > gbt_threshold
///   4. PCAL fiducial   lv > 14.0 AND lw > 14.0
///   5. beta strictly inside (0.9, 1.1)
///
/// \param px,py,pz  photon momentum, GeV. The energy is E = |p|: photons are
///        massless and no calorimeter energy correction is applied anywhere in
///        this chain (the note records the fitted peak sitting ~4.4% low as a
///        consequence).
/// \param gbt_score THE SIGMOID-MAPPED CLASSIFIER PROBABILITY, in [0, 1] -- i.e.
///        1/(1 + exp(-model_output)), not the raw model output. cuts.json's
///        gbt_threshold (0.78) is a probability, so handing this the raw score
///        would silently apply a different cut. pi0::photonid does the mapping.
/// \param beta      REC::Particle.beta.
/// \param cuts      the loaded configuration.
[[nodiscard]] bool pass_photon(int pid,
                               double px,
                               double py,
                               double pz,
                               double pcal_energy_gev,
                               double pcal_lv_cm,
                               double pcal_lw_cm,
                               double beta,
                               double gbt_score,
                               const Cuts& cuts);

/// As pass_photon(), but the score is produced by a callable that is invoked
/// ONLY if the pid check and the pre-filter have already passed.
///
/// This is the form the skim should use: the pre-filter exists precisely so the
/// GBT is not evaluated on the majority of clusters, and the eager overload
/// above throws that away by requiring the score up front.
///
/// \param score_fn  any callable returning double, e.g.
///        `[&] { return photonid::score(run, row, ...); }`. Same contract as
///        `gbt_score`: it must return the sigmoid-mapped probability.
template <typename ScoreFn>
[[nodiscard]] bool pass_photon_scored(int pid,
                                      double px,
                                      double py,
                                      double pz,
                                      double pcal_energy_gev,
                                      double pcal_lv_cm,
                                      double pcal_lw_cm,
                                      double beta,
                                      ScoreFn&& score_fn,
                                      const Cuts& cuts) {
    if (pid != kPdgPhoton) return false;

    const double e_gamma_gev = std::sqrt(px * px + py * py + pz * pz);  // massless: E = |p|
    if (!passes_gbt_prefilter(e_gamma_gev, pcal_energy_gev, photon_theta_deg(px, py, pz), cuts)) {
        return false;
    }

    if (!(static_cast<double>(score_fn()) > cuts.photon.gbt_threshold)) return false;

    // 14.0 cm, harder than the electron's 9.0 cm. Intentional -- do not unify.
    if (!(pcal_lv_cm > cuts.photon.pcal_lv_min_cm && pcal_lw_cm > cuts.photon.pcal_lw_min_cm)) {
        return false;
    }

    // Strict at both ends.
    return beta > cuts.photon.beta_min && beta < cuts.photon.beta_max;
}

}  // namespace pi0::selection
