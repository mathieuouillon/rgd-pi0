#pragma once

/// \file Constants.hpp
/// \brief The single source of truth for physics constants used across the
///        analysis.
///
/// Every constant lives here exactly once. Do NOT re-define any of these in a
/// translation unit, a program, or a config file: a duplicated constant that
/// drifts out of sync is the classic way an analysis silently changes physics.
///
/// Unit conventions for the whole project:
///   * energies and momenta are in GeV,
///   * angles are in DEGREES at every public API boundary (a parameter that
///     wants radians must be named `*_rad` and say so in its doc comment).

#define _USE_MATH_DEFINES  // MSVC needs this before <cmath> to expose M_PI.
#include <cmath>

#ifndef M_PI
/// Fallback for toolchains that do not expose the POSIX `M_PI` from <cmath>.
#define M_PI 3.14159265358979323846
#endif

namespace pi0 {

/// Nominal CLAS12 RG-D beam energy, GeV.
inline constexpr double kBeamEnergyGeV = 10.53;

// ---------------------------------------------------------------------------
// Particle masses, GeV.
//
// ONE definition of each, here. Two consumers need them and neither may keep a
// private copy: pi0::kin (the scattered electron's energy) and pi0::photonid
// (the neighbour energies E = sqrt(p^2 + m^2) that feed the GBT's dE features).
// The photonid port originally stashed four of these in an anonymous namespace
// in Features.cpp because this header did not own them; that is exactly the
// duplication this header's preamble warns about, so they moved here.
// ---------------------------------------------------------------------------

/// Proton mass, GeV. Used for the DIS variables even on nuclear targets
/// (see pi0::kin::compute_dis for why).
inline constexpr double kProtonMassGeV = 0.938272;

/// Neutron mass, GeV.
inline constexpr double kNeutronMassGeV = 0.93956536;

/// Neutral pion mass, GeV (PDG).
inline constexpr double kPi0MassGeV = 0.1349768;

/// Electron mass, GeV.
///
/// The value the superseded tree used (clas12/include/clas12/Constants.hpp
/// ELECTRONMASS), kept bit-identical for the same reason kPionMassGeV is: it
/// enters the GBT's feature vector.
inline constexpr double kElectronMassGeV = 0.00051099891;

/// Charged-pion mass, GeV.
///
/// BIT-IDENTICAL WART -- DO NOT "CORRECT" THIS. 0.1396 is NOT the PDG
/// charged-pion mass (0.13957039). It is what the old Constants::PIMASS said,
/// it enters E = sqrt(p^2 + m^2) for every pi+/pi- neighbour of a photon, and
/// therefore enters the dE_ch features the CatBoost models were TRAINED on.
/// Changing it does not fix a number, it feeds the classifier a feature vector
/// that means something different from what the model expects.
inline constexpr double kPionMassGeV = 0.1396;  // sic -- see above

/// Charged-kaon mass, GeV.
inline constexpr double kKaonMassGeV = 0.493677;

/// Relativistic energy of a particle of mass `mass_gev` and momentum `p_gev`:
/// E = sqrt(p^2 + m^2). GeV throughout.
///
/// A function rather than an expression repeated at each site, so that "which
/// mass did this energy use" is answered by reading one argument. Photons are
/// MASSLESS by convention in this analysis (E_gamma = |p_gamma|, see
/// config/cuts.json's /photon block) -- do not route them through this with a
/// zero mass to look symmetric; sqrt(p^2) is not |p| for free.
[[nodiscard]] inline double energy_from_p(double p_gev, double mass_gev) {
    return std::sqrt(p_gev * p_gev + mass_gev * mass_gev);
}

/// Multiply radians by this to get degrees.
inline constexpr double kRadToDeg = 180.0 / M_PI;

/// Multiply degrees by this to get radians.
inline constexpr double kDegToRad = M_PI / 180.0;

}  // namespace pi0
