#pragma once

/// \file Kinematics.hpp
/// \brief DIS and SIDIS kinematics. Pure functions, no ROOT, no HIPO, no state.
///
/// Units: momenta and energies in GeV, angles in DEGREES at every boundary.
/// The beam is assumed to travel along +z, and the INCIDENT electron to be
/// massless (E_beam = |k|, exact to a relative 1e-9 at 10.53 GeV). The
/// SCATTERED electron's energy is whatever the caller passes as `e_energy` --
/// this file makes no mass assumption about it. stageA_skim passes
/// sqrt(p^2 + m_e^2); see pi0::energy_from_p in core/Constants.hpp.
///
/// Every division in the implementation is guarded. Where the old analysis
/// would have produced a NaN (division by a zero nu, a zero |q|, an acos of
/// 1 + 1e-16), these functions return a defined, documented sentinel instead --
/// NaNs propagate silently into histograms and are a pain to trace back.

#include "core/Types.hpp"

namespace pi0::kin {

/// Compute the inclusive DIS kinematics of an event from the scattered electron.
///
/// \param ex,ey,ez   scattered-electron momentum components, GeV.
/// \param e_energy   scattered-electron energy E', GeV. Used directly in the
///                   Q2 and nu formulae (NOT recomputed from |p|), so that a
///                   caller supplying sqrt(p^2 + m_e^2) -- as stageA_skim does
///                   -- or a calorimeter-corrected energy gets the energy it
///                   asked for. This function assumes NOTHING about the
///                   scattered electron's mass; the choice is the caller's and
///                   is visible at the call site.
/// \param beam_energy  incident beam energy E_beam, GeV, along +z.
///
/// Definitions:
///   Q2 = 4 * E_beam * E' * sin^2(theta_e / 2)   -- the angle form, matching the
///                                                  analysis this replaces.
///
/// THE ANGLE FORM IS ITSELF THE MASSLESS ONE, and passing a massive E' does not
/// change that: Q2 = -(k - k')^2 = 2*(E_beam*E' - k.k') - m_e^2 - m_e^2 reduces
/// to 4*E_beam*E'*sin^2(theta/2) only when both electrons are massless. The
/// residual is O(m_e^2) ~ 2.6e-7 GeV^2, which is why cuts.json's /dis block
/// records the angle form as "equivalent to better than the resolution at these
/// energies" rather than as an identity. Carrying m_e in E' is a separate,
/// smaller correction (relative ~5e-9); neither is visible against a 1 GeV^2
/// cut, and both are recorded rather than hidden.
///   nu = E_beam - E'
///   xb = Q2 / (2 * M_p * nu)
///   W  = sqrt(M_p^2 + 2 * M_p * nu - Q2), 0 if the radicand is <= 0
///   y  = nu / E_beam
///
/// theta_e is the polar angle of the scattered electron w.r.t. the +z beam
/// direction, taken from the momentum vector (cos theta_e = ez / |p_e|).
///
/// WHY THE PROTON MASS, EVEN FOR NUCLEAR TARGETS: xb here is the PER-NUCLEON
/// Bjorken variable. RG-D measures the nuclear-to-deuteron ratio
/// R_A = (multiplicity on A) / (multiplicity on D) binned in these variables.
/// Using M_p for every target means the binning variable is defined identically
/// for A and for D, so the definition cancels in the ratio. Substituting a
/// nuclear mass M_A would rescale xb per target and destroy that cancellation --
/// the bins would no longer line up between numerator and denominator. This is
/// the standard convention in nuclear SIDIS (HERMES, CLAS, EIC), and it is a
/// choice of variable, not an approximation about the struck nucleon.
///
/// Guards: |p_e| <= 0 gives theta_e = 0 (hence Q2 = 0); nu <= 0 gives xb = 0;
/// beam_energy <= 0 gives y = 0; a non-positive W^2 gives W = 0.
[[nodiscard]] DisKin compute_dis(double ex, double ey, double ez, double e_energy, double beam_energy);

/// Compute the SIDIS kinematics of one hadron relative to the virtual photon.
///
/// \param hx,hy,hz   hadron momentum components, GeV.
/// \param he         hadron energy, GeV.
/// \param dis        the event DIS kinematics, for nu (from compute_dis).
/// \param ex,ey,ez   scattered-electron momentum components, GeV.
/// \param e_energy   scattered-electron energy, GeV. Accepted for symmetry with
///                   compute_dis and for a stable call signature; the gamma*
///                   frame is fixed by the three-momenta and nu alone, so this
///                   value is not read.
/// \param beam_energy  incident beam energy, GeV, along +z.
///
/// The gamma* frame:
///   q      = k_beam - k_e'          (three-vector)
///   q_hat  = q / |q|
///   y_hat  = (q x k_beam) / |q x k_beam|
///   x_hat  = y_hat x q_hat
/// (x_hat, y_hat, q_hat) is right-handed, and x_hat points along the component
/// of the INCOMING lepton transverse to q -- because
/// (q x k) x q = |q|^2 * (k - (k.q_hat) q_hat) -- which is what puts phi_h = 0
/// in the lepton scattering plane on the lepton side.
///
/// Then:
///   z      = E_h / nu
///   p_perp = p_h - (p_h . q_hat) q_hat
///   pt2    = |p_perp|^2
///   phi_h  = atan2(p_perp . y_hat, p_perp . x_hat), in DEGREES, [-180, 180]
///
/// THIS IS THE TRENTO CONVENTION (Bacchetta, D'Alesio, Diehl, Miller,
/// PRD 70, 117504 (2004)).
///
/// VERIFICATION (these numbers were measured, not asserted; both checks are
/// described completely enough below to re-derive from this comment alone,
/// beam along +z at kBeamEnergyGeV, theta_e in [2, 35] deg, hadron over 4pi):
///
/// 1. Cross-check against the canonical sign/arccos form
///
///      cos phi_h = (q_hat x k) . (q_hat x p_h) / (|q_hat x k| |q_hat x p_h|)
///      sin phi_h = ((k x p_h) . q_hat)         / (|q_hat x k| |q_hat x p_h|)
///      phi_h     = sign(sin phi_h) * arccos(cos phi_h)
///
///    Over 2e7 random configurations the two forms agree to a maximum absolute
///    deviation of 5e-12 degrees for |phi_h| in [1, 179].
///
/// 2. Round-trip against a constructed truth: build a hadron at a KNOWN phi_h
///    in the (x_hat, y_hat, q_hat) frame and check it is recovered. Over 2e7
///    configurations, max deviation 5e-13 degrees, uniform in phi_h. This test
///    does not involve the arccos form at all, so it bounds THIS function's
///    error rather than the reference's.
///
/// A WARNING ABOUT QUOTING A SINGLE AGREEMENT NUMBER HERE: taken over ALL
/// phi_h, the deviation in check 1 is NOT bounded -- it is 1.5e-8 degrees over
/// 2e6 configurations and 1.3e-7 over 2e7, growing with the sample simply
/// because a larger sample lands closer to phi_h = 0 and +-180. That residual
/// is the ARCCOS reference losing conditioning at its degenerate points (an
/// error eps in cos phi_h becomes an error ~sqrt(2 eps) in phi_h near 0), and
/// check 2 shows this function stays at ~5e-13 degrees there. So the atan2 form
/// implemented here is the better-conditioned of the two, and any all-phi
/// "agreement to X degrees" figure is a statement about the reference and the
/// sample size, not about this code. Do not tighten this comment into a single
/// number without re-reading this paragraph.
///
/// Guards: |q| <= 0 gives an all-zero result; nu <= 0 gives z = 0; an electron
/// exactly collinear with the beam leaves the lepton plane -- and therefore
/// phi_h -- undefined, and gives phi_h_deg = 0 with z and pt2 still valid.
[[nodiscard]] SidisKin compute_sidis(double hx, double hy, double hz, double he, const DisKin& dis, double ex, double ey, double ez,
                                     double e_energy, double beam_energy);

/// Opening angle between two three-vectors, in DEGREES, range [0, 180].
///
/// The acos argument is clamped to [-1, 1] before the call: for nearly parallel
/// vectors the dot/norm quotient can land at 1 + a few ulps and turn acos into
/// a NaN.
///
/// \return 0 if either vector has zero length (the angle is undefined; 0 is
///         returned rather than NaN so callers can cut on it safely).
[[nodiscard]] double angle_between_deg(double px1, double py1, double pz1, double px2, double py2, double pz2);

}  // namespace pi0::kin
