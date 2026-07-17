#pragma once

/// \file Pairing.hpp
/// \brief Building gamma-gamma pairs out of a photon list. No ROOT, no HIPO.
///
/// Units: momenta and energies in GeV, angles in DEGREES.

#include <optional>
#include <vector>

#include "core/Types.hpp"

namespace pi0 {

/// Cut parameters for gamma-gamma pairing. Populated from cuts.json.
struct PairingCuts {
    /// Half-width of the accepted |m_gg - m_pi0| window, GeV. Deliberately WIDE:
    /// the signal is extracted by fitting the m_gg spectrum downstream, so this
    /// window must keep enough sideband to constrain the background shape.
    double mass_window_gev{};

    /// Absolute floor on m_gg, GeV. Kills the split-cluster pile-up at m_gg ~ 0.
    double min_mgg_gev{};

    /// Minimum angle between a photon and the scattered electron, DEGREES.
    /// NOT applied by find_gg_pairs -- see the note on that function.
    double e_gamma_min_angle_deg{};

    /// Minimum-opening-angle curve: a * exp(-b * p) + offset.
    double open_a_deg{};      ///< a, DEGREES
    double open_b_inv_gev{};  ///< b, GeV^-1 (the exponential's decay constant)
    double open_offset_deg{}; ///< offset, DEGREES (the high-momentum asymptote)
};

/// Momentum-dependent minimum opening angle between the two photons, DEGREES.
///
///   theta_min(p) = a * exp(-b * p) + offset
///
/// The shape tracks the calorimeter's two-cluster resolving power: a high-|p|
/// pi0 decays into a narrow pair, so the cut has to relax with momentum, but it
/// cannot relax past `offset`, the point where two clusters stop being
/// separable at all.
///
/// \param pair_momentum_gev  |p_1 + p_2| of the pair, GeV.
[[nodiscard]] double min_opening_angle_deg(double pair_momentum_gev, const PairingCuts& cuts);

/// Test ONE photon pair against the admissibility criteria, and combine it.
///
/// This is find_gg_pairs()'s inner predicate, given a name and made public. A
/// pair is ADMISSIBLE when all three hold:
///   * m_gg > cuts.min_mgg_gev
///   * theta_gg > min_opening_angle_deg(|p_a + p_b|, cuts)
///   * |m_gg - kPi0MassGeV| < cuts.mass_window_gev
///
/// WHY THIS IS PUBLIC, AND WHY THAT IS NOT SCOPE CREEP
/// ---------------------------------------------------
/// find_gg_pairs() takes ONE photon list and pairs it against itself, greedily
/// and EXCLUSIVELY. The event-mixing background (stageB_bin) needs the opposite
/// shape: photons from list A (this event's) against photons from list B (the
/// frozen donor pool's), EXHAUSTIVELY and non-exclusively -- every current
/// photon against every donor, because a donor is a sample of a distribution
/// rather than a physical object that can be consumed, and because the mixed
/// spectrum's job is to trace a shape with as much statistics as it can get.
/// (cuts.json's /mixing/_comment records that the superseded analysis was
/// exhaustive here too, so this is not a change of method.)
///
/// So find_gg_pairs() cannot serve the mixed pass. The alternative was for
/// stageB_bin to re-type these three tests locally -- a SECOND COPY of the
/// admissibility rule, which would agree with this one until the day somebody
/// edited one of them, and then the same-event and mixed spectra would be cut
/// differently while still looking exactly like spectra. That is the defect this
/// project spends more words warning about than any other, and the mixed
/// background is the worst possible place for it: a background shape cut
/// slightly differently from its signal subtracts badly and reports nothing.
///
/// Extracting the predicate makes "same-event and mixed use the SAME cuts"
/// structural instead of a comment. find_gg_pairs() calls this and nothing else
/// tests admissibility, so there is exactly one copy and the existing pairing
/// tests pin it.
///
/// \return the combined pair if admissible, `std::nullopt` otherwise.
///         `GGPair::i` and `::j` are set to 0 and INDEX NOTHING -- this function
///         has no list to index into. find_gg_pairs() overwrites them with the
///         real indices; a mixed-pair caller has no meaningful value for them
///         and must not read them.
[[nodiscard]] std::optional<GGPair> admissible_pair(const Photon& a, const Photon& b, const PairingCuts& cuts);

/// Build gamma-gamma pairs from a photon list by GREEDY EXCLUSIVE pairing.
///
/// Algorithm (matching the old analysis's greedy structure):
///   repeat:
///     among all pairs (i < j) with both photons still unused, a pair is
///     ADMISSIBLE when all three hold:
///       * m_gg > min_mgg_gev
///       * theta_gg > min_opening_angle_deg(|p_i + p_j|)
///       * |m_gg - kPi0MassGeV| < mass_window_gev
///     take the admissible pair minimising |m_gg - kPi0MassGeV|, emit it, and
///     mark both photons used;
///   until no admissible pair remains.
///
/// Each photon is consumed by at most one pair, so N photons yield at most
/// floor(N/2) pairs. Ties in |m_gg - m_pi0| are broken deterministically in
/// favour of the lowest (i, j) in lexicographic order, so the output does not
/// depend on iteration accidents.
///
/// DELIBERATE PHYSICS CHANGE vs. THE OLD ANALYSIS -- A/B TEST THIS:
/// the old analysis applied the opening-angle and mass-floor cuts DOWNSTREAM of
/// pairing. A pair that those cuts later rejected had still CONSUMED both of its
/// photons, so those photons could never pair with anything else. Testing the
/// cuts inside the admissibility check, as here, frees those photons to form
/// other pairs. The two algorithms therefore do not merely differ by a
/// selection: this one can find pi0 candidates the old one structurally could
/// not. It is not a bug fix, and it should be measured rather than assumed.
///
/// A worked case where the two genuinely differ (verified, not hypothetical),
/// with mass_window = 0.06, min_mgg = 0.02, (a, b, offset) = (12, 0.6, 2):
///   photon 0: E = 6.0   GeV at theta = 20.000 deg, phi = 0
///   photon 1: E = 6.0   GeV at theta = 21.289 deg, phi = 0
///   photon 2: E = 0.3   GeV at theta = 25.000 deg, phi = 0
/// Pair (0,1) has m_gg = 0.1350 -- dead on the pi0 mass -- but theta_gg = 1.289
/// deg against a theta_min of 2.009: a real pi0 too energetic for the
/// calorimeter to resolve. Pair (0,2) has m_gg = 0.1170 and theta_gg = 5.0 deg,
/// comfortably admissible. The OLD algorithm picks (0,1) first on mass, consumes
/// photons 0 and 1, then throws the pair away downstream on the opening angle --
/// orphaning photon 2 and yielding ZERO pairs. This one rejects (0,1) at the
/// admissibility test and yields ONE pair, (0,2).
///
/// Note what that example implies: the freed photons get re-paired with whatever
/// is left, which here means a genuine pi0's decay photon is recombined with an
/// unrelated soft photon into a combinatorial "candidate" at m_gg = 0.117. So
/// this change does not simply recover efficiency -- it also feeds the m_gg fit
/// a background it never used to see. Which of the two effects wins is an
/// empirical question about the fit, and is exactly why this is flagged for an
/// A/B test rather than presented as an improvement.
///
/// THE e_gamma CUT IS NOT APPLIED HERE. It is a property of a single photon,
/// not of a pair, so it belongs upstream: filter the photon list on
/// `Photon::e_gamma_deg > cuts.e_gamma_min_angle_deg` before calling this.
/// `PairingCuts::e_gamma_min_angle_deg` is carried in this struct only so that
/// the pairing block of cuts.json maps onto one struct.
///
/// m_gg is computed from the massless four-vector sum:
///   m_gg^2 = (E_1 + E_2)^2 - |p_1 + p_2|^2 = 2 E_1 E_2 (1 - cos theta_12),
/// with the tiny negative m_gg^2 that rounding produces for a nearly collinear
/// pair clamped to 0.
///
/// \param photons  the photon list; `GGPair::i` and `::j` index into it.
/// \return the emitted pairs, in the order the greedy loop found them (i.e.
///         best |m_gg - m_pi0| first).
[[nodiscard]] std::vector<GGPair> find_gg_pairs(const std::vector<Photon>& photons, const PairingCuts& cuts);

}  // namespace pi0
