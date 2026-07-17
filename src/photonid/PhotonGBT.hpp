#pragma once

/// \file PhotonGBT.hpp
/// \brief Turning a GBT model's raw margin into a photon probability, and that
///        probability into a yes/no.
///
/// Ports `classify_photon` from the old clas-analysis-1 PhotonCutsService,
/// split into its two halves so that the score can be histogrammed, cached and
/// tested without re-deciding the cut, and so that the threshold is visibly a
/// parameter rather than a constant buried next to the arithmetic.

#include <vector>

#include "photonid/Features.hpp"          // kNumFeatures
#include "photonid/RunRangeModelMap.hpp"  // ModelFn

namespace pi0::photonid {

/// Photon probability from the classifier.
///
///     score = 1 / (1 + exp(-model(feats)))
///
/// The CatBoost applicators return a raw log-odds margin, not a probability;
/// the logistic is what the old code applied on top, and it is what the
/// threshold in cuts.json is calibrated against.
///
/// \param feats  exactly kNumFeatures (45) floats from build_features().
/// \param model  from select_model().
///
/// \throws std::invalid_argument if `model` is null or `feats.size()` is not
///         kNumFeatures. The generated applicators loop `for i in 0..45` over
///         `features[i]` with no bounds check, so a short vector is a silent
///         heap overread producing a plausible-looking score -- exactly the
///         failure this check exists to prevent.
///
/// \return a probability in [0, 1].
[[nodiscard]] double score(const std::vector<float>& feats, ModelFn model);

/// \param photon_score  from score().
/// \param threshold     Cuts::photon.gbt_threshold (config/cuts.json
///                      photon.gbt_threshold). NEVER hard-code it here: the
///                      old code's `static constexpr double AI_THRESHOLD = 0.78`
///                      is precisely the kind of scattered constant cuts.json
///                      exists to abolish.
///
/// \return true iff the candidate is classified as a photon.
///
/// BIT-IDENTICAL: strictly greater, matching `prediction > m_ai_threshold`.
/// `>=` would differ for a score landing exactly on the threshold.
[[nodiscard]] bool passes(double photon_score, double threshold);

}  // namespace pi0::photonid
