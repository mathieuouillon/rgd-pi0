#include "photonid/PhotonGBT.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace pi0::photonid {

double score(const std::vector<float>& feats, ModelFn model) {
    if (model == nullptr) {
        throw std::invalid_argument("pi0::photonid::score: null model function");
    }
    if (feats.size() != kNumFeatures) {
        throw std::invalid_argument(
            "pi0::photonid::score: got " + std::to_string(feats.size()) + " features, the model needs "
            "exactly " + std::to_string(kNumFeatures));
    }

    const double margin = model(feats);
    return 1.0 / (1.0 + std::exp(-margin));
}

bool passes(double photon_score, double threshold) {
    return photon_score > threshold;
}

}  // namespace pi0::photonid
