#include "stageB_bin/PoolGrid.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace pi0 {
namespace {

/// Validate one pool-grid axis. Grid1D does not validate itself (see its own
/// header: it is hot-path code, and Binning::load() validates the grids it
/// loads). This grid never goes through Binning::load(), so its axes are checked
/// here or nowhere.
void validate_axis(const Grid1D& g, const char* which) {
    if (g.edges.size() < 2) {
        throw std::invalid_argument("PoolGridCuts: the " + std::string(which) +
                                    " axis needs at least 2 edges (1 bin), got " +
                                    std::to_string(g.edges.size()) + ".");
    }
    for (std::size_t i = 0; i < g.edges.size(); ++i) {
        if (!std::isfinite(g.edges[i])) {
            throw std::invalid_argument("PoolGridCuts: the " + std::string(which) + " axis edge[" +
                                        std::to_string(i) +
                                        "] is not finite. An axis with a NaN or infinite edge has "
                                        "no bin anyone can name.");
        }
    }
    for (std::size_t i = 1; i < g.edges.size(); ++i) {
        if (!(g.edges[i] > g.edges[i - 1])) {
            std::ostringstream os;
            os.precision(17);
            os << "PoolGridCuts: the " << which << " axis edges must be STRICTLY increasing, but "
               << "edge[" << (i - 1) << "] = " << g.edges[i - 1] << " and edge[" << i << "] = "
               << g.edges[i]
               << ". Equal edges make a zero-width bin that can never be filled; decreasing edges "
                  "make a lookup that files photons into the wrong pool bin without ever erroring.";
            throw std::invalid_argument(os.str());
        }
    }
}

}  // namespace

bool MultClass::contains(std::size_t n_photons) const {
    // min is validated >= 1, so a zero-photon event matches nothing. The cast is
    // safe for any count a real event can carry, but a count above INT_MAX is
    // not a physics case -- it is a corrupt read -- and it must not wrap onto a
    // negative number that then compares below `min` and quietly matches
    // nothing. An open-ended class is the only honest home for it.
    if (n_photons > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        return !max.has_value();
    }
    const int n = static_cast<int>(n_photons);
    if (n < min) return false;
    return !max.has_value() || n <= *max;
}

int PoolGridCuts::mult_class(std::size_t n_photons) const {
    for (std::size_t i = 0; i < n_photons_classes.size(); ++i) {
        if (n_photons_classes[i].contains(n_photons)) return static_cast<int>(i);
    }
    return -1;
}

std::size_t PoolGridCuts::n_bins() const {
    if (q2.nbins() <= 0 || xb.nbins() <= 0) return 0;
    return static_cast<std::size_t>(q2.nbins()) * static_cast<std::size_t>(xb.nbins()) *
           n_photons_classes.size();
}

void PoolGridCuts::validate() const {
    validate_axis(q2, "Q^2");
    validate_axis(xb, "x_B");

    if (n_photons_classes.empty()) {
        throw std::invalid_argument(
            "PoolGridCuts: there are no photon-multiplicity classes, so no event can enter the "
            "pool and the mixed background would be empty.");
    }

    for (std::size_t i = 0; i < n_photons_classes.size(); ++i) {
        const MultClass& c = n_photons_classes[i];
        const std::string at = "class[" + std::to_string(i) + "] ('" + c.label + "')";

        if (c.min < 1) {
            throw std::invalid_argument(
                "PoolGridCuts: " + at + " has min = " + std::to_string(c.min) +
                ", but a multiplicity class must start at 1 or above: an event with zero photons "
                "produces no pairs and must enter no class.");
        }
        if (c.max.has_value() && *c.max < c.min) {
            throw std::invalid_argument("PoolGridCuts: " + at +
                                        " has max = " + std::to_string(*c.max) + " below min = " +
                                        std::to_string(c.min) + ", so it contains nothing.");
        }
        // An open-ended class absorbs every higher multiplicity, so anything
        // after it can never be reached. A config that lists one is a config
        // whose later classes are dead -- and a dead class that still counts
        // toward n_bins() silently allocates pool bins nothing can ever fill.
        if (!c.max.has_value() && i + 1 != n_photons_classes.size()) {
            throw std::invalid_argument(
                "PoolGridCuts: " + at +
                " is open-ended (max = null) but is not the last class. It absorbs every higher "
                "multiplicity, so every class after it is unreachable.");
        }
        if (i > 0) {
            const MultClass& prev = n_photons_classes[i - 1];
            // prev.max has a value: only the last class may be open-ended, and
            // prev is not the last one.
            if (c.min <= *prev.max) {
                throw std::invalid_argument(
                    "PoolGridCuts: " + at + " starts at " + std::to_string(c.min) +
                    ", which is not above the previous class's max of " + std::to_string(*prev.max) +
                    ". Classes must be ascending and disjoint, or mult_class() returns whichever "
                    "one happens to be listed first and the pool key stops meaning what it says.");
            }
        }
    }
}

}  // namespace pi0
