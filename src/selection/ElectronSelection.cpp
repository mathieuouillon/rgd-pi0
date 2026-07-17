#include "selection/ElectronSelection.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include "selection/SamplingFraction.hpp"

namespace pi0::selection {
namespace {

/// Render a threshold the way a human wants to read it: "2" and not "2.000000".
/// Default float formatting, which drops trailing zeros and keeps enough digits
/// to be unambiguous for the values in cuts.json.
[[nodiscard]] std::string fmt(double value) {
    std::ostringstream os;
    os << value;
    return os.str();
}

[[nodiscard]] bool stage_is(const char* stage, const char* which) {
    // Compare by value, not by pointer: a caller may well have round-tripped
    // the identifier through a std::string on its way into a cutflow map.
    return stage != nullptr && std::strcmp(stage, which) == 0;
}

}  // namespace

std::optional<std::size_t> find_trigger_electron(const ROOT::VecOps::RVec<int>& pid,
                                                 const ROOT::VecOps::RVec<short>& status,
                                                 const ROOT::VecOps::RVec<float>& px,
                                                 const ROOT::VecOps::RVec<float>& py,
                                                 const ROOT::VecOps::RVec<float>& pz,
                                                 const Cuts& cuts) {
    // Never read past the shortest column. In a well-formed REC::Particle these
    // are all the same length; if they are not, the bank is broken and the
    // right answer is "no electron", not an out-of-bounds read.
    const std::size_t n =
        std::min({pid.size(), status.size(), px.size(), py.size(), pz.size()});

    for (std::size_t i = 0; i < n; ++i) {
        if (pid[i] != kPdgElectron) continue;

        const int st = static_cast<int>(status[i]);
        if (st >= 0) continue;  // must be negative: the row is trigger-flagged

        const int abs_status = std::abs(st);
        if (abs_status < cuts.electron.status_min) continue;
        if (abs_status >= cuts.electron.status_max) continue;  // EXCLUSIVE upper bound

        // "Non-zero momentum" means the vector is not identically zero. Exact
        // comparison against 0 is deliberate and is not a tolerance question:
        // it filters rows the reconstruction left unfilled, whose components
        // are exactly 0.0f, not merely small.
        if (px[i] == 0.0f && py[i] == 0.0f && pz[i] == 0.0f) continue;

        return i;  // FIRST match wins. Ordering dependence, preserved on purpose.
    }
    return std::nullopt;
}

int electron_stage_index(const char* stage) {
    if (stage == nullptr) return -1;  // ElectronCutResult::passed -- no stage rejected it

    for (std::size_t i = 0; i < kElectronStages.size(); ++i) {
        if (stage_is(stage, kElectronStages[i])) return static_cast<int>(i);
    }

    throw std::invalid_argument(
        std::string("pi0::selection::electron_stage_index: unknown stage '") + stage + "'");
}

std::string electron_cutflow_label(const char* stage, const Cuts& cuts) {
    const auto& e = cuts.electron;

    if (stage_is(stage, electron_stage::kChi2Pid)) {
        return "chi2pid in (" + fmt(e.chi2pid_min) + ", " + fmt(e.chi2pid_max) + ")";
    }
    if (stage_is(stage, electron_stage::kMomentum)) {
        return "Momentum > " + fmt(e.min_momentum_gev) + " GeV";
    }
    if (stage_is(stage, electron_stage::kVertex)) {
        // No number: the window is target-dependent and lives in pi0::vertex.
        // Better to say nothing than to say a number this module cannot know.
        return "Vertex z (target window)";
    }
    if (stage_is(stage, electron_stage::kSamplingFraction)) {
        return "Sampling fraction within " + fmt(e.sf_n_sigma) + " sigma";
    }
    if (stage_is(stage, electron_stage::kPcalFiducial)) {
        return "PCAL fiducial lv > " + fmt(e.pcal_lv_min_cm) + " cm, lw > " +
               fmt(e.pcal_lw_min_cm) + " cm";
    }
    if (stage_is(stage, electron_stage::kDcEdge)) {
        return "DC edge R1 > " + fmt(e.dc_edge_r1_cm) + " cm, R2 > " + fmt(e.dc_edge_r2_cm) +
               " cm, R3 > " + fmt(e.dc_edge_r3_cm) + " cm";
    }

    throw std::invalid_argument(
        std::string("pi0::selection::electron_cutflow_label: unknown stage '") +
        (stage != nullptr ? stage : "(null)") + "'");
}

ElectronCutResult pass_electron(double chi2pid,
                                double p_gev,
                                bool vertex_passed,
                                double sampling_fraction,
                                int pcal_sector,
                                double pcal_lv_cm,
                                double pcal_lw_cm,
                                double dc_edge_r1_cm,
                                double dc_edge_r2_cm,
                                double dc_edge_r3_cm,
                                const Cuts& cuts) {
    const auto& e = cuts.electron;

    // 1. chi2pid -- strictly inside the open interval.
    if (!(e.chi2pid_min < chi2pid && chi2pid < e.chi2pid_max)) {
        return {false, electron_stage::kChi2Pid};
    }

    // 2. momentum. This is 2.0 GeV. The old cutflow's "> 0.8 GeV" row label was
    //    a stale string; the applied cut was always this one.
    if (!(p_gev > e.min_momentum_gev)) {
        return {false, electron_stage::kMomentum};
    }

    // 3. vertex -- the caller's verdict, computed in pi0::vertex.
    if (!vertex_passed) {
        return {false, electron_stage::kVertex};
    }

    // 4. sampling fraction. `pass` returns false for a sector outside [1, 6],
    //    so an unplaceable track is rejected here rather than throwing.
    if (!pass(sampling_fraction, p_gev, pcal_sector, polarity_from_string(cuts.beam.polarity),
              e.sf_n_sigma)) {
        return {false, electron_stage::kSamplingFraction};
    }

    // 5. PCAL fiducial. 9.0 cm -- the LOOSE level. Photons use 14.0 cm and the
    //    difference is intentional (cuts.json says so twice); do not unify.
    if (!(pcal_lv_cm > e.pcal_lv_min_cm && pcal_lw_cm > e.pcal_lw_min_cm)) {
        return {false, electron_stage::kPcalFiducial};
    }

    // 6. DC edge distance at regions 1/2/3 (REC::Traj layers 6/18/36).
    if (!(dc_edge_r1_cm > e.dc_edge_r1_cm && dc_edge_r2_cm > e.dc_edge_r2_cm &&
          dc_edge_r3_cm > e.dc_edge_r3_cm)) {
        return {false, electron_stage::kDcEdge};
    }

    return {true, nullptr};
}

}  // namespace pi0::selection
