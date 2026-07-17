#include "photonid/RunRangeModelMap.hpp"

#include <stdexcept>
#include <string>

// The generated CatBoost applicators. VERBATIM copies from clas-framework --
// never edit them. Each is `inline`, so including all five here costs one
// translation unit's compile time and nothing at link time.
#include "photonid/models/RGA_inbending_pass1.hpp"
#include "photonid/models/RGA_inbending_pass2.hpp"
#include "photonid/models/RGA_outbending_pass1.hpp"
#include "photonid/models/RGA_outbending_pass2.hpp"
#include "photonid/models/RGC_Summer2022_pass1.hpp"

namespace pi0::photonid {

namespace models {

// Each ApplyCatboostModel_* name is overloaded (one taking just the float
// vector, one also taking an ignored vector<string>). The static_cast to
// ModelFn is what picks the first overload; without it the address is
// ambiguous.
ModelFn rga_inbending_pass1() { return static_cast<ModelFn>(&ApplyCatboostModel_RGA_inbending_pass1); }
ModelFn rga_inbending_pass2() { return static_cast<ModelFn>(&ApplyCatboostModel_RGA_inbending_pass2); }
ModelFn rga_outbending_pass1() { return static_cast<ModelFn>(&ApplyCatboostModel_RGA_outbending_pass1); }
ModelFn rga_outbending_pass2() { return static_cast<ModelFn>(&ApplyCatboostModel_RGA_outbending_pass2); }
ModelFn rgc_summer2022_pass1() { return static_cast<ModelFn>(&ApplyCatboostModel_RGC_Summer2022_pass1); }

}  // namespace models

namespace {

struct ModelRange {
    int run_min;
    int run_max;
    int pass;
    ModelFn (*fn)();
    std::string_view name;
};

// ---------------------------------------------------------------------------
// The run -> model map.
//
// TRANSCRIBED ENTRY-FOR-ENTRY from PhotonCutsService.hpp lines 173-190 (16
// entries, runs 5032..16772), in the same source order. The old container was a
// std::map keyed by {run_min, run_max, pass} and searched linearly; the ranges
// are disjoint, so a flat table searched linearly is the same lookup.
//
// Reading the table: the run ranges are NOT contiguous and NOT sorted (6616
// precedes 6156, exactly as upstream), and the gaps between them are real --
// runs 5667..6155, 6604..6615, 6784..11092, 11301..11322, 11572..16041 have no
// trained model. So does everything after 16772, which includes all of RG-D.
//
// Note the last two rows: RG-C Summer 2022 pass 2 is deliberately mapped to the
// pass1 model, because no pass2 RG-C model was ever trained. That is upstream's
// choice, preserved.
// ---------------------------------------------------------------------------
constexpr ModelRange kModelTable[] = {
    {5032, 5332, 1, &models::rga_inbending_pass1, "RGA_inbending_pass1"},
    {5032, 5332, 2, &models::rga_inbending_pass2, "RGA_inbending_pass2"},
    {5333, 5666, 1, &models::rga_outbending_pass1, "RGA_outbending_pass1"},
    {5333, 5666, 2, &models::rga_outbending_pass2, "RGA_outbending_pass2"},
    {6616, 6783, 1, &models::rga_inbending_pass1, "RGA_inbending_pass1"},
    {6616, 6783, 2, &models::rga_inbending_pass2, "RGA_inbending_pass2"},
    {6156, 6603, 1, &models::rga_inbending_pass1, "RGA_inbending_pass1"},
    {6156, 6603, 2, &models::rga_inbending_pass2, "RGA_inbending_pass2"},
    {11093, 11283, 1, &models::rga_outbending_pass1, "RGA_outbending_pass1"},
    {11093, 11283, 2, &models::rga_outbending_pass2, "RGA_outbending_pass2"},
    {11284, 11300, 1, &models::rga_inbending_pass1, "RGA_inbending_pass1"},
    {11284, 11300, 2, &models::rga_inbending_pass2, "RGA_inbending_pass2"},
    {11323, 11571, 1, &models::rga_inbending_pass1, "RGA_inbending_pass1"},
    {11323, 11571, 2, &models::rga_inbending_pass2, "RGA_inbending_pass2"},
    {16042, 16772, 1, &models::rgc_summer2022_pass1, "RGC_Summer2022_pass1"},
    {16042, 16772, 2, &models::rgc_summer2022_pass1, "RGC_Summer2022_pass1"},
};

static_assert(sizeof(kModelTable) / sizeof(kModelTable[0]) == 16,
              "the old s_model_map had exactly 16 entries; if this count changed, the port drifted");

const ModelRange* lookup(int run, int pass) {
    for (const auto& entry : kModelTable) {
        if (run >= entry.run_min && run <= entry.run_max && pass == entry.pass) return &entry;
    }
    return nullptr;
}

std::string no_model_message(int run, int pass) {
    std::string msg =
        "pi0::photonid::select_model: no GBT photon classifier is trained for run " + std::to_string(run) +
        " (pass " + std::to_string(pass) + ").";

    if (run >= kRgdRunMin && run <= kRgdRunMax) {
        msg += " Run " + std::to_string(run) + " is RG-D (runs " + std::to_string(kRgdRunMin) + "-" +
               std::to_string(kRgdRunMax) +
               "), and NO model has ever been trained on RG-D data: the five available models cover RG-A "
               "(2018-2019) and RG-C (Summer 2022) only.";
    } else {
        msg +=
            " The five available models cover RG-A runs 5032-5666, 6156-6783, 11093-11571 and RG-C "
            "Summer 2022 runs 16042-16772 only.";
    }

    msg +=
        " The superseded analysis SILENTLY fell back to the RGA-inbending-pass1 model here, which is"
        " WRONG for this data -- the model encodes the torus polarity, target and beam conditions it"
        " was trained under, so its photon-vs-background boundary does not transfer, and every"
        " downstream yield inherits an unquantified and unrecorded bias. Refusing rather than"
        " guessing."
        " To take the fallback anyway, set photon.allow_rga_fallback = true in config/cuts.json;"
        " fallback_used(run, pass) will then report true so the choice is stamped into the output"
        " provenance instead of vanishing.";

    return msg;
}

}  // namespace

bool fallback_used(int run, int pass) {
    return lookup(run, pass) == nullptr;
}

ModelFn select_model(int run, int pass, bool allow_rga_fallback) {
    if (const ModelRange* entry = lookup(run, pass); entry != nullptr) {
        return entry->fn();
    }

    if (!allow_rga_fallback) {
        throw std::runtime_error(no_model_message(run, pass));
    }
    return models::rga_inbending_pass1();
}

std::string_view model_name(ModelFn fn) {
    if (fn == nullptr) return "null";
    for (const auto& entry : kModelTable) {
        if (entry.fn() == fn) return entry.name;
    }
    return "unknown";
}

}  // namespace pi0::photonid
