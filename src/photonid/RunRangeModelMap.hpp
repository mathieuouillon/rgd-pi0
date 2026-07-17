#pragma once

/// \file RunRangeModelMap.hpp
/// \brief Which trained GBT photon classifier applies to which run and pass.
///
/// Ports `s_model_map` and `get_model_function` from the old
/// clas-analysis-1 PhotonCutsService, with ONE deliberate behaviour change --
/// see select_model().
///
/// The five models in photonid/models/ cover RG-A (2018-2019) and RG-C
/// (Summer 2022). They do NOT cover RG-D, which is the data this project
/// analyses. That gap is the entire reason this file is not a straight copy.

#include <string_view>
#include <vector>

namespace pi0::photonid {

/// A CatBoost applicator: 45 floats in, a raw (pre-sigmoid) margin out.
/// See PhotonGBT.hpp for turning that into a probability.
using ModelFn = double (*)(const std::vector<float>&);

/// Direct handles on the five trained models.
///
/// These exist so that callers -- and tests -- can name a model without
/// including a 5000-line generated header. Prefer select_model(); reach for
/// these only to compare against what select_model() returned.
namespace models {
[[nodiscard]] ModelFn rga_inbending_pass1();
[[nodiscard]] ModelFn rga_inbending_pass2();
[[nodiscard]] ModelFn rga_outbending_pass1();
[[nodiscard]] ModelFn rga_outbending_pass2();
[[nodiscard]] ModelFn rgc_summer2022_pass1();
}  // namespace models

/// The RG-D run range, for diagnostics. NOT a cut: nothing is selected on it.
/// It exists so the error from select_model() can say why it fired.
inline constexpr int kRgdRunMin = 18305;
inline constexpr int kRgdRunMax = 19131;

/// \return true iff NO trained model covers (run, pass) -- i.e. iff
///         select_model() would have to fall back to RGA inbending pass1.
///
/// Companion to select_model(): call it to stamp the fallback into the run's
/// provenance record, so that a plot made with the wrong model is identifiable
/// as such after the fact.
[[nodiscard]] bool fallback_used(int run, int pass);

/// Pick the GBT model trained for this run and reconstruction pass.
///
/// *** DELIBERATE BEHAVIOUR CHANGE vs the old code. ***
/// The old `get_model_function` ended with an unconditional
///
///     return [](auto const& d) { return ApplyCatboostModel_RGA_inbending_pass1(d); };
///
/// so ANY unrecognised run -- every RG-D run among them -- was silently scored
/// by a model trained on RG-A inbending data. No log line, no flag in the
/// output, nothing downstream could tell. Upstream Iguana at least emits a
/// warning. Neither is good enough for a published multiplicity ratio, so:
///
///   * allow_rga_fallback == false (config/cuts.json photon.allow_rga_fallback,
///     which is `false`): an unmatched run THROWS std::runtime_error naming the
///     run and explaining the consequence.
///   * allow_rga_fallback == true: returns the RGA inbending pass1 model, the
///     old behaviour, but fallback_used(run, pass) then returns true so the
///     caller can record it.
///
/// Matched runs are unaffected: same 16 (run range, pass) entries, same models,
/// same result as before.
///
/// \throws std::runtime_error if no model matches and allow_rga_fallback is false.
[[nodiscard]] ModelFn select_model(int run, int pass, bool allow_rga_fallback);

/// \return the model's name, e.g. "RGA_inbending_pass1", for provenance and
///         error messages. "unknown" for a pointer this map did not produce.
[[nodiscard]] std::string_view model_name(ModelFn fn);

}  // namespace pi0::photonid
