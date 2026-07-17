#pragma once

/// \file SamplingFraction.hpp
/// \brief Per-sector, per-polarity sampling-fraction bands for electron ID.
///
/// The sampling fraction is
///
///     SF = (E_PCAL + E_ECIN + E_ECOUT) / p
///
/// and an electron is required to sit inside a band around a momentum-dependent
/// mean:
///
///     mu(p) - n_sigma * sigma(p)  <  SF  <  mu(p) + n_sigma * sigma(p)
///
/// with mu and sigma both third-order polynomials in p. This is the ONLY e/pi
/// separation in the analysis.
///
/// WHY THE COEFFICIENTS ARE NOT IN cuts.json: they are a detector calibration,
/// not a cut. cuts.json owns the one knob that selects (`n_sigma`); the
/// polynomials describe how the calorimeter responds and would belong in a
/// calibration file, which this project does not have yet. See the
/// `_coefficients_comment` in cuts.json, which states this explicitly. If a
/// calibration file is ever added, move the tables in the .cpp into it -- but
/// do not move them into cuts.json.
///
/// The coefficients are TARGET-AVERAGED: one set per (sector, polarity),
/// applied to LD2, C, Cu and Sn alike.
///
/// Units: momenta in GeV. SF is dimensionless. No angles appear here.

#include <string>

namespace pi0::selection {

/// Torus polarity. Spelled "inbending"/"outbending" in cuts.json's
/// `beam.polarity`; this project ships outbending, but both tables are present
/// because they were both in the source and dropping one would be a silent
/// edit of a data table.
enum class Polarity { Inbending, Outbending };

/// Parse cuts.json's `beam.polarity` spelling.
/// \throws std::runtime_error on anything but "inbending"/"outbending".
[[nodiscard]] Polarity polarity_from_string(const std::string& polarity);

/// The cuts.json spelling of `polarity`, i.e. the inverse of the above.
[[nodiscard]] const char* to_string(Polarity polarity);

/// A third-order polynomial, a + b*x + c*x^2 + d*x^3.
struct Poly3 {
    double a{}, b{}, c{}, d{};

    [[nodiscard]] double eval(double x) const { return a + b * x + c * x * x + d * x * x * x; }
};

/// The mu(p) and sigma(p) curves for one sector at one polarity.
struct SectorSf {
    Poly3 mu{};
    Poly3 sigma{};
};

/// Calibration coefficients for `sector`.
/// \param sector CLAS12 sector, 1-6 (PCAL sector, i.e. REC::Calorimeter.sector
///               for the PCAL layer -- not the DC sector).
/// \throws std::out_of_range if `sector` is outside [1, 6]. Callers holding an
///         unvalidated sector should use pass(), which returns false instead.
[[nodiscard]] const SectorSf& sf_params(int sector, Polarity polarity);

/// Mean sampling fraction at momentum `p_gev`. \throws as sf_params().
[[nodiscard]] double sf_mu(double p_gev, int sector, Polarity polarity);

/// Sampling-fraction width at momentum `p_gev`. \throws as sf_params().
[[nodiscard]] double sf_sigma(double p_gev, int sector, Polarity polarity);

/// The sampling-fraction cut itself:
///
///     mu(p) - n_sigma*sigma(p) < sf < mu(p) + n_sigma*sigma(p)
///
/// STRICT at both ends, and symmetric about mu by construction.
///
/// \param sf       the measured (E_PCAL + E_ECIN + E_ECOUT) / p.
/// \param p_gev    the track momentum used in that ratio, GeV.
/// \param sector   PCAL sector. Returns FALSE (does not throw) outside [1, 6]:
///                 a track with no valid PCAL sector has no band to be tested
///                 against, and "no band" means "not an identified electron".
/// \param n_sigma  the band half-width, from cuts.electron.sf_n_sigma. Never
///                 hard-code it.
[[nodiscard]] bool pass(double sf, double p_gev, int sector, Polarity polarity, double n_sigma);

}  // namespace pi0::selection
