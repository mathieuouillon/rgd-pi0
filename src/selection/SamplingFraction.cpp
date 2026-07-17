#include "selection/SamplingFraction.hpp"

#include <array>
#include <stdexcept>
#include <string>

namespace pi0::selection {
namespace {

constexpr int kNumSectors = 6;

// ---------------------------------------------------------------------------
// CALIBRATION DATA -- TRANSCRIBED VERBATIM.
//
// Source: clas-analysis-1/clas12/include/clas12/services/SamplingFractionService.hpp,
//         SamplingFractionService::initialize_data(), the m_data["OB"] and
//         m_data["IB"] blocks. Six sectors each, eight coefficients each
//         (mu: a b c d, then sigma: a b c d), in sector order 1..6.
//
// These are MEASUREMENTS, not code. Every digit is significant and none has
// been rounded, reformatted or "cleaned up": the E-notation and the digit
// count below match the source character-for-character so that a diff against
// it is meaningful. Do not reflow this table. Do not fold a common factor out
// of it. If it must change, it must change because the calorimeter was
// recalibrated, and then the whole table is replaced from the new fit output.
//
// Layout note: the source keyed these on the strings "OB"/"IB" in an
// unordered_map. Here they are two flat arrays indexed by (sector - 1), which
// is the same data with the string lookup removed.
// ---------------------------------------------------------------------------

// Target-averaged Outbending parameters (6 sectors)
constexpr std::array<SectorSf, kNumSectors> kOutbending = {{
    {{2.300896E-01, 4.862074E-03, -6.441645E-04, 1.247456E-05},
     {2.226749E-02, -3.186742E-03, 3.087387E-04, -1.247691E-05}},
    {{2.286452E-01, 4.985634E-03, -3.668113E-04, -1.110160E-05},
     {2.150934E-02, -2.549354E-03, 1.855833E-04, -4.394923E-06}},
    {{2.258053E-01, 9.073840E-03, -1.240745E-03, 3.579321E-05},
     {2.399538E-02, -4.208668E-03, 4.503779E-04, -1.783519E-05}},
    {{2.281016E-01, 6.073295E-03, -5.916479E-04, 2.569914E-06},
     {2.217527E-02, -2.788589E-03, 2.266746E-04, -7.717066E-06}},
    {{2.302059E-01, 3.020237E-03, -4.307581E-05, -2.216133E-05},
     {2.277357E-02, -3.941980E-03, 4.455085E-04, -1.912118E-05}},
    {{2.262583E-01, 6.849767E-03, -8.027577E-04, 1.288497E-05},
     {2.253050E-02, -3.334553E-03, 3.374127E-04, -1.413483E-05}},
}};

// Target-averaged Inbending parameters (6 sectors)
//
// NOT USED BY THIS PROJECT TODAY -- cuts.json is outbending-only, and
// beam.polarity is the single switch. Carried because it is part of the same
// transcribed table; selecting it requires changing the config, not this file.
constexpr std::array<SectorSf, kNumSectors> kInbending = {{
    {{2.311088E-01, 4.309933E-03, -5.451627E-04, 4.965823E-06},
     {2.235239E-02, -3.117249E-03, 2.887205E-04, -1.130317E-05}},
    {{2.304337E-01, 3.585716E-03, -1.672223E-04, -7.317175E-06},
     {2.170118E-02, -2.428037E-03, 1.302147E-04, -2.503248E-09}},
    {{2.269123E-01, 8.797582E-03, -1.171416E-03, 3.051342E-05},
     {2.423467E-02, -4.301928E-03, 4.742880E-04, -1.983022E-05}},
    {{2.288870E-01, 6.786520E-03, -8.022602E-04, 1.557251E-05},
     {2.313626E-02, -3.260271E-03, 3.018634E-04, -1.156770E-05}},
    {{2.288645E-01, 4.637362E-03, -3.662892E-04, -4.233509E-06},
     {2.310895E-02, -3.892876E-03, 4.115173E-04, -1.656769E-05}},
    {{2.272093E-01, 6.581981E-03, -8.197177E-04, 1.632178E-05},
     {2.321748E-02, -3.665685E-03, 3.804588E-04, -1.587616E-05}},
}};

[[nodiscard]] bool is_valid_sector(int sector) { return sector >= 1 && sector <= kNumSectors; }

}  // namespace

Polarity polarity_from_string(const std::string& polarity) {
    if (polarity == "outbending") return Polarity::Outbending;
    if (polarity == "inbending") return Polarity::Inbending;
    throw std::runtime_error("pi0::selection::polarity_from_string: unknown polarity '" + polarity +
                             "'; expected 'inbending' or 'outbending' (cuts.json beam.polarity)");
}

const char* to_string(Polarity polarity) {
    return polarity == Polarity::Outbending ? "outbending" : "inbending";
}

const SectorSf& sf_params(int sector, Polarity polarity) {
    if (!is_valid_sector(sector)) {
        throw std::out_of_range("pi0::selection::sf_params: sector " + std::to_string(sector) +
                                " is outside [1, 6]");
    }
    const auto& table = (polarity == Polarity::Outbending) ? kOutbending : kInbending;
    return table[static_cast<std::size_t>(sector - 1)];
}

double sf_mu(double p_gev, int sector, Polarity polarity) {
    return sf_params(sector, polarity).mu.eval(p_gev);
}

double sf_sigma(double p_gev, int sector, Polarity polarity) {
    return sf_params(sector, polarity).sigma.eval(p_gev);
}

bool pass(double sf, double p_gev, int sector, Polarity polarity, double n_sigma) {
    // A track with no valid PCAL sector has no band to be tested against. The
    // old service returned false here too; keep that, and keep it BEFORE
    // sf_params() so this function never throws.
    if (!is_valid_sector(sector)) return false;

    const SectorSf& params = sf_params(sector, polarity);
    const double mu = params.mu.eval(p_gev);
    const double half_width = n_sigma * params.sigma.eval(p_gev);

    // Strict at both ends, matching every other comparison in this analysis.
    return (mu - half_width) < sf && sf < (mu + half_width);
}

}  // namespace pi0::selection
