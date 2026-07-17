#include "selection/PhotonSelection.hpp"

namespace pi0::selection {

bool passes_gbt_prefilter(double e_gamma_gev,
                          double pcal_energy_gev,
                          double theta_deg,
                          const Cuts& cuts) {
    const auto& g = cuts.photon;

    // NOTE THE INCLUSIVE BOUNDS. Almost every other cut in this analysis is
    // strict; these three are not, and they are not by accident -- cuts.json
    // calls the energy floor out explicitly ("this bound is INCLUSIVE (>=),
    // unlike most cuts in this analysis") and the theta window "INCLUSIVE at
    // both ends". Transcribed from the old pass_pid_purity_filter, which is
    // what ran.
    //
    // E_PCAL > 0 is strict, and is not a tuned threshold: it is the structural
    // requirement that the shower started in the preshower at all. cuts.json
    // spells it `require_pcal_energy: true`, which Cuts::load() verifies.
    return e_gamma_gev >= g.min_energy_gev       // INCLUSIVE
           && pcal_energy_gev > 0.0              // strict, structural
           && theta_deg >= g.theta_min_deg       // INCLUSIVE
           && theta_deg <= g.theta_max_deg;      // INCLUSIVE
}

bool pass_photon(int pid,
                 double px,
                 double py,
                 double pz,
                 double pcal_energy_gev,
                 double pcal_lv_cm,
                 double pcal_lw_cm,
                 double beta,
                 double gbt_score,
                 const Cuts& cuts) {
    // One implementation of the cut order, in the template; this overload just
    // hands it a score that is already known. Duplicating the order here is how
    // the two forms would drift apart.
    return pass_photon_scored(pid, px, py, pz, pcal_energy_gev, pcal_lv_cm, pcal_lw_cm, beta,
                              [gbt_score] { return gbt_score; }, cuts);
}

}  // namespace pi0::selection
