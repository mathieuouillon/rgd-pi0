/// \file test_photonid.cpp
/// \brief Unit tests for src/photonid -- the ported CLAS12 GBT photon classifier.
///
/// The port's contract is "bit-identical score to the old PhotonCutsService for
/// the same photon", so these tests pin the things a refactor could silently
/// move: the length of the feature vector, the ORDER of its 45 elements, which
/// neighbours are eligible, and the run -> model mapping. They deliberately do
/// NOT assert on any model's numeric output -- that is the models' business and
/// they are verbatim upstream code.
///
/// No data file, no HIPO: the RVec columns are hand-built here.

#include <cmath>
#include <string>
#include <vector>

#include <ROOT/RVec.hxx>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "photonid/Features.hpp"
#include "photonid/PhotonGBT.hpp"
#include "photonid/RunRangeModelMap.hpp"

using Catch::Approx;
using namespace pi0::photonid;

namespace {

template <typename T>
using RVec = ROOT::VecOps::RVec<T>;

constexpr double kDeg = M_PI / 180.0;

/// The neighbour window, matching config/cuts.json photon.{min_energy_gev,
/// theta_min_deg, theta_max_deg}. Spelled out here rather than loaded because
/// these tests pin the CODE's behaviour for a GIVEN window; the loading of the
/// window is Cuts::load's job and is tested there.
constexpr FeatureCuts kCuts{0.2, 5.0, 35.0};

/// A synthetic event, assembled row by row.
///
/// GEOMETRY WARNING: the momentum direction and the calorimeter position of a
/// row are set INDEPENDENTLY here, and are not physically consistent with each
/// other. That is intentional and it is exactly what the code under test does:
/// the theta window is applied to the MOMENTUM, while the neighbour angle R is
/// computed from the CALORIMETER POSITION. Decoupling them lets each be pinned
/// without the other moving.
struct EventBuilder {
    RVec<int> pid;
    RVec<float> px, py, pz;

    RVec<short> calo_pindex, calo_layer, calo_sector;
    RVec<float> calo_e, calo_x, calo_y, calo_z, calo_m2u, calo_m2v, calo_lu, calo_lv, calo_lw;

    /// \param momentum_gev  |p|; the direction is fixed at theta_deg, phi = 0.
    /// \return the row index just added.
    std::size_t add_particle(int particle_pid, double momentum_gev, double theta_deg) {
        pid.push_back(particle_pid);
        px.push_back(static_cast<float>(momentum_gev * std::sin(theta_deg * kDeg)));
        py.push_back(0.0F);
        pz.push_back(static_cast<float>(momentum_gev * std::cos(theta_deg * kDeg)));
        return pid.size() - 1;
    }

    void add_calo(std::size_t row, int layer, double e, double x, double z, double m2u, double m2v) {
        calo_pindex.push_back(static_cast<short>(row));
        calo_layer.push_back(static_cast<short>(layer));
        calo_sector.push_back(1);
        calo_e.push_back(static_cast<float>(e));
        calo_x.push_back(static_cast<float>(x));
        calo_y.push_back(0.0F);
        calo_z.push_back(static_cast<float>(z));
        calo_m2u.push_back(static_cast<float>(m2u));
        calo_m2v.push_back(static_cast<float>(m2v));
        calo_lu.push_back(20.0F);
        calo_lv.push_back(20.0F);
        calo_lw.push_back(20.0F);
    }

    [[nodiscard]] CaloMap map() const {
        return CaloMap::build(calo_pindex, calo_layer, calo_sector, calo_e, calo_x, calo_y, calo_z,
                              calo_m2u, calo_m2v, calo_lu, calo_lv, calo_lw);
    }
};

/// Calorimeter face distance, cm, and the candidate's cluster x, cm.
///
/// kXPoi IS NOT ALLOWED TO BE ZERO, and that is not a detail of the test: the
/// ported calo_position() decides "this particle has no PCAL cluster" by
/// testing `pcal.x == 0`, so a candidate parked on the beam axis falls through
/// PCAL -> ECIN -> ECOUT to (0,0,0), and every neighbour angle collapses to
/// acos(0) = pi/2. That is faithful old-code behaviour (see the dedicated test
/// below); here we just keep the fixture out of that corner.
constexpr double kZ = 700.0;
constexpr double kXPoi = 70.0;

/// The expected angle, RADIANS, between the candidate's cluster at
/// (kXPoi, 0, kZ) and a neighbour's at (x, 0, kZ).
///
/// Computed independently of the code under test: both vectors lie in the x-z
/// half-plane with z = kZ > 0, so each makes a polar angle atan(x / kZ) with
/// the z axis and the angle between them is the difference.
double expected_R(double x) {
    return std::abs(std::atan(x / kZ) - std::atan(kXPoi / kZ));
}

}  // namespace

// ===========================================================================
// Feature vector: length and ordering
// ===========================================================================

TEST_CASE("feature vector is exactly 45 long", "[photonid][features]") {
    EventBuilder ev;
    const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
    ev.add_calo(poi, kLayerPcal, 0.9, kXPoi, kZ, 5.5, 6.5);

    // A lone photon with no neighbours at all still yields the full vector,
    // padded with the "empty slot" zeros.
    const auto feats = build_features(poi, ev.map(), ev.pid, ev.px, ev.py, ev.pz, kCuts);

    REQUIRE(feats.size() == 45);
    REQUIRE(feats.size() == kNumFeatures);

    for (std::size_t i = 5; i < feats.size(); ++i) {
        INFO("neighbour slot " << i << " must be zero when the event has no neighbours");
        CHECK(feats[i] == 0.0F);
    }
}

TEST_CASE("feature vector ordering is pinned element by element", "[photonid][features]") {
    EventBuilder ev;

    // --- the candidate -----------------------------------------------------
    // |p| = 2 GeV at theta = 10 deg, PCAL cluster at (kXPoi, 0, kZ).
    const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
    ev.add_calo(poi, kLayerPcal, 0.9, kXPoi, kZ, 5.5, 6.5);

    // --- neighbour photons, added OUT of R order, so the insertion sort has
    //     something to do. Expected sorted order is A (nearest), B, C, with
    //     R = 0.0099, 0.152, 0.293 rad -- straddling the 0.1 / 0.2 / 0.35
    //     counter boundaries one each.
    const std::size_t gam_b = ev.add_particle(22, 0.5, 10.0);  // R ~ 0.152
    ev.add_calo(gam_b, kLayerPcal, 0.4, 180.0, kZ, 2.1, 2.2);

    const std::size_t gam_c = ev.add_particle(22, 0.3, 10.0);  // R ~ 0.293
    ev.add_calo(gam_c, kLayerPcal, 0.3, 290.0, kZ, 3.1, 3.2);

    const std::size_t gam_a = ev.add_particle(22, 1.0, 10.0);  // R ~ 0.0099
    ev.add_calo(gam_a, kLayerPcal, 0.5, 77.0, kZ, 1.1, 1.2);

    // --- neighbour photons that must be EXCLUDED ---------------------------
    // All three sit CLOSER than A, so a leak is loud: it would take slot 0 and
    // bump the num_photons counters.
    const std::size_t gam_soft = ev.add_particle(22, 0.1, 10.0);  // E < min_energy_gev
    ev.add_calo(gam_soft, kLayerPcal, 0.05, 73.0, kZ, 9.1, 9.2);

    const std::size_t gam_nopcal = ev.add_particle(22, 1.0, 10.0);  // Epcal == 0
    ev.add_calo(gam_nopcal, kLayerEcin, 0.7, 74.0, kZ, 9.3, 9.4);   // ECIN only

    const std::size_t gam_fwd = ev.add_particle(22, 1.0, 40.0);  // theta > theta_max_deg
    ev.add_calo(gam_fwd, kLayerPcal, 0.6, 75.0, kZ, 9.5, 9.6);

    // --- charged hadrons, out of order: expected [ch_near, ch_far] ---------
    const std::size_t ch_far = ev.add_particle(2212, 1.5, 10.0);  // R ~ 0.107
    ev.add_calo(ch_far, kLayerPcal, 0.25, 147.0, kZ, 4.1, 4.2);

    const std::size_t ch_near = ev.add_particle(211, 0.8, 10.0);  // R ~ 0.0198
    ev.add_calo(ch_near, kLayerPcal, 0.15, 84.0, kZ, 5.1, 5.2);

    // --- neutral hadrons, out of order: expected [nh_near, nh_far] ---------
    const std::size_t nh_far = ev.add_particle(-2112, 1.2, 10.0);  // R ~ 0.192
    ev.add_calo(nh_far, kLayerPcal, 0.35, 210.0, kZ, 6.1, 6.2);

    const std::size_t nh_near = ev.add_particle(2112, 0.9, 10.0);  // R ~ 0.0492
    ev.add_calo(nh_near, kLayerPcal, 0.45, 105.0, kZ, 7.1, 7.2);

    // --- two electrons: features [5]/[6] take the NEARER one ---------------
    const std::size_t e_far = ev.add_particle(11, 4.0, 10.0);  // R ~ 0.441
    ev.add_calo(e_far, kLayerPcal, 1.8, 420.0, kZ, 8.1, 8.2);

    const std::size_t e_near = ev.add_particle(11, 3.0, 10.0);  // R ~ 0.0296
    ev.add_calo(e_near, kLayerPcal, 1.5, 91.0, kZ, 8.3, 8.4);

    // --- rows that must be skipped entirely --------------------------------
    const std::size_t unknown = ev.add_particle(45, 1.0, 10.0);  // pid -> NaN mass
    ev.add_calo(unknown, kLayerPcal, 0.5, 71.0, kZ, 9.7, 9.8);

    ev.add_particle(2212, 1.0, 10.0);  // no calorimeter row at all -> skipped

    const auto feats = build_features(poi, ev.map(), ev.pid, ev.px, ev.py, ev.pz, kCuts);
    REQUIRE(feats.size() == kNumFeatures);

    SECTION("[0..4] candidate's own quantities") {
        CHECK(feats[0] == Approx(2.0).margin(1e-6));            // gE = |p|
        CHECK(feats[1] == Approx(0.9).margin(1e-6));            // gEpcal
        CHECK(feats[2] == Approx(10.0 * kDeg).margin(1e-6));    // gTheta, RADIANS not degrees
        CHECK(feats[3] == Approx(5.5).margin(1e-6));            // gm2u
        CHECK(feats[4] == Approx(6.5).margin(1e-6));            // gm2v
    }

    SECTION("[2] gTheta really is radians") {
        // The single most valuable assertion in this file: 10 deg is 0.1745 rad,
        // and the models were trained on the radian. Guard both directions.
        CHECK(feats[2] == Approx(0.17453293).margin(1e-6));
        CHECK(feats[2] != Approx(10.0).margin(1e-3));
    }

    SECTION("[5..6] nearest electron, of ANY pid==11 row") {
        CHECK(feats[5] == Approx(expected_R(91.0)).margin(1e-6));  // e_near, not e_far
        CHECK(feats[6] == Approx(2.0 - 3.0).margin(1e-5));         // gE - E_e_near
    }

    SECTION("[7..21] photon neighbours, grouped by variable and sorted by R") {
        // R_gamma[0..2] -- ascending, and the three excluded photons are absent.
        CHECK(feats[7] == Approx(expected_R(77.0)).margin(1e-6));
        CHECK(feats[8] == Approx(expected_R(180.0)).margin(1e-6));
        CHECK(feats[9] == Approx(expected_R(290.0)).margin(1e-6));
        CHECK(feats[7] < feats[8]);
        CHECK(feats[8] < feats[9]);

        // dE_gamma[0..2] = gE - E_neighbour
        CHECK(feats[10] == Approx(2.0 - 1.0).margin(1e-5));
        CHECK(feats[11] == Approx(2.0 - 0.5).margin(1e-5));
        CHECK(feats[12] == Approx(2.0 - 0.3).margin(1e-5));

        // Epcal_gamma[0..2]
        CHECK(feats[13] == Approx(0.5).margin(1e-6));
        CHECK(feats[14] == Approx(0.4).margin(1e-6));
        CHECK(feats[15] == Approx(0.3).margin(1e-6));

        // m2u_gamma[0..2]
        CHECK(feats[16] == Approx(1.1).margin(1e-6));
        CHECK(feats[17] == Approx(2.1).margin(1e-6));
        CHECK(feats[18] == Approx(3.1).margin(1e-6));

        // m2v_gamma[0..2]
        CHECK(feats[19] == Approx(1.2).margin(1e-6));
        CHECK(feats[20] == Approx(2.2).margin(1e-6));
        CHECK(feats[21] == Approx(3.2).margin(1e-6));
    }

    SECTION("[7..21] are grouped BY VARIABLE, not by neighbour") {
        // If the layout were ever "flattened" to {R,dE,Epcal,m2u,m2v} per
        // neighbour, feats[8] would be dE of the nearest photon (1.0) rather
        // than R of the second-nearest (0.149). This is the assertion that
        // catches that regression.
        CHECK(feats[8] == Approx(expected_R(180.0)).margin(1e-6));
        CHECK(feats[8] != Approx(1.0).margin(1e-3));
    }

    SECTION("[22..31] charged-hadron neighbours (pi+, p), sorted by R") {
        CHECK(feats[22] == Approx(expected_R(84.0)).margin(1e-6));   // ch_near
        CHECK(feats[23] == Approx(expected_R(147.0)).margin(1e-6));   // ch_far
        // dE uses E = sqrt(p^2 + m^2): pi+ at |p| = 0.8 with m = 0.1396,
        // proton at |p| = 1.5 with m = 0.938272.
        CHECK(feats[24] == Approx(2.0 - std::hypot(0.8, 0.1396)).margin(1e-5));
        CHECK(feats[25] == Approx(2.0 - std::hypot(1.5, 0.938272)).margin(1e-5));
        CHECK(feats[26] == Approx(0.15).margin(1e-6));
        CHECK(feats[27] == Approx(0.25).margin(1e-6));
        CHECK(feats[28] == Approx(5.1).margin(1e-6));
        CHECK(feats[29] == Approx(4.1).margin(1e-6));
        CHECK(feats[30] == Approx(5.2).margin(1e-6));
        CHECK(feats[31] == Approx(4.2).margin(1e-6));
    }

    SECTION("[32..41] neutral-hadron neighbours (n, nbar), sorted by R") {
        CHECK(feats[32] == Approx(expected_R(105.0)).margin(1e-6));   // nh_near
        CHECK(feats[33] == Approx(expected_R(210.0)).margin(1e-6));  // nh_far
        CHECK(feats[34] == Approx(2.0 - std::hypot(0.9, 0.93956536)).margin(1e-5));
        CHECK(feats[35] == Approx(2.0 - std::hypot(1.2, 0.93956536)).margin(1e-5));
        CHECK(feats[36] == Approx(0.45).margin(1e-6));
        CHECK(feats[37] == Approx(0.35).margin(1e-6));
        CHECK(feats[38] == Approx(7.1).margin(1e-6));
        CHECK(feats[39] == Approx(6.1).margin(1e-6));
        CHECK(feats[40] == Approx(7.2).margin(1e-6));
        CHECK(feats[41] == Approx(6.2).margin(1e-6));
    }

    SECTION("[42..44] photon counters count only eligible photons") {
        // R = 0.0100, 0.1489, 0.2915 rad for A, B, C.
        CHECK(feats[42] == Approx(1.0).margin(1e-6));  // R < 0.10 : A
        CHECK(feats[43] == Approx(2.0).margin(1e-6));  // R < 0.20 : A, B
        CHECK(feats[44] == Approx(3.0).margin(1e-6));  // R < 0.35 : A, B, C
    }
}

TEST_CASE("charged and neutral hadrons are NOT subject to the photon purity filter",
          "[photonid][features]") {
    // The old code applies the full purity filter (E, Epcal, theta) to photon
    // neighbours but only the theta window to hadrons. A hadron with no PCAL
    // energy therefore still occupies a neighbour slot. Asymmetric, and pinned
    // here so that "tidying" the two branches into one is caught.
    EventBuilder ev;
    const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
    ev.add_calo(poi, kLayerPcal, 0.9, kXPoi, kZ, 5.5, 6.5);

    const std::size_t ch = ev.add_particle(211, 0.8, 10.0);
    ev.add_calo(ch, kLayerEcin, 0.5, 84.0, kZ, 1.0, 2.0);  // ECIN only -> Epcal = 0

    const auto feats = build_features(poi, ev.map(), ev.pid, ev.px, ev.py, ev.pz, kCuts);

    CHECK(feats[22] == Approx(expected_R(84.0)).margin(1e-6));  // slot filled
    CHECK(feats[26] == 0.0F);                                   // ...with Epcal = 0
}

TEST_CASE("calorimeter position falls back PCAL -> ECIN -> ECOUT", "[photonid][features]") {
    EventBuilder ev;
    const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
    ev.add_calo(poi, kLayerPcal, 0.9, kXPoi, kZ, 5.5, 6.5);

    // Neighbour with ECIN and ECOUT but no PCAL: R must come from ECIN.
    const std::size_t nh = ev.add_particle(2112, 0.9, 10.0);
    ev.add_calo(nh, kLayerEcin, 0.5, 105.0, kZ, 1.0, 2.0);
    ev.add_calo(nh, kLayerEcout, 0.5, 999.0, kZ, 3.0, 4.0);

    const auto feats = build_features(poi, ev.map(), ev.pid, ev.px, ev.py, ev.pz, kCuts);
    CHECK(feats[32] == Approx(expected_R(105.0)).margin(1e-6));
}

TEST_CASE("a cluster at x == 0 is read as 'no cluster in this layer'", "[photonid][features]") {
    // A WART, PORTED ON PURPOSE. get_particle_calo_vector() in the old code
    // picks the layer by testing `pcal_x == 0`, not by asking whether the layer
    // exists. So a candidate whose PCAL cluster genuinely sits at x = 0 falls
    // through PCAL -> ECIN -> ECOUT; with no ECIN/ECOUT either, its position
    // vector is (0,0,0), CosTheta takes its ptot2 <= 0 branch, and EVERY
    // neighbour lands at acos(0) = pi/2 exactly.
    //
    // Physically x is a float from reconstruction and never lands on exactly
    // 0.0, so this costs nothing on real data -- but it is what the models were
    // trained against, so it stays, and it is pinned here rather than left as a
    // trap for the next person to rediscover (which is how this test was born).
    EventBuilder ev;
    const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
    ev.add_calo(poi, kLayerPcal, 0.9, 0.0, kZ, 5.5, 6.5);  // x == 0, no other layer

    const std::size_t gam = ev.add_particle(22, 1.0, 10.0);
    ev.add_calo(gam, kLayerPcal, 0.5, 77.0, kZ, 1.1, 1.2);

    const auto feats = build_features(poi, ev.map(), ev.pid, ev.px, ev.py, ev.pz, kCuts);

    CHECK(feats[7] == Approx(M_PI / 2.0).margin(1e-6));
    // ...and pi/2 > 0.35, so the counters see nothing despite the neighbour
    // being 1 cm away on the calorimeter face.
    CHECK(feats[42] == 0.0F);
    CHECK(feats[43] == 0.0F);
    CHECK(feats[44] == 0.0F);

    // The candidate's OWN features are unaffected: they come from the momentum
    // and the PCAL energy, neither of which routes through calo_position().
    CHECK(feats[0] == Approx(2.0).margin(1e-6));
    CHECK(feats[1] == Approx(0.9).margin(1e-6));
}

TEST_CASE("the PID purity pre-filter gates the candidate", "[photonid][features]") {
    SECTION("a good candidate passes") {
        EventBuilder ev;
        const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
        ev.add_calo(poi, kLayerPcal, 0.9, kXPoi, kZ, 5.5, 6.5);
        CHECK(passes_pid_purity(poi, ev.map(), ev.px, ev.py, ev.pz, kCuts));
    }

    SECTION("no calorimeter data at all -> rejected, not thrown") {
        EventBuilder ev;
        const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
        CHECK_FALSE(passes_pid_purity(poi, ev.map(), ev.px, ev.py, ev.pz, kCuts));
        // ...and build_features on it is a caller bug, so it throws.
        CHECK_THROWS_AS(build_features(poi, ev.map(), ev.pid, ev.px, ev.py, ev.pz, kCuts),
                        std::invalid_argument);
    }

    SECTION("below min_energy_gev -> rejected") {
        EventBuilder ev;
        const std::size_t poi = ev.add_particle(22, 0.1, 10.0);
        ev.add_calo(poi, kLayerPcal, 0.05, 0.0, kZ, 5.5, 6.5);
        CHECK_FALSE(passes_pid_purity(poi, ev.map(), ev.px, ev.py, ev.pz, kCuts));
    }

    SECTION("no PCAL energy -> rejected") {
        EventBuilder ev;
        const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
        ev.add_calo(poi, kLayerEcin, 0.9, 0.0, kZ, 5.5, 6.5);
        CHECK_FALSE(passes_pid_purity(poi, ev.map(), ev.px, ev.py, ev.pz, kCuts));
    }

    SECTION("outside the theta window -> rejected on both sides") {
        for (const double theta_deg : {2.0, 40.0}) {
            EventBuilder ev;
            const std::size_t poi = ev.add_particle(22, 2.0, theta_deg);
            ev.add_calo(poi, kLayerPcal, 0.9, kXPoi, kZ, 5.5, 6.5);
            INFO("theta = " << theta_deg << " deg");
            CHECK_FALSE(passes_pid_purity(poi, ev.map(), ev.px, ev.py, ev.pz, kCuts));
        }
    }
}

TEST_CASE("CaloMap distinguishes absent from present-but-no-PCAL", "[photonid][features]") {
    EventBuilder ev;
    const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
    ev.add_calo(poi, kLayerEcin, 0.9, 0.0, kZ, 5.5, 6.5);
    ev.add_particle(22, 2.0, 10.0);  // row 1, no calorimeter row

    const CaloMap map = ev.map();
    CHECK(map.size() == 1);
    CHECK(map.contains(0));
    CHECK_FALSE(map.contains(1));
    REQUIRE(map.find(0) != nullptr);
    CHECK(map.find(0)->pcal.e == 0.0);   // present, but no PCAL
    CHECK(map.find(0)->ecin.e == Approx(0.9));
    CHECK(map.find(1) == nullptr);
}

TEST_CASE("ragged input columns are rejected loudly", "[photonid][features]") {
    EventBuilder ev;
    const std::size_t poi = ev.add_particle(22, 2.0, 10.0);
    ev.add_calo(poi, kLayerPcal, 0.9, kXPoi, kZ, 5.5, 6.5);

    RVec<float> short_pz{};  // one row too few
    CHECK_THROWS_AS(build_features(poi, ev.map(), ev.pid, ev.px, ev.py, short_pz, kCuts),
                    std::invalid_argument);
    CHECK_THROWS_AS(build_features(99, ev.map(), ev.pid, ev.px, ev.py, ev.pz, kCuts),
                    std::invalid_argument);
}

// ===========================================================================
// Run -> model map
// ===========================================================================

TEST_CASE("select_model returns the model trained for each run range", "[photonid][model]") {
    struct Expectation {
        int run;
        int pass;
        ModelFn expected;
        const char* name;
    };

    // One representative run from the interior of each of the 16 entries in the
    // old s_model_map.
    const Expectation cases[] = {
        {5100, 1, models::rga_inbending_pass1(), "RGA inbending, 5032-5332 pass1"},
        {5100, 2, models::rga_inbending_pass2(), "RGA inbending, 5032-5332 pass2"},
        {5500, 1, models::rga_outbending_pass1(), "RGA outbending, 5333-5666 pass1"},
        {5500, 2, models::rga_outbending_pass2(), "RGA outbending, 5333-5666 pass2"},
        {6700, 1, models::rga_inbending_pass1(), "RGA inbending, 6616-6783 pass1"},
        {6700, 2, models::rga_inbending_pass2(), "RGA inbending, 6616-6783 pass2"},
        {6300, 1, models::rga_inbending_pass1(), "RGA inbending, 6156-6603 pass1"},
        {6300, 2, models::rga_inbending_pass2(), "RGA inbending, 6156-6603 pass2"},
        {11200, 1, models::rga_outbending_pass1(), "RGA outbending, 11093-11283 pass1"},
        {11200, 2, models::rga_outbending_pass2(), "RGA outbending, 11093-11283 pass2"},
        {11290, 1, models::rga_inbending_pass1(), "RGA inbending, 11284-11300 pass1"},
        {11290, 2, models::rga_inbending_pass2(), "RGA inbending, 11284-11300 pass2"},
        {11400, 1, models::rga_inbending_pass1(), "RGA inbending, 11323-11571 pass1"},
        {11400, 2, models::rga_inbending_pass2(), "RGA inbending, 11323-11571 pass2"},
        {16500, 1, models::rgc_summer2022_pass1(), "RGC Summer2022, 16042-16772 pass1"},
        // Upstream maps RG-C pass2 to the pass1 model -- no pass2 model exists.
        {16500, 2, models::rgc_summer2022_pass1(), "RGC Summer2022, 16042-16772 pass2 -> pass1 model"},
    };

    for (const auto& c : cases) {
        INFO(c.name << " (run " << c.run << ", pass " << c.pass << ")");
        CHECK(select_model(c.run, c.pass, false) == c.expected);
        CHECK_FALSE(fallback_used(c.run, c.pass));
    }
}

TEST_CASE("select_model pins the exact range boundaries", "[photonid][model]") {
    // Inclusive on both ends, and the neighbouring run is a different model --
    // an off-by-one in a boundary silently rescores a whole run.
    CHECK(select_model(5032, 1, false) == models::rga_inbending_pass1());
    CHECK(select_model(5332, 1, false) == models::rga_inbending_pass1());
    CHECK(select_model(5333, 1, false) == models::rga_outbending_pass1());
    CHECK(select_model(5666, 1, false) == models::rga_outbending_pass1());
    CHECK(select_model(16042, 1, false) == models::rgc_summer2022_pass1());
    CHECK(select_model(16772, 1, false) == models::rgc_summer2022_pass1());

    // The gaps between ranges are real: no model, so it must refuse.
    CHECK(fallback_used(5031, 1));   // below the first range
    CHECK(fallback_used(5667, 1));   // 5667-6155 gap
    CHECK(fallback_used(6604, 1));   // 6604-6615 gap
    CHECK(fallback_used(11301, 1));  // 11301-11322 gap
    CHECK(fallback_used(16773, 1));  // above the last range

    // Pass 3 exists for no range at all.
    CHECK(fallback_used(5100, 3));
    CHECK_THROWS_AS(select_model(5100, 3, false), std::runtime_error);
}

TEST_CASE("select_model REFUSES an RG-D run when the fallback is not allowed", "[photonid][model]") {
    // THE headline behaviour change vs the old code, which returned the
    // RGA-inbending-pass1 model here without a word. RG-D is the data this
    // project analyses, so this is not a corner case -- it is the default path.
    CHECK_THROWS_AS(select_model(18500, 1, false), std::runtime_error);

    CHECK_THROWS_WITH(
        select_model(18500, 1, false),
        Catch::Matchers::ContainsSubstring("18500") &&
            Catch::Matchers::ContainsSubstring("RG-D") &&
            Catch::Matchers::ContainsSubstring("allow_rga_fallback"));

    // Every RG-D run, not just the representative one.
    CHECK_THROWS_AS(select_model(kRgdRunMin, 1, false), std::runtime_error);
    CHECK_THROWS_AS(select_model(kRgdRunMax, 1, false), std::runtime_error);
    CHECK_THROWS_AS(select_model(18500, 2, false), std::runtime_error);
}

TEST_CASE("select_model takes the RGA fallback when allowed, and reports it", "[photonid][model]") {
    ModelFn fn{};
    CHECK_NOTHROW(fn = select_model(18500, 1, true));

    CHECK(fn == models::rga_inbending_pass1());
    CHECK(model_name(fn) == "RGA_inbending_pass1");

    // The point of the change: the caller can SEE that it happened and stamp it
    // into the output's provenance.
    CHECK(fallback_used(18500, 1));
    CHECK(fallback_used(18500, 2));

    // A matched run reports no fallback, so the flag means what it says.
    CHECK_FALSE(fallback_used(5100, 1));
    CHECK(select_model(5100, 1, true) == models::rga_inbending_pass1());  // same fn, but matched
}

TEST_CASE("model_name reports each model for provenance", "[photonid][model]") {
    CHECK(model_name(models::rga_inbending_pass1()) == "RGA_inbending_pass1");
    CHECK(model_name(models::rga_inbending_pass2()) == "RGA_inbending_pass2");
    CHECK(model_name(models::rga_outbending_pass1()) == "RGA_outbending_pass1");
    CHECK(model_name(models::rga_outbending_pass2()) == "RGA_outbending_pass2");
    CHECK(model_name(models::rgc_summer2022_pass1()) == "RGC_Summer2022_pass1");
    CHECK(model_name(nullptr) == "null");
}

// ===========================================================================
// Scoring
// ===========================================================================

TEST_CASE("score applies the logistic to the model's margin", "[photonid][gbt]") {
    // A stub model standing in for a CatBoost applicator, so the sigmoid is
    // tested independently of any trained weights.
    struct Stub {
        static double zero(const std::vector<float>&) { return 0.0; }
        static double big_positive(const std::vector<float>&) { return 40.0; }
        static double big_negative(const std::vector<float>&) { return -40.0; }
        static double one(const std::vector<float>&) { return 1.0; }
        static double minus_one(const std::vector<float>&) { return -1.0; }
    };

    const std::vector<float> feats(kNumFeatures, 0.0F);

    CHECK(score(feats, &Stub::zero) == Approx(0.5));
    CHECK(score(feats, &Stub::one) == Approx(1.0 / (1.0 + std::exp(-1.0))));
    CHECK(score(feats, &Stub::one) == Approx(0.7310585786));
    CHECK(score(feats, &Stub::minus_one) == Approx(1.0 - 0.7310585786));

    // Saturating, but never outside [0, 1] and never NaN.
    CHECK(score(feats, &Stub::big_positive) == Approx(1.0).margin(1e-12));
    CHECK(score(feats, &Stub::big_negative) == Approx(0.0).margin(1e-12));
    CHECK(score(feats, &Stub::big_positive) <= 1.0);
    CHECK(score(feats, &Stub::big_negative) >= 0.0);

    // sigmoid(-x) == 1 - sigmoid(x)
    CHECK(score(feats, &Stub::one) + score(feats, &Stub::minus_one) == Approx(1.0));

    SECTION("a wrong-length feature vector is rejected, not overread") {
        CHECK_THROWS_AS(score(std::vector<float>(44, 0.0F), &Stub::zero), std::invalid_argument);
        CHECK_THROWS_AS(score(std::vector<float>(46, 0.0F), &Stub::zero), std::invalid_argument);
        CHECK_THROWS_AS(score({}, &Stub::zero), std::invalid_argument);
        CHECK_THROWS_AS(score(feats, nullptr), std::invalid_argument);
    }
}

TEST_CASE("passes compares strictly above the threshold", "[photonid][gbt]") {
    // The threshold is a parameter, never a constant in this library: it comes
    // from Cuts::photon.gbt_threshold. 0.78 appears here only as test data.
    constexpr double threshold = 0.78;

    CHECK(passes(0.9, threshold));
    CHECK_FALSE(passes(0.5, threshold));

    // Strictly greater, matching the old `prediction > m_ai_threshold`.
    CHECK_FALSE(passes(threshold, threshold));
    CHECK(passes(std::nextafter(threshold, 1.0), threshold));
    CHECK_FALSE(passes(std::nextafter(threshold, 0.0), threshold));

    // The threshold really is honoured, not ignored in favour of a baked-in one.
    CHECK(passes(0.6, 0.5));
    CHECK_FALSE(passes(0.6, 0.7));
}
