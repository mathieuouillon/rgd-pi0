// Unit tests for the core physics of the RG-D pi0 analysis.
//
// These tests exist to make three classes of error impossible, because the
// analysis this code replaces shipped all three:
//
//   1. degree/radian confusion at a public boundary (twice, independently);
//   2. a Trento phi_h frame that was asserted rather than verified;
//   3. a bin-abscissa bug (out of scope for this file, see tests for Binning).
//
// Every test here is written so that it fails for a specific, named reason.
// In particular the reference values are computed independently *in the test*
// -- either as literals worked out offline, or from first principles with
// local arithmetic -- never by calling the function under test.
//
// ---------------------------------------------------------------------------
// API RECONCILIATION, 2026-07-16
// ---------------------------------------------------------------------------
// This file was authored in parallel with src/core, against an API that was
// guessed rather than read, and it drifted. The headers are the authority and
// the code under test was NOT changed; the call sites here were. What moved:
//
//   compute_dis(E', theta_deg)                  -> compute_dis(ex, ey, ez, e_energy, beam_energy)
//       The electron is now passed as a momentum VECTOR and an energy, and
//       theta_e is taken from the vector. E' and |p_e| are separate arguments
//       on purpose (a calorimeter-corrected E' must not be re-derived from the
//       tracked |p|), so the local helper `dis_at` ties them together for the
//       massless case and one test below pins them apart deliberately.
//
//   compute_sidis(ex,ey,ez, hx,hy,hz, h_e)      -> compute_sidis(hx, hy, hz, he, dis,
//                                                                ex, ey, ez, e_energy, beam_energy)
//       Hadron first, and the event DisKin is threaded in for nu rather than
//       recomputed. See the local helper `sidis_of`.
//
//   pi0::pairing::min_opening_angle_deg(p)      -> pi0::min_opening_angle_deg(p, cuts)
//   pi0::pairing::find_gg_pairs(...)            -> pi0::find_gg_pairs(...)
//       There is no `pi0::pairing` namespace; both live directly in `pi0`. The
//       opening-angle curve's (a, b, offset) are CUT VALUES carried in
//       PairingCuts, not constants baked into the function.
//
//   pi0::pairing::PairCuts                      -> pi0::PairingCuts
//       Renamed, and it has NO default member values -- it is documented as
//       "Populated from cuts.json". A default-constructed PairingCuts{} is all
//       zeros, which disables every cut, so every test below states the shipped
//       numbers explicitly via `shipped_cuts()`. Those numbers are transcribed
//       from config/cuts.json ("pairing" block), which is the single source of
//       truth; if they ever disagree with that file, that file wins.
//
//   PairCuts::min_e_gamma_deg                   -> PairingCuts::e_gamma_min_angle_deg
//       AND, more than a rename: find_gg_pairs does not apply this cut at all.
//       See the e-gamma test at the bottom of this file for the full story.
//
// The physics assertions are unchanged by any of the above -- they were
// independent of the spelling, which is what made this reconciliation
// mechanical rather than a rewrite.
//
// Physics reference: analysis note sections 03 (DIS) and 04 (pi0 reco / SIDIS).
// ---------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "core/Constants.hpp"
#include "core/Kinematics.hpp"
#include "core/Pairing.hpp"
#include "core/Types.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// ---------------------------------------------------------------------------
// Test-local vector algebra. Deliberately duplicated rather than imported:
// these tests must not share an implementation with the code under test.
// ---------------------------------------------------------------------------
struct V3 {
    double x{}, y{}, z{};
};

constexpr V3 operator+(const V3& a, const V3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
constexpr V3 operator-(const V3& a, const V3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
constexpr V3 operator*(double s, const V3& a) { return {s * a.x, s * a.y, s * a.z}; }

constexpr double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

constexpr V3 cross(const V3& a, const V3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double mag(const V3& a) { return std::sqrt(dot(a, a)); }

V3 unit(const V3& a) {
    const double m = mag(a);
    return {a.x / m, a.y / m, a.z / m};
}

/// Spherical -> Cartesian. theta/phi in DEGREES (this is a test helper; the
/// naming rule for radians does not bite because nothing here is radians).
V3 from_sph(double p, double theta_deg, double phi_deg) {
    const double t = theta_deg * pi0::kDegToRad;
    const double f = phi_deg * pi0::kDegToRad;
    return {p * std::sin(t) * std::cos(f), p * std::sin(t) * std::sin(f), p * std::cos(t)};
}

pi0::Photon photon(const V3& v, double e_gamma_deg = 30.0) {
    return pi0::Photon{v.x, v.y, v.z, e_gamma_deg};
}

pi0::Photon photon_sph(double e, double theta_deg, double phi_deg, double e_gamma_deg = 30.0) {
    return photon(from_sph(e, theta_deg, phi_deg), e_gamma_deg);
}

// ---------------------------------------------------------------------------
// Adapters onto the real core API.
// ---------------------------------------------------------------------------

/// compute_dis for a MASSLESS electron of energy `e_prime_gev` at polar angle
/// `theta_deg`: |p_e| = E', so the momentum vector and the energy argument
/// agree. The azimuth is arbitrary for every DIS variable and is fixed at 0.
///
/// This helper exists only to keep the DIS tests reading in (E', theta), which
/// is how the note states them. It deliberately COUPLES |p_e| to E'; the test
/// "compute_dis uses the energy it is given" below bypasses it to pin the fact
/// that the two are independent arguments.
pi0::DisKin dis_at(double e_prime_gev, double theta_deg, double beam_gev = pi0::kBeamEnergyGeV) {
    const V3 pe = from_sph(e_prime_gev, theta_deg, 0.0);
    return pi0::kin::compute_dis(pe.x, pe.y, pe.z, e_prime_gev, beam_gev);
}

/// compute_sidis for a massless scattered electron `ep` and a hadron `h` of
/// energy `h_e`. Threads the event DisKin through exactly as production must:
/// compute_sidis reads `dis.nu`, and nu = E_beam - E' with E' = |p_e| here.
pi0::SidisKin sidis_of(const V3& ep, const V3& h, double h_e, double beam_gev = pi0::kBeamEnergyGeV) {
    const double e_energy = mag(ep);  // massless electron
    const pi0::DisKin dis = pi0::kin::compute_dis(ep.x, ep.y, ep.z, e_energy, beam_gev);
    return pi0::kin::compute_sidis(h.x, h.y, h.z, h_e, dis, ep.x, ep.y, ep.z, e_energy, beam_gev);
}

/// The pairing cuts this analysis actually ships, transcribed from
/// config/cuts.json ("pairing" block -- the single source of truth for every
/// cut value in this project).
///
/// PairingCuts has no default member values on purpose: it is documented as
/// "Populated from cuts.json", and a struct default would be a second home for
/// a cut value, which is exactly what cuts.json exists to prevent. So the tests
/// must supply the numbers, and PairingCuts{} (all zeros, every cut disabled)
/// is never what a test wants.
pi0::PairingCuts shipped_cuts() {
    pi0::PairingCuts c{};
    c.mass_window_gev = 0.2;          // |m_gg - m_pi0| < 0.2  -> m_gg < 0.335, no lower bound
    c.min_mgg_gev = 0.001;            // m_gg > 0.001
    c.e_gamma_min_angle_deg = 8.0;    // theta_e_gamma > 8.0   -- NOT applied by find_gg_pairs
    c.open_a_deg = 17.561;            // theta_min(p) = a exp(-b p) + offset
    c.open_b_inv_gev = 0.756;
    c.open_offset_deg = 1.0;
    return c;
}

/// The single-photon e-gamma filter that Pairing.hpp REQUIRES the caller to
/// apply BEFORE find_gg_pairs (`Photon::e_gamma_deg > cuts.e_gamma_min_angle_deg`).
/// find_gg_pairs does not apply it -- see the e-gamma test at the bottom.
std::vector<pi0::Photon> filter_e_gamma(const std::vector<pi0::Photon>& in,
                                        const pi0::PairingCuts& cuts) {
    std::vector<pi0::Photon> out;
    for (const pi0::Photon& p : in) {
        if (p.e_gamma_deg > cuts.e_gamma_min_angle_deg) out.push_back(p);
    }
    return out;
}

/// Opening angle, in DEGREES, that two massless photons of energy e1,e2 must
/// subtend to give invariant mass m:  m^2 = 2 e1 e2 (1 - cos t12).
double theta_for_mass_deg(double m, double e1, double e2) {
    return std::acos(1.0 - m * m / (2.0 * e1 * e2)) * pi0::kRadToDeg;
}

/// Invariant mass of two massless photons, computed locally.
double mgg_local(const pi0::Photon& a, const pi0::Photon& b) {
    const double e = a.p() + b.p();
    const V3 s{a.px + b.px, a.py + b.py, a.pz + b.pz};
    const double m2 = e * e - dot(s, s);
    return m2 > 0.0 ? std::sqrt(m2) : 0.0;
}

/// The CANONICAL Trento definition, from scratch, with l = the INCOMING lepton:
///   phi_h = sign[(q x l) . P_h] * acos( ((q x l).(q x P_h)) / (|q x l||q x P_h|) )
/// Returned in DEGREES.
///
/// Aside, worth knowing before you "fix" an implementation that uses the
/// scattered lepton here: since q = l - l', we have
///     q x l = (l - l') x l = -l' x l = l x l'
///     q x l' = (l - l') x l' = l x l'
/// so q x l and q x l' are the SAME vector and the incoming/scattered choice
/// is not a real degree of freedom. What does matter, and what the tests below
/// pin, is the ORDER: cross(q, l) vs cross(l, q) flips yhat and hence phi_h.
double phi_h_canonical_deg(const V3& e_prime, const V3& p_h, double beam_gev) {
    const V3 l{0.0, 0.0, beam_gev};  // incoming lepton, along +z
    const V3 q = l - e_prime;

    const V3 qxl = cross(q, l);
    const V3 qxh = cross(q, p_h);

    double c = dot(qxl, qxh) / (mag(qxl) * mag(qxh));
    c = std::clamp(c, -1.0, 1.0);  // acos is undefined a hair outside [-1,1]

    const double s = (dot(qxl, p_h) >= 0.0) ? 1.0 : -1.0;
    return s * std::acos(c) * pi0::kRadToDeg;
}

/// Signed difference of two angles in DEGREES, wrapped into [-180, 180]. Needed
/// because phi_h = -180 and phi_h = +180 are the same direction: a raw
/// subtraction there reports a 360 degree "error" that does not exist.
double angle_diff_deg(double a, double b) { return std::remainder(a - b, 360.0); }

/// True if `p` joins photons i and j, in either field order.
bool joins(const pi0::GGPair& p, std::size_t i, std::size_t j) {
    return (p.i == i && p.j == j) || (p.i == j && p.j == i);
}

}  // namespace

// ===========================================================================
// Kinematics: compute_dis
// ===========================================================================

TEST_CASE("compute_dis reproduces a hand-computed DIS point", "[kin][dis]") {
    // E_beam = 10.53, E' = 4.0 GeV, theta_e = 20 degrees.
    //
    // Worked out offline to full double precision:
    //   Q2 = 4 * 10.53 * 4.0 * sin^2(10 deg)      = 5.0802936249950763
    //   nu = 10.53 - 4.0                          = 6.53
    //   y  = 6.53 / 10.53                         = 0.62013295346628683
    //   xb = Q2 / (2 * 0.938272 * 6.53)           = 0.41458814616740869
    //   W  = sqrt(Mp^2 + 2*Mp*nu - Q2)            = 2.8379381672243884
    const pi0::DisKin k = dis_at(4.0, 20.0);

    REQUIRE_THAT(k.q2, WithinRel(5.0802936249950763, 1e-12));
    REQUIRE_THAT(k.nu, WithinRel(6.53, 1e-12));
    REQUIRE_THAT(k.y, WithinRel(0.62013295346628683, 1e-12));
    REQUIRE_THAT(k.xb, WithinRel(0.41458814616740869, 1e-12));
    REQUIRE_THAT(k.w, WithinRel(2.8379381672243884, 1e-12));

    // ...and again from the defining formulae, spelled out here so this test
    // does not merely re-assert a number copied out of the implementation.
    const double half_rad = 0.5 * 20.0 * pi0::kDegToRad;
    const double q2 = 4.0 * pi0::kBeamEnergyGeV * 4.0 * std::sin(half_rad) * std::sin(half_rad);
    const double nu = pi0::kBeamEnergyGeV - 4.0;
    const double y = nu / pi0::kBeamEnergyGeV;
    const double xb = q2 / (2.0 * pi0::kProtonMassGeV * nu);
    const double w2 =
        pi0::kProtonMassGeV * pi0::kProtonMassGeV + 2.0 * pi0::kProtonMassGeV * nu - q2;

    CHECK_THAT(k.q2, WithinRel(q2, 1e-12));
    CHECK_THAT(k.nu, WithinRel(nu, 1e-12));
    CHECK_THAT(k.y, WithinRel(y, 1e-12));
    CHECK_THAT(k.xb, WithinRel(xb, 1e-12));
    CHECK_THAT(k.w, WithinRel(std::sqrt(w2), 1e-12));
}

TEST_CASE("compute_dis takes theta in DEGREES, not radians", "[kin][dis][units]") {
    // Analytic anchors that are unambiguous in degrees:
    //   theta = 180 deg -> sin^2(90 deg) = 1   -> Q2 = 4 E E'
    //   theta =  90 deg -> sin^2(45 deg) = 1/2 -> Q2 = 2 E E'
    // Were the angle silently treated as radians, sin^2(90 rad) = 0.79923...
    // and sin^2(45 rad) = 0.72659..., both grossly different.
    //
    // theta_e now reaches compute_dis as the direction of the electron momentum
    // vector rather than as a scalar in some unit, so a radian/degree mix-up can
    // no longer happen at THIS boundary -- the vector carries no unit to
    // confuse. The test is kept anyway: it still pins the sin^2(theta/2) form
    // and the half-angle identity the implementation uses in its place
    // (sin^2(theta/2) = (1 - cos theta)/2), which is where the same error would
    // now have to live.
    const double e_prime = 1.0;

    CHECK_THAT(dis_at(e_prime, 180.0).q2, WithinRel(4.0 * pi0::kBeamEnergyGeV * e_prime, 1e-12));
    CHECK_THAT(dis_at(e_prime, 90.0).q2, WithinRel(2.0 * pi0::kBeamEnergyGeV * e_prime, 1e-12));

    // theta = 0 is exactly forward: no momentum transfer at all.
    CHECK_THAT(dis_at(e_prime, 0.0).q2, WithinAbs(0.0, 1e-15));
}

TEST_CASE("compute_dis uses the energy it is given, not |p_e|", "[kin][dis][units]") {
    // Not reachable through dis_at, which ties |p_e| to E' for the massless
    // case. This pins the documented reason the two are separate arguments: a
    // caller passing a calorimeter-corrected E' must get THAT energy in Q2 and
    // nu, not one silently re-derived from the tracked momentum.
    //
    // Electron direction fixed at theta = 20 deg but with |p_e| = 8.0, while the
    // energy argument says E' = 4.0. Every DIS variable must come out exactly as
    // in the hand-computed point above, which used |p_e| = E' = 4.0: only the
    // DIRECTION of the vector may be read.
    const V3 pe = from_sph(8.0, 20.0, 0.0);
    const pi0::DisKin k = pi0::kin::compute_dis(pe.x, pe.y, pe.z, 4.0, pi0::kBeamEnergyGeV);

    CHECK_THAT(k.q2, WithinRel(5.0802936249950763, 1e-12));
    CHECK_THAT(k.nu, WithinRel(6.53, 1e-12));
    CHECK_THAT(k.xb, WithinRel(0.41458814616740869, 1e-12));
    CHECK_THAT(k.w, WithinRel(2.8379381672243884, 1e-12));
}

// ---------------------------------------------------------------------------
// The electron mass
// ---------------------------------------------------------------------------

TEST_CASE("energy_from_p carries the electron mass", "[kin][dis][mass]") {
    // The scattered electron's energy is E' = sqrt(p^2 + m_e^2), NOT |p|. The
    // change is tiny and therefore easy to make invisible -- to apply it and
    // have nothing anywhere show that it happened, or to revert it and have
    // nothing complain. So the SIZE of the shift is pinned, not merely its
    // sign: these assertions fail both if the mass is dropped and if a
    // different mass creeps in.
    constexpr double p = 5.0;
    const double m = pi0::kElectronMassGeV;
    const double e_prime = pi0::energy_from_p(p, m);

    SECTION("the value is the one the old code used") {
        // clas12/include/clas12/Constants.hpp ELECTRONMASS.
        CHECK_THAT(pi0::kElectronMassGeV, WithinAbs(0.00051099891, 1e-17));
    }

    SECTION("E' exceeds |p| by m_e^2 / 2p") {
        // sqrt(p^2 + m^2) = p * sqrt(1 + (m/p)^2) ~ p + m^2/(2p) for m << p.
        // At p = 5 GeV that is (0.000511)^2 / 10 = 2.61e-8 GeV.
        const double expected_excess = m * m / (2.0 * p);
        CHECK_THAT(expected_excess, WithinRel(2.6112e-8, 1e-4));  // the order, stated

        CHECK(e_prime > p);                       // nonzero, and the right sign
        CHECK(e_prime != p);                      // not lost to rounding
        CHECK_THAT(e_prime - p, WithinRel(expected_excess, 1e-6));
    }

    SECTION("Q2 shifts by the same relative amount, ~5e-9") {
        // Q2 = 4 E_beam E' sin^2(theta/2) is LINEAR in E', so dQ2/Q2 = dE'/E'
        // = m^2/(2p^2) exactly to leading order. The electron direction is held
        // fixed and only the energy argument changes, so this isolates the mass.
        const V3 pe = from_sph(p, 20.0, 0.0);
        const double q2_massless = pi0::kin::compute_dis(pe.x, pe.y, pe.z, p, pi0::kBeamEnergyGeV).q2;
        const double q2_massive =
            pi0::kin::compute_dis(pe.x, pe.y, pe.z, e_prime, pi0::kBeamEnergyGeV).q2;

        CHECK(q2_massive > q2_massless);
        CHECK(q2_massive != q2_massless);

        const double rel = (q2_massive - q2_massless) / q2_massless;
        CHECK_THAT(rel, WithinRel(m * m / (2.0 * p * p), 1e-6));
        CHECK_THAT(rel, WithinRel(5.2224e-9, 1e-4));  // the order, stated

        // Well clear of double rounding: ~2e7 ulps, so this is a real shift
        // being measured and not noise.
        CHECK(rel > 1e-12);
    }

    SECTION("nu falls by the same absolute amount") {
        // nu = E_beam - E', so carrying the mass LOWERS nu by exactly the excess.
        const V3 pe = from_sph(p, 20.0, 0.0);
        const double nu_massless = pi0::kin::compute_dis(pe.x, pe.y, pe.z, p, pi0::kBeamEnergyGeV).nu;
        const double nu_massive =
            pi0::kin::compute_dis(pe.x, pe.y, pe.z, e_prime, pi0::kBeamEnergyGeV).nu;
        CHECK(nu_massive < nu_massless);
        CHECK_THAT(nu_massless - nu_massive, WithinRel(m * m / (2.0 * p), 1e-6));
    }

    SECTION("photons are NOT routed through this") {
        // E_gamma = |p_gamma| is the analysis's convention (cuts.json /photon),
        // and a massless energy_from_p call is not the same expression: it is a
        // sqrt of a square. They agree here, which is exactly why this must be
        // stated rather than left as a thing someone might "unify".
        CHECK_THAT(pi0::energy_from_p(3.0, 0.0), WithinAbs(3.0, 1e-15));
    }
}

TEST_CASE("y = nu/E_beam lies in (0,1) for physical inputs", "[kin][dis]") {
    for (double e_prime = 0.25; e_prime < pi0::kBeamEnergyGeV; e_prime += 0.25) {
        for (double theta_deg : {5.0, 12.5, 20.0, 35.0}) {
            const pi0::DisKin k = dis_at(e_prime, theta_deg);
            INFO("E' = " << e_prime << " GeV, theta = " << theta_deg << " deg");
            CHECK(k.y > 0.0);
            CHECK(k.y < 1.0);
            CHECK_THAT(k.y, WithinRel(k.nu / pi0::kBeamEnergyGeV, 1e-14));
        }
    }
}

TEST_CASE("compute_dis guards the xb division when nu = 0", "[kin][dis][guard]") {
    // E' = E_beam -> nu = 0 -> xb = Q2 / 0. Must not be inf or NaN.
    const pi0::DisKin k = dis_at(pi0::kBeamEnergyGeV, 20.0);

    REQUIRE_THAT(k.nu, WithinAbs(0.0, 1e-15));
    CHECK(std::isfinite(k.xb));
    CHECK_FALSE(std::isnan(k.xb));
    CHECK_THAT(k.y, WithinAbs(0.0, 1e-15));

    // Approaching nu = 0 from below must stay finite too.
    for (double eps : {1e-6, 1e-9, 1e-12}) {
        const pi0::DisKin kk = dis_at(pi0::kBeamEnergyGeV - eps, 20.0);
        INFO("nu = " << eps);
        CHECK(std::isfinite(kk.xb));
    }
}

TEST_CASE("compute_dis returns W = 0 rather than NaN when W^2 <= 0", "[kin][dis][guard]") {
    // E' = 10.0 GeV at theta = 60 deg: Q2 = 105.3, nu = 0.53
    //   -> W^2 = Mp^2 + 2 Mp nu - Q2 = -103.425... < 0
    // This is unphysical, but the guard must hold: sqrt of a negative is NaN.
    const pi0::DisKin k = dis_at(10.0, 60.0);

    const double w2 = pi0::kProtonMassGeV * pi0::kProtonMassGeV +
                      2.0 * pi0::kProtonMassGeV * k.nu - k.q2;
    REQUIRE(w2 < 0.0);  // the premise of the test

    CHECK_FALSE(std::isnan(k.w));
    CHECK(std::isfinite(k.w));
    CHECK_THAT(k.w, WithinAbs(0.0, 1e-15));
}

// ===========================================================================
// Trento phi_h -- the one that matters
// ===========================================================================

TEST_CASE("compute_sidis phi_h stays inside [-180, 180] over a scan", "[kin][sidis][trento]") {
    double max_abs = 0.0;

    for (double theta_e : {5.0, 15.0, 30.0}) {
        for (double phi_e : {-150.0, -60.0, 0.0, 75.0, 170.0}) {
            const V3 ep = from_sph(4.0, theta_e, phi_e);
            for (double theta_h : {2.0, 10.0, 25.0, 50.0}) {
                for (double phi_h_lab = -180.0; phi_h_lab < 180.0; phi_h_lab += 15.0) {
                    const V3 h = from_sph(1.5, theta_h, phi_h_lab);
                    const pi0::SidisKin s = sidis_of(ep, h, mag(h));

                    INFO("theta_e=" << theta_e << " phi_e=" << phi_e << " theta_h=" << theta_h
                                    << " phi_h_lab=" << phi_h_lab);
                    CHECK(std::isfinite(s.phi_h_deg));
                    CHECK(s.phi_h_deg >= -180.0);
                    CHECK(s.phi_h_deg <= 180.0);
                    max_abs = std::max(max_abs, std::abs(s.phi_h_deg));
                }
            }
        }
    }

    // A radian-returning implementation would also satisfy [-180,180] -- it
    // would just never exceed pi. Pin the units by demanding real coverage.
    CHECK(max_abs > 90.0);
}

TEST_CASE("compute_sidis phi_h matches the canonical Trento definition away from the poles",
          "[kin][sidis][trento]") {
    // The frame construction (note eq:frame) is
    //     qhat = q/|q|,  yhat = (q x k_beam)/|q x k_beam|,  xhat = yhat x qhat
    //     phi_h = atan2(pperp.yhat, pperp.xhat)
    // and is only *claimed* to be Trento. Check it against the canonical
    // sign/acos form, computed from scratch above, over pseudo-random
    // configurations. FIXED SEED: this test is deterministic.
    //
    // THIS TEST IS RESTRICTED TO |phi_h| in [1, 179], AND THAT IS NOT A DODGE.
    // The arccos reference is ill-conditioned at its degenerate points: near
    // phi_h = 0 and +-180 an error eps in cos phi_h emerges as an error of
    // ~sqrt(2 eps) in phi_h. So an "agreement over all phi" figure measures the
    // REFERENCE and the sample size, not this code -- it is empirically
    // unbounded in N, degrading to ~1.3e-7 deg at N = 2e7 (Kinematics.hpp
    // documents the measurement). Excluding the poles is what makes the number
    // below a statement about compute_sidis.
    //
    // The poles are NOT left untested: "compute_sidis recovers a phi_h built
    // into the hadron by construction" below covers them, uniformly in phi_h,
    // without an arccos anywhere. The two tests are complementary and neither
    // is redundant -- this one verifies the CONVENTION against an independent
    // definition, that one bounds the NUMERICAL ERROR.
    std::mt19937_64 rng(20250716ULL);
    std::uniform_real_distribution<double> u_eprime(2.0, 10.0);
    std::uniform_real_distribution<double> u_theta_e(5.0, 35.0);
    std::uniform_real_distribution<double> u_phi(-180.0, 180.0);
    std::uniform_real_distribution<double> u_ph(0.3, 5.0);
    std::uniform_real_distribution<double> u_theta_h(1.0, 60.0);

    constexpr int kConfigs = 2000;  // >= 500 required
    double worst = 0.0;
    int compared = 0;

    for (int n = 0; n < kConfigs; ++n) {
        const V3 ep = from_sph(u_eprime(rng), u_theta_e(rng), u_phi(rng));
        const V3 h = from_sph(u_ph(rng), u_theta_h(rng), u_phi(rng));

        const double want = phi_h_canonical_deg(ep, h, pi0::kBeamEnergyGeV);

        // Skip where the REFERENCE, not this code, loses its conditioning.
        if (std::abs(want) < 1.0 || std::abs(want) > 179.0) continue;
        ++compared;

        const pi0::SidisKin s = sidis_of(ep, h, mag(h));

        const double d = std::abs(angle_diff_deg(s.phi_h_deg, want));
        worst = std::max(worst, d);

        INFO("config " << n << ": got " << s.phi_h_deg << " want " << want);
        REQUIRE(d < 1e-9);
    }

    // The pole exclusion must not have quietly emptied the scan.
    INFO("compared " << compared << " of " << kConfigs << " configs");
    CHECK(compared > kConfigs / 2);

    // Aggregate bound. Do NOT report this with WARN: Catch2's TAP reporter
    // emits warnings as lines meson's TAP parser rejects, which turns a fully
    // passing run into a test ERROR. INFO only prints when the CHECK fails.
    INFO("max |phi_h - canonical| over " << compared << " in-range configs = " << worst << " deg");
    CHECK(worst < 1e-9);
}

TEST_CASE("compute_sidis recovers a phi_h built into the hadron by construction",
          "[kin][sidis][trento]") {
    // The round trip the arccos cross-check above cannot do: give the hadron a
    // phi_h that is true BY CONSTRUCTION and check it comes back. Build
    //     p_h = pt (cos phi xhat + sin phi yhat) + pl qhat
    // in the frame Kinematics.hpp documents, and compute_sidis must return that
    // phi. No arccos is involved on either side, so the deviation here bounds
    // THIS function's error -- and it holds uniformly in phi_h, including at
    // 0 and +-180 where the arccos reference falls apart.
    //
    // Note precisely what this does and does not prove. The frame is built here
    // from the same definition the implementation uses, so this test CANNOT
    // catch a wrong convention (a flipped yhat, the scattered lepton in place of
    // the beam) -- it would follow the error. Catching that is the job of the
    // canonical cross-check above, which uses an independent definition. This
    // test is the numerical-conditioning half of the pair, and the pair is only
    // meaningful together.
    std::mt19937_64 rng(20250716ULL);
    std::uniform_real_distribution<double> u_eprime(2.0, 10.0);
    std::uniform_real_distribution<double> u_theta_e(5.0, 35.0);
    std::uniform_real_distribution<double> u_phi_e(-180.0, 180.0);
    std::uniform_real_distribution<double> u_phi_h(-180.0, 180.0);
    std::uniform_real_distribution<double> u_pt(0.05, 2.0);
    std::uniform_real_distribution<double> u_pl(-2.0, 5.0);

    constexpr int kConfigs = 2000;
    double worst = 0.0;
    double worst_near_pole = 0.0;
    int near_pole = 0;

    for (int n = 0; n < kConfigs; ++n) {
        const V3 ep = from_sph(u_eprime(rng), u_theta_e(rng), u_phi_e(rng));

        // The Trento frame, from Kinematics.hpp's definition, built locally.
        const V3 l{0.0, 0.0, pi0::kBeamEnergyGeV};
        const V3 q = l - ep;
        const V3 q_hat = unit(q);
        const V3 y_hat = unit(cross(q, l));
        const V3 x_hat = cross(y_hat, q_hat);

        const double phi_true = u_phi_h(rng);
        const double pt = u_pt(rng);
        const double pl = u_pl(rng);

        const double phi_rad = phi_true * pi0::kDegToRad;
        const V3 h = (pt * std::cos(phi_rad)) * x_hat + (pt * std::sin(phi_rad)) * y_hat +
                     pl * q_hat;

        const pi0::SidisKin s = sidis_of(ep, h, mag(h));

        const double d = std::abs(angle_diff_deg(s.phi_h_deg, phi_true));
        worst = std::max(worst, d);
        if (std::abs(phi_true) < 1.0 || std::abs(phi_true) > 179.0) {
            ++near_pole;
            worst_near_pole = std::max(worst_near_pole, d);
        }

        INFO("config " << n << ": built phi_h = " << phi_true << ", got " << s.phi_h_deg);
        REQUIRE(d < 1e-9);

        // pt was put in by construction; pt2 must come back as its square.
        CHECK_THAT(s.pt2, WithinRel(pt * pt, 1e-9));
    }

    // The poles are the whole point of this test -- make sure they were sampled.
    INFO("sampled " << near_pole << " configs with |phi_h| outside [1, 179]");
    CHECK(near_pole > 0);

    INFO("max |phi_h - constructed| over " << kConfigs << " configs = " << worst
                                           << " deg (near poles: " << worst_near_pole << ")");
    CHECK(worst < 1e-9);
}

TEST_CASE("compute_sidis phi_h ~ 0 in-plane on the scattered-electron side", "[kin][sidis][trento]") {
    // Scattered electron at theta = 20 deg, phi = 0 -> it lies on the +x side
    // of the beam. xhat is proportional to k_perp, i.e. it points to that same
    // side, so an in-plane hadron on the +x side has phi_h = 0 and one on the
    // -x side has phi_h = 180.
    //
    // (Kinematics.hpp derives xhat from the INCOMING lepton's transverse
    // component. That is the same side: q_perp = 0 by construction, and
    // q = k - k', so k_perp = k'_perp identically. Incoming or scattered is not
    // a degree of freedom here -- see phi_h_canonical_deg's note.)
    const V3 ep = from_sph(4.0, 20.0, 0.0);

    SECTION("same side as the scattered electron -> phi_h ~ 0") {
        for (double theta_h : {5.0, 20.0, 45.0}) {
            const V3 h = from_sph(2.0, theta_h, 0.0);  // +x side, in plane
            const pi0::SidisKin s = sidis_of(ep, h, mag(h));
            INFO("theta_h = " << theta_h);
            CHECK_THAT(s.phi_h_deg, WithinAbs(0.0, 1e-9));
        }
    }

    SECTION("opposite side -> phi_h ~ 180") {
        for (double theta_h : {20.0, 40.0, 60.0}) {
            const V3 h = from_sph(2.0, theta_h, 180.0);  // -x side, in plane
            const pi0::SidisKin s = sidis_of(ep, h, mag(h));
            INFO("theta_h = " << theta_h);
            // atan2(+0, negative) = +pi, but do not depend on the sign of zero.
            CHECK_THAT(std::abs(s.phi_h_deg), WithinAbs(180.0, 1e-9));
        }
    }
}

TEST_CASE("compute_sidis phi_h flips sign when the hadron's yhat component flips",
          "[kin][sidis][trento]") {
    // Reflecting the hadron through the lepton scattering plane sends
    // (pperp.yhat) -> -(pperp.yhat) and leaves (pperp.xhat) alone, so
    // phi_h -> -phi_h. Build the reflection from yhat itself rather than from
    // "flip py", so the test does not assume the electron sits at phi = 0.
    const struct {
        double theta_e, phi_e, theta_h, phi_h;
    } cases[] = {
        {20.0, 0.0, 15.0, 60.0},
        {20.0, 0.0, 15.0, -35.0},
        {8.0, 130.0, 30.0, 10.0},
        {33.0, -95.0, 5.0, 145.0},
    };

    for (const auto& c : cases) {
        const V3 ep = from_sph(4.0, c.theta_e, c.phi_e);
        const V3 h = from_sph(2.0, c.theta_h, c.phi_h);

        const V3 l{0.0, 0.0, pi0::kBeamEnergyGeV};
        const V3 q = l - ep;
        const V3 yh = unit(cross(q, l));

        // Mirror the hadron across the scattering plane: p -> p - 2(p.yhat)yhat
        const double b = dot(h, yh);
        const V3 h_mirror = h - (2.0 * b) * yh;

        const pi0::SidisKin s0 = sidis_of(ep, h, mag(h));
        const pi0::SidisKin s1 = sidis_of(ep, h_mirror, mag(h_mirror));

        INFO("theta_e=" << c.theta_e << " phi_e=" << c.phi_e << " theta_h=" << c.theta_h
                        << " phi_h=" << c.phi_h << " -> phi " << s0.phi_h_deg << " / "
                        << s1.phi_h_deg);

        // Guard against a degenerate config where the flip is a no-op.
        REQUIRE(std::abs(s0.phi_h_deg) > 1e-3);
        REQUIRE(std::abs(s0.phi_h_deg) < 180.0 - 1e-3);

        CHECK_THAT(s1.phi_h_deg, WithinAbs(-s0.phi_h_deg, 1e-9));

        // The reflection is an isometry about the q axis: z and pt2 are untouched.
        CHECK_THAT(s1.pt2, WithinRel(s0.pt2, 1e-12));
        CHECK_THAT(s1.z, WithinRel(s0.z, 1e-12));
    }
}

TEST_CASE("compute_sidis z and pt2 follow their definitions", "[kin][sidis]") {
    // z = E_h / nu ; pt2 = |p_h - (p_h.qhat) qhat|^2  (note eq:z, eq:pT2)
    const V3 ep = from_sph(4.0, 20.0, 35.0);
    const V3 h = from_sph(2.0, 12.0, -70.0);
    const double h_e = mag(h);

    const pi0::SidisKin s = sidis_of(ep, h, h_e);

    const double nu = pi0::kBeamEnergyGeV - mag(ep);
    CHECK_THAT(s.z, WithinRel(h_e / nu, 1e-12));

    const V3 l{0.0, 0.0, pi0::kBeamEnergyGeV};
    const V3 qh = unit(l - ep);
    const V3 pperp = h - dot(h, qh) * qh;
    CHECK_THAT(s.pt2, WithinRel(dot(pperp, pperp), 1e-12));

    // pt2 is a squared magnitude: never negative, and bounded by |p_h|^2.
    CHECK(s.pt2 >= 0.0);
    CHECK(s.pt2 <= dot(h, h) + 1e-12);
}

// ===========================================================================
// Pairing: min_opening_angle_deg
// ===========================================================================

TEST_CASE("min_opening_angle_deg follows 17.561*exp(-0.756*p) + 1.0", "[pairing][opening]") {
    // note eq:theta-min. Threshold is in DEGREES and p is the PAIR momentum.
    //
    // The three coefficients are cut values carried in PairingCuts, not
    // constants inside the function, so what is pinned here is that the shipped
    // cuts.json numbers produce the note's curve AND that each lands in its
    // documented role: open_a_deg is the amplitude, open_b_inv_gev the decay
    // constant, open_offset_deg the asymptote. Swapping any two would survive
    // the p = 0 anchor alone, which is why the scan and the asymptote are here.
    const pi0::PairingCuts cuts = shipped_cuts();

    SECTION("p = 0 gives a + offset = 18.561 deg") {
        CHECK_THAT(pi0::min_opening_angle_deg(0.0, cuts), WithinRel(18.561, 1e-12));
    }

    SECTION("the note's quoted waypoint: 4.87 deg at p = 2 GeV") {
        CHECK_THAT(pi0::min_opening_angle_deg(2.0, cuts), WithinRel(4.8716490763290485, 1e-12));
    }

    SECTION("matches the closed form across the range") {
        for (double p = 0.0; p <= 12.0; p += 0.1) {
            INFO("p = " << p);
            CHECK_THAT(pi0::min_opening_angle_deg(p, cuts),
                       WithinRel(17.561 * std::exp(-0.756 * p) + 1.0, 1e-12));
        }
    }

    SECTION("strictly decreasing in p") {
        double prev = pi0::min_opening_angle_deg(0.0, cuts);
        for (double p = 0.05; p <= 12.0; p += 0.05) {
            const double cur = pi0::min_opening_angle_deg(p, cuts);
            INFO("p = " << p << ": " << cur << " should be < " << prev);
            CHECK(cur < prev);
            prev = cur;
        }
    }

    SECTION("tends to the 1.0 deg offset at large p") {
        CHECK(pi0::min_opening_angle_deg(50.0, cuts) > 1.0);
        CHECK_THAT(pi0::min_opening_angle_deg(50.0, cuts), WithinAbs(1.0, 1e-9));
        CHECK_THAT(pi0::min_opening_angle_deg(1e3, cuts), WithinAbs(1.0, 1e-12));
    }
}

// ===========================================================================
// Pairing: find_gg_pairs
// ===========================================================================

TEST_CASE("find_gg_pairs returns nothing for degenerate inputs", "[pairing][gg]") {
    const pi0::PairingCuts cuts = shipped_cuts();

    SECTION("no photons") {
        const std::vector<pi0::Photon> none;
        CHECK(pi0::find_gg_pairs(none, cuts).empty());
    }

    SECTION("one photon cannot pair with itself") {
        const std::vector<pi0::Photon> one{photon_sph(1.0, 0.0, 0.0)};
        CHECK(pi0::find_gg_pairs(one, cuts).empty());
    }
}

TEST_CASE("find_gg_pairs recovers a pair sitting exactly on the pi0 mass", "[pairing][gg]") {
    // Two 1 GeV photons split symmetrically about +z by exactly the angle that
    // makes m_gg = m_pi0 (7.7395 deg), which clears theta_min(~1.998 GeV)=4.885.
    const double t = theta_for_mass_deg(pi0::kPi0MassGeV, 1.0, 1.0);
    const std::vector<pi0::Photon> photons{
        photon_sph(1.0, +0.5 * t, 0.0),
        photon_sph(1.0, -0.5 * t, 0.0),  // negative theta = phi flipped by 180
    };

    const auto pairs = pi0::find_gg_pairs(photons, shipped_cuts());

    REQUIRE(pairs.size() == 1);
    CHECK(joins(pairs[0], 0, 1));
    CHECK_THAT(pairs[0].mgg, WithinAbs(pi0::kPi0MassGeV, 1e-9));

    // The stored four-momentum must be the plain sum, massless photons.
    CHECK_THAT(pairs[0].e, WithinRel(2.0, 1e-12));
    CHECK_THAT(pairs[0].px, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pairs[0].py, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pairs[0].pz, WithinRel(2.0 * std::cos(0.5 * t * pi0::kDegToRad), 1e-12));

    // ...and mgg must be consistent with that four-momentum.
    const double m2 = pairs[0].e * pairs[0].e -
                      (pairs[0].px * pairs[0].px + pairs[0].py * pairs[0].py +
                       pairs[0].pz * pairs[0].pz);
    CHECK_THAT(pairs[0].mgg, WithinAbs(std::sqrt(m2), 1e-9));
}

TEST_CASE("find_gg_pairs is exclusive: no photon is used twice", "[pairing][gg][exclusivity]") {
    // Four photons in one cluster. All six pairs are admissible here, so a
    // non-exclusive implementation would happily emit more than two.
    const double t = theta_for_mass_deg(pi0::kPi0MassGeV, 1.0, 1.0);
    const std::vector<pi0::Photon> photons{
        photon_sph(1.0, 0.0, 0.0),
        photon_sph(1.0, t, 0.0),
        photon_sph(1.2, 9.0, 120.0),
        photon_sph(0.9, 11.0, 240.0),
    };

    const auto pairs = pi0::find_gg_pairs(photons, shipped_cuts());

    REQUIRE_FALSE(pairs.empty());               // else the test proves nothing
    CHECK(pairs.size() <= photons.size() / 2);  // floor(N/2) = 2

    std::vector<std::size_t> seen;
    for (const auto& p : pairs) {
        CHECK(p.i != p.j);
        CHECK(p.i < photons.size());
        CHECK(p.j < photons.size());
        seen.push_back(p.i);
        seen.push_back(p.j);
    }
    std::sort(seen.begin(), seen.end());
    CHECK(std::adjacent_find(seen.begin(), seen.end()) == seen.end());  // all distinct
}

TEST_CASE("find_gg_pairs is greedy best-first, not first-match", "[pairing][gg][greedy]") {
    // Photon 0 on +z; photons 1 and 2 each pair admissibly with it. The pair
    // closest to m_pi0 must win regardless of index order, and the loser's
    // partner must then be left unpaired (only one photon remains).
    const double t_pi0 = theta_for_mass_deg(pi0::kPi0MassGeV, 1.0, 1.0);  // 7.7395 deg
    const double t_far = theta_for_mass_deg(0.25, 1.0, 1.0);              // 14.3615 deg

    SECTION("(0,1) is the closer pair -> (0,1) wins, photon 2 unpaired") {
        const std::vector<pi0::Photon> photons{
            photon_sph(1.0, 0.0, 0.0),
            photon_sph(1.0, t_pi0, 0.0),   // m(0,1) = m_pi0      -> |dm| = 0
            photon_sph(1.0, t_far, 90.0),  // m(0,2) = 0.25       -> |dm| = 0.115
        };

        const auto pairs = pi0::find_gg_pairs(photons, shipped_cuts());

        REQUIRE(pairs.size() == 1);
        CHECK(joins(pairs[0], 0, 1));
        CHECK_THAT(pairs[0].mgg, WithinAbs(pi0::kPi0MassGeV, 1e-9));
    }

    SECTION("(0,2) is the closer pair -> (0,2) wins even though (0,1) comes first") {
        // Same geometry, roles swapped. A first-match implementation that
        // simply scans i<j and takes the first admissible pair would return
        // (0,1) here and fail. This is what makes the section above meaningful.
        const std::vector<pi0::Photon> photons{
            photon_sph(1.0, 0.0, 0.0),
            photon_sph(1.0, t_far, 0.0),   // m(0,1) = 0.25   -> |dm| = 0.115
            photon_sph(1.0, t_pi0, 90.0),  // m(0,2) = m_pi0  -> |dm| = 0
        };

        const auto pairs = pi0::find_gg_pairs(photons, shipped_cuts());

        REQUIRE(pairs.size() == 1);
        CHECK(joins(pairs[0], 0, 2));
        CHECK_THAT(pairs[0].mgg, WithinAbs(pi0::kPi0MassGeV, 1e-9));
    }
}

TEST_CASE("find_gg_pairs applies the opening-angle cut at the PAIR momentum",
          "[pairing][gg][opening]") {
    // Isolate theta_min: both pairs below clear the mass floor and the window,
    // so only the opening-angle cut can separate them.
    //
    //   two 1 GeV photons, 3 deg apart -> m = 0.0524, pair p = 1.9993
    //                                     theta_min(1.9993) = 4.8737 > 3  -> REJECT
    //   two 1 GeV photons, 8 deg apart -> m = 0.1395, pair p = 1.9951
    //                                     theta_min(1.9951) = 4.8859 < 8  -> KEEP
    const pi0::PairingCuts cuts = shipped_cuts();

    SECTION("merged-cluster pair (3 deg) is rejected") {
        const std::vector<pi0::Photon> photons{
            photon_sph(1.0, +1.5, 0.0),
            photon_sph(1.0, -1.5, 0.0),
        };
        const double m = mgg_local(photons[0], photons[1]);
        REQUIRE(m > cuts.min_mgg_gev);                                   // floor passes
        REQUIRE(std::abs(m - pi0::kPi0MassGeV) < cuts.mass_window_gev);  // window passes
        CHECK(pi0::find_gg_pairs(photons, cuts).empty());
    }

    SECTION("resolved pair (8 deg) is kept") {
        const std::vector<pi0::Photon> photons{
            photon_sph(1.0, +4.0, 0.0),
            photon_sph(1.0, -4.0, 0.0),
        };
        CHECK(pi0::find_gg_pairs(photons, cuts).size() == 1);
    }

    SECTION("the threshold is evaluated at the pair momentum, not a photon's") {
        // Same 3 deg opening, but each photon carries 5 GeV, so the PAIR
        // momentum is ~10 GeV and theta_min(10) = 1.009 deg -> now admissible.
        // Evaluating theta_min at a single photon's momentum (5 GeV,
        // theta_min = 1.40 deg) would also admit it, so pair the check with
        // the low-momentum case above, where the two disagree sharply
        // (theta_min(2.0) = 4.87 vs theta_min(1.0) = 9.25).
        const std::vector<pi0::Photon> photons{
            photon_sph(5.0, +1.5, 0.0),
            photon_sph(5.0, -1.5, 0.0),
        };
        REQUIRE(pi0::min_opening_angle_deg(10.0, cuts) < 3.0);
        CHECK(pi0::find_gg_pairs(photons, cuts).size() == 1);
    }
}

TEST_CASE("find_gg_pairs mass window is centred on m_pi0 and has NO lower bound",
          "[pairing][gg][window]") {
    // note sec:pairing: the shipped window is +-0.2 GeV about m_pi0, so
    //     |m_gg - 0.1349768| < 0.2   <=>   m_gg < 0.335 GeV
    // with no lower bound whatsoever. This is deliberate -- the signal is
    // fitted out of the m_gg spectrum downstream and needs its sidebands.
    //
    // Guard against the tempting misreading `m_gg < mass_window_gev`, which
    // would silently amputate the upper sideband at 0.2 GeV.
    const pi0::PairingCuts cuts = shipped_cuts();

    SECTION("a 0.30 GeV pair is INSIDE the window and must be kept") {
        const double t = theta_for_mass_deg(0.30, 1.0, 1.0);  // 17.25 deg
        const std::vector<pi0::Photon> photons{
            photon_sph(1.0, +0.5 * t, 0.0),
            photon_sph(1.0, -0.5 * t, 0.0),
        };
        REQUIRE(0.30 > cuts.mass_window_gev);  // the misreading would reject it
        const auto pairs = pi0::find_gg_pairs(photons, cuts);
        REQUIRE(pairs.size() == 1);
        CHECK_THAT(pairs[0].mgg, WithinAbs(0.30, 1e-9));
    }

    SECTION("a 0.34 GeV pair is OUTSIDE the window and must be dropped") {
        const double t = theta_for_mass_deg(0.34, 1.0, 1.0);  // 19.58 deg
        const std::vector<pi0::Photon> photons{
            photon_sph(1.0, +0.5 * t, 0.0),
            photon_sph(1.0, -0.5 * t, 0.0),
        };
        CHECK(pi0::find_gg_pairs(photons, cuts).empty());
    }
}

TEST_CASE("find_gg_pairs leaves photons of a window-rejected pair AVAILABLE",
          "[pairing][gg][window]") {
    // This pins the deliberate semantic change vs the old code: failing the
    // mass window removes the PAIR, not its photons.
    //
    //   photon 0 on +z
    //   photon 1 at 28.955 deg, phi = 0   -> m(0,1) = 0.500  -> |dm| = 0.365 > 0.2  REJECTED
    //   photon 2 at  7.740 deg, phi = 180 -> m(0,2) = m_pi0  -> |dm| = 0            OK
    //   m(1,2) = 0.6296                   -> |dm| = 0.495 > 0.2                     REJECTED
    //
    // Photon 0 appears in a rejected pair, yet must still be consumed by (0,2).
    const double t_pi0 = theta_for_mass_deg(pi0::kPi0MassGeV, 1.0, 1.0);
    const double t_out = theta_for_mass_deg(0.500, 1.0, 1.0);

    const std::vector<pi0::Photon> photons{
        photon_sph(1.0, 0.0, 0.0),
        photon_sph(1.0, t_out, 0.0),
        photon_sph(1.0, t_pi0, 180.0),
    };

    // Premise: (0,1) really is outside the window.
    REQUIRE(std::abs(mgg_local(photons[0], photons[1]) - pi0::kPi0MassGeV) > 0.2);

    const auto pairs = pi0::find_gg_pairs(photons, shipped_cuts());

    REQUIRE(pairs.size() == 1);
    CHECK(joins(pairs[0], 0, 2));  // photon 0 survived the rejected pair (0,1)
    CHECK_THAT(pairs[0].mgg, WithinAbs(pi0::kPi0MassGeV, 1e-9));

    // No pair may carry photon 1: nothing admissible remains for it.
    for (const auto& p : pairs) {
        CHECK(p.i != 1u);
        CHECK(p.j != 1u);
    }
}

TEST_CASE("find_gg_pairs rejects collinear photons via the mass floor", "[pairing][gg][mass-floor]") {
    // Exactly collinear -> theta_12 = 0 -> m_gg = 0. Note that the mass WINDOW
    // does not catch this: |0 - 0.135| = 0.135 < 0.2. Only min_mgg_gev does.
    const pi0::PairingCuts cuts = shipped_cuts();
    const std::vector<pi0::Photon> photons{
        photon(V3{0.0, 0.0, 1.0}),
        photon(V3{0.0, 0.0, 2.0}),
    };

    CHECK_THAT(mgg_local(photons[0], photons[1]), WithinAbs(0.0, 1e-12));
    // The window would have let it through -- this is the point of the floor.
    CHECK(std::abs(mgg_local(photons[0], photons[1]) - pi0::kPi0MassGeV) < cuts.mass_window_gev);

    CHECK(pi0::find_gg_pairs(photons, cuts).empty());

    SECTION("nearly collinear is rejected too") {
        const std::vector<pi0::Photon> near{
            photon_sph(1.0, 0.0, 0.0),
            photon_sph(2.0, 0.02, 0.0),  // 0.02 deg apart -> m ~ 5e-4 GeV
        };
        CHECK(mgg_local(near[0], near[1]) < cuts.min_mgg_gev);
        CHECK(pi0::find_gg_pairs(near, cuts).empty());
    }
}

TEST_CASE("find_gg_pairs: min_mgg_gev is wired up independently of the opening angle",
          "[pairing][gg][mass-floor]") {
    // The collinear test above does NOT isolate the mass floor: a collinear
    // pair also fails the opening-angle cut, so deleting the floor entirely
    // leaves that test green. Worse, for any *physical* photon (E >= 0.2 GeV,
    // note tab:photon-cuts) the opening-angle cut already implies
    // m_gg >= ~0.049 GeV, ~50x the floor -- the floor cannot fire in
    // production at all.
    //
    // So to prove the floor is actually connected we need a synthetic pair,
    // with energies far below the photon threshold, that clears the opening
    // angle and is rejected by the floor alone:
    //
    //   E1 = E2 = 0.003 GeV, opening 19 deg
    //     pair p    = 0.005918 GeV -> theta_min = 18.4826 deg  -> 19 > 18.48  PASSES
    //     m_gg      = 0.00099029 GeV                           -> < 0.001     REJECTED
    //     |m - m_pi0| = 0.134 < 0.2                            -> window      PASSES
    //
    // If this ever starts failing because the floor was removed as dead code,
    // that is a legitimate decision -- delete this test with it, deliberately.
    const double e = 0.003;
    const double t = 19.0;
    const std::vector<pi0::Photon> photons{
        photon_sph(e, +0.5 * t, 0.0),
        photon_sph(e, -0.5 * t, 0.0),
    };

    const pi0::PairingCuts cuts = shipped_cuts();

    // Premises: only the floor should be able to reject this pair.
    const double m = mgg_local(photons[0], photons[1]);
    REQUIRE(m < cuts.min_mgg_gev);                                  // floor rejects
    REQUIRE(std::abs(m - pi0::kPi0MassGeV) < cuts.mass_window_gev);  // window does not
    const double p_pair = mag(V3{photons[0].px + photons[1].px, photons[0].py + photons[1].py,
                                 photons[0].pz + photons[1].pz});
    REQUIRE(t > pi0::min_opening_angle_deg(p_pair, cuts));  // opening angle does not

    CHECK(pi0::find_gg_pairs(photons, cuts).empty());

    // Lowering the floor below this pair's mass -- and changing NOTHING else --
    // must let the pair through. That pins the rejection to the floor rather
    // than to some other cut that happens to fire on the same pair.
    pi0::PairingCuts loose = shipped_cuts();
    loose.min_mgg_gev = 1e-6;
    CHECK(pi0::find_gg_pairs(photons, loose).size() == 1);
}

TEST_CASE("the photon-electron angle cut lives UPSTREAM of find_gg_pairs",
          "[pairing][gg][e-gamma]") {
    // note tab:pair-cuts: theta_e_gamma > 8.0 deg, to kill bremsstrahlung
    // photons sitting on top of the scattered electron.
    //
    // API RECONCILIATION, and a real change of contract rather than a rename:
    // this test used to assert that find_gg_pairs applied this cut. IT DOES
    // NOT, deliberately and documented (Pairing.hpp): theta_e_gamma is a
    // property of a SINGLE photon, not of a pair, so find_gg_pairs ignores
    // Photon::e_gamma_deg entirely and the caller must filter the photon list
    // on `e_gamma_deg > cuts.e_gamma_min_angle_deg` first.
    // PairingCuts::e_gamma_min_angle_deg is carried only so the pairing block
    // of cuts.json maps onto one struct.
    //
    // The original intent -- a brems photon must never reach a pi0 candidate,
    // and BOTH photons of a pair must be tested, strictly -- is preserved
    // below by testing the documented pipeline (filter, then pair) instead of
    // find_gg_pairs alone.
    //
    // KNOWN GAP, worth stating rather than hiding: `filter_e_gamma` is a TEST
    // helper. The library ships no upstream filter and there is no production
    // caller yet (src/tools holds only dump_columns), so nothing in src/
    // currently applies this cut. These sections pin the recipe and the
    // boundary; they cannot pin an implementation that does not exist. When the
    // real photon-selection step lands, point these sections at it.
    const double t = theta_for_mass_deg(pi0::kPi0MassGeV, 1.0, 1.0);
    const pi0::PairingCuts cuts = shipped_cuts();

    SECTION("find_gg_pairs itself ignores e_gamma_deg -- the documented contract") {
        // Both photons sit right on top of the electron. Handed straight to
        // find_gg_pairs, with no upstream filter, the pair MUST still come out:
        // if this ever starts failing, someone has moved the cut inside, and
        // then the cut has two homes and the cutflow can double-count it.
        const std::vector<pi0::Photon> unfiltered{
            photon_sph(1.0, +0.5 * t, 0.0, 0.5),  // 0.5 deg from the electron
            photon_sph(1.0, -0.5 * t, 0.0, 0.5),
        };
        CHECK(pi0::find_gg_pairs(unfiltered, cuts).size() == 1);
    }

    SECTION("both photons well clear of the electron -> pair survives") {
        const std::vector<pi0::Photon> photons{
            photon_sph(1.0, +0.5 * t, 0.0, 30.0),
            photon_sph(1.0, -0.5 * t, 0.0, 30.0),
        };
        CHECK(pi0::find_gg_pairs(filter_e_gamma(photons, cuts), cuts).size() == 1);
    }

    // Both operands must be checked. Testing only one side lets a filter that
    // examines a single photon of the pair survive.
    SECTION("SECOND photon too close to the electron -> pair dropped") {
        const std::vector<pi0::Photon> photons{
            photon_sph(1.0, +0.5 * t, 0.0, 30.0),
            photon_sph(1.0, -0.5 * t, 0.0, 2.0),  // 2 deg < 8 deg
        };
        CHECK(pi0::find_gg_pairs(filter_e_gamma(photons, cuts), cuts).empty());
    }

    SECTION("FIRST photon too close to the electron -> pair dropped") {
        const std::vector<pi0::Photon> photons{
            photon_sph(1.0, +0.5 * t, 0.0, 2.0),  // 2 deg < 8 deg
            photon_sph(1.0, -0.5 * t, 0.0, 30.0),
        };
        CHECK(pi0::find_gg_pairs(filter_e_gamma(photons, cuts), cuts).empty());
    }

    SECTION("the cut is strict at exactly 8 deg") {
        // note tab:pair-cuts: "All comparisons are strict" -> 8.0 is NOT kept.
        const std::vector<pi0::Photon> at_edge{
            photon_sph(1.0, +0.5 * t, 0.0, 8.0),
            photon_sph(1.0, -0.5 * t, 0.0, 30.0),
        };
        CHECK(pi0::find_gg_pairs(filter_e_gamma(at_edge, cuts), cuts).empty());

        const std::vector<pi0::Photon> just_over{
            photon_sph(1.0, +0.5 * t, 0.0, 8.0 + 1e-9),
            photon_sph(1.0, -0.5 * t, 0.0, 30.0),
        };
        CHECK(pi0::find_gg_pairs(filter_e_gamma(just_over, cuts), cuts).size() == 1);
    }
}
