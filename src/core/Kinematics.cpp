#include "core/Kinematics.hpp"

#include <algorithm>
#include <cmath>

#include "core/Constants.hpp"

namespace pi0::kin {
namespace {

/// Minimal three-vector helpers, local to this translation unit. We do not want
/// a ROOT TVector3 in src/core (it would drag ROOT into the unit tests), and a
/// general-purpose vector class would be more machinery than five functions.
struct Vec3 {
    double x{}, y{}, z{};
};

[[nodiscard]] constexpr double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

[[nodiscard]] constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

[[nodiscard]] inline double norm(const Vec3& a) {
    return std::sqrt(dot(a, a));
}

/// \return `a` scaled to unit length, or {0,0,0} if it has no length.
///         The caller is expected to have checked; this is belt and braces.
[[nodiscard]] inline Vec3 unit(const Vec3& a) {
    const double n = norm(a);
    if (!(n > 0.0)) return Vec3{};
    return Vec3{a.x / n, a.y / n, a.z / n};
}

/// Clamp to the valid domain of acos. A dot/norm quotient that is
/// mathematically 1 routinely comes out as 1 + 2e-16 in floating point, and
/// std::acos of that is a NaN.
[[nodiscard]] inline double clamped_acos(double c) {
    return std::acos(std::clamp(c, -1.0, 1.0));
}

}  // namespace

double angle_between_deg(double px1, double py1, double pz1, double px2, double py2, double pz2) {
    const Vec3 a{px1, py1, pz1};
    const Vec3 b{px2, py2, pz2};

    const double na = norm(a);
    const double nb = norm(b);

    // A zero-length vector has no direction, so the angle is undefined.
    // Return 0 rather than NaN: a NaN silently fails every comparison a caller
    // might write, so a cut like `theta > min` would drop the event with no
    // trace of why.
    if (!(na > 0.0) || !(nb > 0.0)) return 0.0;

    return clamped_acos(dot(a, b) / (na * nb)) * kRadToDeg;
}

DisKin compute_dis(double ex, double ey, double ez, double e_energy, double beam_energy) {
    DisKin k{};

    const Vec3 pe{ex, ey, ez};
    const double pe_mag = norm(pe);

    // theta_e w.r.t. the +z beam axis. With no momentum there is no direction,
    // so fall back to theta_e = 0, which gives Q2 = 0.
    const double cos_theta = (pe_mag > 0.0) ? std::clamp(ez / pe_mag, -1.0, 1.0) : 1.0;

    // Q2 = 4 E_beam E' sin^2(theta/2), with the half-angle identity
    // sin^2(theta/2) = (1 - cos theta) / 2 so we never round-trip through
    // acos and back.
    const double sin2_half = 0.5 * (1.0 - cos_theta);
    k.q2 = 4.0 * beam_energy * e_energy * sin2_half;

    k.nu = beam_energy - e_energy;

    // xb = Q2 / (2 M_p nu). See the header for why this is the PROTON mass on
    // every target. nu <= 0 is unphysical for DIS (the electron would have
    // gained energy); guard rather than divide.
    k.xb = (k.nu > 0.0) ? k.q2 / (2.0 * kProtonMassGeV * k.nu) : 0.0;

    // W^2 = (q + P_target)^2 = M_p^2 + 2 M_p nu - Q2, target at rest.
    const double w2 = kProtonMassGeV * kProtonMassGeV + 2.0 * kProtonMassGeV * k.nu - k.q2;
    k.w = (w2 > 0.0) ? std::sqrt(w2) : 0.0;

    k.y = (beam_energy > 0.0) ? k.nu / beam_energy : 0.0;

    return k;
}

SidisKin compute_sidis(double hx, double hy, double hz, double he, const DisKin& dis, double ex, double ey, double ez,
                       [[maybe_unused]] double e_energy, double beam_energy) {
    // e_energy is intentionally unused: the gamma* frame is fixed by the beam
    // and scattered-electron THREE-momenta, and z needs nu, which the caller
    // already computed in `dis`. The parameter stays in the signature so this
    // function reads the same way as compute_dis at every call site.
    SidisKin s{};

    const Vec3 k_beam{0.0, 0.0, beam_energy};  // beam along +z, massless
    const Vec3 k_e{ex, ey, ez};
    const Vec3 p_h{hx, hy, hz};

    // Virtual-photon three-momentum.
    const Vec3 q{k_beam.x - k_e.x, k_beam.y - k_e.y, k_beam.z - k_e.z};
    const double q_mag = norm(q);

    // No virtual photon means no frame to measure anything in.
    if (!(q_mag > 0.0)) return s;

    const Vec3 q_hat{q.x / q_mag, q.y / q_mag, q.z / q_mag};

    // z = E_h / nu. nu <= 0 is unphysical; guard rather than divide.
    s.z = (dis.nu > 0.0) ? he / dis.nu : 0.0;

    // Hadron momentum transverse to the virtual photon.
    const double p_long = dot(p_h, q_hat);
    const Vec3 p_perp{p_h.x - p_long * q_hat.x, p_h.y - p_long * q_hat.y, p_h.z - p_long * q_hat.z};
    s.pt2 = dot(p_perp, p_perp);

    // Trento frame. y_hat is normal to the lepton scattering plane.
    const Vec3 y_raw = cross(q, k_beam);
    const double y_mag = norm(y_raw);

    // y_mag == 0 means the electron scattered exactly along the beam axis, so
    // the lepton plane -- and with it phi_h -- is undefined. z and pt2 above
    // are still perfectly well defined, so keep them and report phi_h = 0.
    if (!(y_mag > 0.0)) {
        s.phi_h_deg = 0.0;
        return s;
    }

    const Vec3 y_hat = unit(y_raw);
    const Vec3 x_hat = cross(y_hat, q_hat);  // already unit: y_hat _|_ q_hat

    // atan2 handles p_perp == 0 by returning 0, which is the right answer for
    // "no transverse component, so no azimuth".
    s.phi_h_deg = std::atan2(dot(p_perp, y_hat), dot(p_perp, x_hat)) * kRadToDeg;

    return s;
}

}  // namespace pi0::kin
