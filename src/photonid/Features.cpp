#include "photonid/Features.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

#include "core/Constants.hpp"  // M_PI and every particle mass get_mass() returns

namespace pi0::photonid {

namespace {

// ---------------------------------------------------------------------------
// Model-structural constants.
//
// These are NOT cut values and do NOT belong in cuts.json. They are part of the
// definition of the 45 features the CatBoost models were trained on -- as fixed
// as the number 45 itself. Changing one does not retune a selection, it makes
// the feature vector mean something different from what the model expects.
// ---------------------------------------------------------------------------

/// Number of neighbour slots per species (m_g / m_ch / m_nh in the original).
constexpr std::size_t kNumPhotonNeighbours = 3;
constexpr std::size_t kNumChargedNeighbours = 2;
constexpr std::size_t kNumNeutralNeighbours = 2;

/// R windows, RADIANS, for the num_photons_0_* counters. The values are spelled
/// into the feature names themselves (num_photons_0_1 == "R < 0.1 rad").
constexpr double kPhotonCountR1 = 0.1;
constexpr double kPhotonCountR2 = 0.2;
constexpr double kPhotonCountR35 = 0.35;

// ---------------------------------------------------------------------------
// Particle masses: core/Constants.hpp, ALL of them.
//
// Ports physics::ParticleUtils::get_mass + Constants.hpp from the old tree.
// This file used to keep the electron, pion, kaon and neutron masses in this
// anonymous namespace because core/Constants.hpp carried only the proton and
// pi0 -- i.e. a second definition of four constants, which is the failure
// core/Constants.hpp's preamble exists to forbid. They now live there with the
// other two, including the deliberate 0.1396 pion wart, and get_mass() below
// reads them by unqualified lookup out of the enclosing pi0 namespace.
// ---------------------------------------------------------------------------

/// \return the mass in GeV, or NaN for a pid the old code did not know.
///
/// The NaN return is the mechanism, not an accident: `filter_photon` skips any
/// neighbour whose mass is NaN, which is how "unknown pid" is spelled.
double get_mass(int pid) {
    switch (pid) {
        case 11:
        case -11:
            return kElectronMassGeV;
        case 22:
            return 0.0;
        case 211:
        case -211:
            return kPionMassGeV;
        case 321:
        case -321:
            return kKaonMassGeV;
        case 2212:
        case -2212:
            return kProtonMassGeV;
        case 2112:
        case -2112:
            return kNeutronMassGeV;
        case 111:
            return kPi0MassGeV;
        default:
            return std::numeric_limits<double>::quiet_NaN();
    }
}

/// A bare 3-vector. Deliberately not ROOT::Math::XYZVector: the only thing the
/// old code used that class for was VectorUtil::Angle, reproduced below.
struct Vec3 {
    double x{0.0}, y{0.0}, z{0.0};
};

/// The calorimeter hit position used for the neighbour angles.
///
/// BIT-IDENTICAL: priority is PCAL -> ECIN -> ECOUT, decided by testing
/// `x == 0`, exactly as `get_particle_calo_vector` did. This is NOT the same as
/// "first layer present in the map": a real cluster sitting at x == 0.0 falls
/// through to the next layer, and a particle with no layers at all yields
/// (0,0,0). Both behaviours are inherited on purpose.
Vec3 calo_position(const CaloRowData& calo) {
    if (calo.pcal.x == 0) {
        if (calo.ecin.x == 0) {
            return Vec3{calo.ecout.x, calo.ecout.y, calo.ecout.z};
        }
        return Vec3{calo.ecin.x, calo.ecin.y, calo.ecin.z};
    }
    return Vec3{calo.pcal.x, calo.pcal.y, calo.pcal.z};
}

/// The 3D angle between two vectors, RADIANS.
///
/// BIT-IDENTICAL: this is a verbatim transcription of
/// ROOT::Math::VectorUtil::Angle + ::CosTheta for double-precision vectors
/// (ROOT 6.40, Math/GenVector/VectorUtil.h lines 143-173), which is what the
/// old code called on ROOT::Math::XYZVector. It is reproduced rather than
/// called so that the operation order -- and the clamp, and the ptot2 <= 0
/// branch that returns acos(0) = pi/2 for a zero-length vector -- are pinned
/// here and cannot drift with a ROOT upgrade.
///
/// core/Kinematics.hpp's angle_between_deg is deliberately NOT used: it returns
/// degrees, and a rad->deg->rad round trip is not bit-preserving.
double angle_between_rad(const Vec3& v1, const Vec3& v2) {
    double arg = 0.0;
    const double v1_r2 = v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
    const double v2_r2 = v2.x * v2.x + v2.y * v2.y + v2.z * v2.z;
    const double ptot2 = v1_r2 * v2_r2;
    if (ptot2 <= 0) {
        arg = 0.0;
    } else {
        const double pdot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
        arg = pdot / std::sqrt(ptot2);
        if (arg > 1.0) arg = 1.0;
        if (arg < -1.0) arg = -1.0;
    }
    return std::acos(arg);
}

/// The PID-purity filter, on a RADIAN angle.
///
/// BIT-IDENTICAL: the degree conversion is spelled `theta_rad * 180.0 / M_PI`
/// -- two operations -- and NOT `theta_rad * kRadToDeg`, which folds the
/// division into the constant and can land one ulp away. The result feeds a
/// `>=` against the theta window, so one ulp is the difference between keeping
/// and dropping a photon that sits exactly on 5.000 deg. This is the only place
/// in the project that may spell a rad->deg conversion by hand.
bool pass_pid_purity_filter(double e_gev, double epcal_gev, double theta_rad, const FeatureCuts& cuts) {
    const double theta_deg = theta_rad * 180.0 / M_PI;
    return e_gev >= cuts.min_energy_gev && epcal_gev > 0 && theta_deg >= cuts.theta_min_deg &&
           theta_deg <= cuts.theta_max_deg;
}

/// The N nearest neighbours of one species, kept sorted by R ascending.
///
/// BIT-IDENTICAL: `insert` is the original's hand-rolled insertion sort, quirks
/// intact. Note `R[i] == 0` doubles as "slot empty", so a genuine neighbour at
/// exactly R == 0 (a particle sharing the candidate's calorimeter position)
/// never displaces anything and is dropped. That is what the models were
/// trained against.
template <std::size_t N>
struct NeighbourSlots {
    double R[N]{};
    double dE[N]{};
    double Epcal[N]{};
    double m2u[N]{};
    double m2v[N]{};

    void insert(double r, double de, double epcal, double mu, double mv) {
        for (std::size_t i = 0; i < N; ++i) {
            if (r < R[i] || R[i] == 0) {
                for (std::size_t j = N - 1; j > i; --j) {
                    R[j] = R[j - 1];
                    dE[j] = dE[j - 1];
                    Epcal[j] = Epcal[j - 1];
                    m2u[j] = m2u[j - 1];
                    m2v[j] = m2v[j - 1];
                }
                R[i] = r;
                dE[i] = de;
                Epcal[i] = epcal;
                m2u[i] = mu;
                m2v[i] = mv;
                return;
            }
        }
    }
};

bool is_charged_hadron(int pid) {
    return pid == 2212 || pid == -2212 || pid == 211 || pid == -211 || pid == 321 || pid == -321;
}

bool is_neutral_hadron(int pid) {
    return pid == 2112 || pid == -2112;
}

/// Check the momentum columns line up, and that cand_row is addressable.
void check_momentum_columns(std::size_t cand_row, std::size_t n_px, std::size_t n_py, std::size_t n_pz) {
    if (n_px != n_py || n_px != n_pz) {
        throw std::invalid_argument(
            "pi0::photonid: REC::Particle momentum columns disagree in length (px=" +
            std::to_string(n_px) + ", py=" + std::to_string(n_py) + ", pz=" + std::to_string(n_pz) + ")");
    }
    if (cand_row >= n_px) {
        throw std::invalid_argument(
            "pi0::photonid: candidate row " + std::to_string(cand_row) + " is out of range for " +
            std::to_string(n_px) + " REC::Particle rows");
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// CaloMap
// ---------------------------------------------------------------------------

CaloMap CaloMap::build(
    const ROOT::VecOps::RVec<short>& pindex,
    const ROOT::VecOps::RVec<short>& layer,
    const ROOT::VecOps::RVec<short>& sector,
    const ROOT::VecOps::RVec<float>& energy,
    const ROOT::VecOps::RVec<float>& x,
    const ROOT::VecOps::RVec<float>& y,
    const ROOT::VecOps::RVec<float>& z,
    const ROOT::VecOps::RVec<float>& m2u,
    const ROOT::VecOps::RVec<float>& m2v,
    const ROOT::VecOps::RVec<float>& lu,
    const ROOT::VecOps::RVec<float>& lv,
    const ROOT::VecOps::RVec<float>& lw) {
    const std::size_t n = pindex.size();
    const bool ragged = layer.size() != n || sector.size() != n || energy.size() != n || x.size() != n ||
                        y.size() != n || z.size() != n || m2u.size() != n || m2v.size() != n ||
                        lu.size() != n || lv.size() != n || lw.size() != n;
    if (ragged) {
        throw std::invalid_argument(
            "pi0::photonid::CaloMap::build: REC::Calorimeter columns disagree in length (pindex=" +
            std::to_string(n) + ")");
    }

    CaloMap map;
    for (std::size_t row = 0; row < n; ++row) {
        CaloLayerData data;
        data.e = energy[row];
        data.x = x[row];
        data.y = y[row];
        data.z = z[row];
        data.m2u = m2u[row];
        data.m2v = m2v[row];
        data.lu = lu[row];
        data.lv = lv[row];
        data.lw = lw[row];
        data.sector = static_cast<int>(sector[row]);

        // operator[] default-constructs on first sight of this pindex, which is
        // the old code's `if (find == end) calo_map[pindex] = CaloRowData();`.
        CaloRowData& target = map.m_rows[static_cast<int>(pindex[row])];

        switch (static_cast<int>(layer[row])) {
            case kLayerPcal:
                target.pcal = data;
                break;
            case kLayerEcin:
                target.ecin = data;
                break;
            case kLayerEcout:
                target.ecout = data;
                break;
            default:
                // Unknown layer: the pindex is still registered (the entry was
                // created above), matching the old code, which also created the
                // entry before switching on the layer. That matters -- a
                // particle with only an unknown layer is "has calo data, no
                // PCAL", not "absent".
                break;
        }
    }
    return map;
}

const CaloRowData* CaloMap::find(std::size_t pindex) const {
    const auto it = m_rows.find(static_cast<int>(pindex));
    return it == m_rows.end() ? nullptr : &it->second;
}

// ---------------------------------------------------------------------------
// passes_pid_purity
// ---------------------------------------------------------------------------

bool passes_pid_purity(
    std::size_t cand_row,
    const CaloMap& calo,
    const ROOT::VecOps::RVec<float>& px,
    const ROOT::VecOps::RVec<float>& py,
    const ROOT::VecOps::RVec<float>& pz,
    const FeatureCuts& cuts) {
    check_momentum_columns(cand_row, px.size(), py.size(), pz.size());

    const CaloRowData* calo_poi = calo.find(cand_row);
    if (calo_poi == nullptr) return false;

    // BIT-IDENTICAL: the candidate's kinematics are computed in DOUBLE (the old
    // code cast each float component to double before squaring), unlike the
    // neighbours' -- see build_features.
    const double g_px = static_cast<double>(px[cand_row]);
    const double g_py = static_cast<double>(py[cand_row]);
    const double g_pz = static_cast<double>(pz[cand_row]);
    const double g_e = std::sqrt(g_px * g_px + g_py * g_py + g_pz * g_pz);
    const double g_theta_rad = std::acos(g_pz / g_e);

    return pass_pid_purity_filter(g_e, calo_poi->pcal.e, g_theta_rad, cuts);
}

// ---------------------------------------------------------------------------
// build_features
// ---------------------------------------------------------------------------

std::vector<float> build_features(
    std::size_t cand_row,
    const CaloMap& calo,
    const ROOT::VecOps::RVec<int>& pid,
    const ROOT::VecOps::RVec<float>& px,
    const ROOT::VecOps::RVec<float>& py,
    const ROOT::VecOps::RVec<float>& pz,
    const FeatureCuts& cuts) {
    check_momentum_columns(cand_row, px.size(), py.size(), pz.size());
    if (pid.size() != px.size()) {
        throw std::invalid_argument(
            "pi0::photonid::build_features: REC::Particle pid column has " + std::to_string(pid.size()) +
            " rows but the momentum columns have " + std::to_string(px.size()));
    }

    const CaloRowData* calo_poi_ptr = calo.find(cand_row);
    if (calo_poi_ptr == nullptr) {
        throw std::invalid_argument(
            "pi0::photonid::build_features: candidate row " + std::to_string(cand_row) +
            " has no REC::Calorimeter data; call passes_pid_purity() first, which reports this case "
            "as a rejection rather than an error");
    }
    const CaloRowData& calo_poi = *calo_poi_ptr;

    // --- candidate ("photon of interest") ------------------------------------
    //
    // BIT-IDENTICAL: double arithmetic throughout, because the old code declared
    // gPx/gPy/gPz as `double` via static_cast before squaring them. The
    // neighbour loop below does the same algebra in FLOAT. That inconsistency is
    // in the original and is preserved: the same photon gets a very slightly
    // different energy depending on whether it is the candidate or a neighbour.
    const double g_px = static_cast<double>(px[cand_row]);
    const double g_py = static_cast<double>(py[cand_row]);
    const double g_pz = static_cast<double>(pz[cand_row]);

    const double g_e = std::sqrt(g_px * g_px + g_py * g_py + g_pz * g_pz);
    const double g_theta_rad = std::acos(g_pz / g_e);  // RADIANS -- feature [2]

    const double g_epcal = calo_poi.pcal.e;
    const double g_m2u = calo_poi.pcal.m2u;
    const double g_m2v = calo_poi.pcal.m2v;

    // --- neighbour accumulators ---------------------------------------------
    double R_e = 0;
    double dE_e = 0;

    NeighbourSlots<kNumPhotonNeighbours> gam;
    NeighbourSlots<kNumChargedNeighbours> ch;
    NeighbourSlots<kNumNeutralNeighbours> nh;

    double num_photons_0_1 = 0;
    double num_photons_0_2 = 0;
    double num_photons_0_35 = 0;

    const Vec3 v_poi = calo_position(calo_poi);

    // --- neighbour loop, in REC::Particle row order --------------------------
    const std::size_t n_rows = px.size();
    for (std::size_t inner_row = 0; inner_row < n_rows; ++inner_row) {
        if (inner_row == cand_row) continue;

        const CaloRowData* calo_part_ptr = calo.find(inner_row);
        if (calo_part_ptr == nullptr) continue;
        const CaloRowData& calo_part = *calo_part_ptr;

        const int part_pid = pid[inner_row];
        const double mass = get_mass(part_pid);
        if (std::isnan(mass)) continue;  // unknown pid

        // BIT-IDENTICAL: FLOAT arithmetic. The old code took `auto` from
        // hipo::bank::getFloat, so `px*px + py*py + pz*pz` and its sqrt were
        // evaluated in float, and so was acos(pz/p). Only E promotes to double,
        // because `mass` is a double. Widening any of this to double would move
        // the dE_* features and the theta cut by ~1e-7 -- small, but not zero,
        // and the requirement is bit-identical.
        const float n_px = px[inner_row];
        const float n_py = py[inner_row];
        const float n_pz = pz[inner_row];
        const float p = std::sqrt(n_px * n_px + n_py * n_py + n_pz * n_pz);
        const double e_neighbour = std::sqrt(p * p + mass * mass);
        const float th_rad = std::acos(n_pz / p);

        // Forward-detector window. Applies to EVERY species, including the
        // electron branch below.
        const double th_deg = th_rad * 180.0 / M_PI;  // see pass_pid_purity_filter
        if (th_deg < cuts.theta_min_deg || th_deg > cuts.theta_max_deg) continue;

        const Vec3 v_part = calo_position(calo_part);
        const double R = angle_between_rad(v_poi, v_part);

        if (part_pid == 22) {
            // BIT-IDENTICAL: photon neighbours get the FULL purity filter
            // (energy AND PCAL energy AND theta), while charged/neutral hadron
            // neighbours get only the theta window applied above. Asymmetric in
            // the original; preserved.
            if (!pass_pid_purity_filter(e_neighbour, calo_part.pcal.e, th_rad, cuts)) continue;

            if (R < kPhotonCountR1) num_photons_0_1++;
            if (R < kPhotonCountR2) num_photons_0_2++;
            if (R < kPhotonCountR35) num_photons_0_35++;

            gam.insert(R, g_e - e_neighbour, calo_part.pcal.e, calo_part.pcal.m2u, calo_part.pcal.m2v);
        } else if (part_pid == 11) {
            // Nearest electron of ANY pid==11 row -- NOT a caller-designated
            // scattered electron. See the note on build_features in the header.
            if (R < R_e || R_e == 0) {
                R_e = R;
                dE_e = g_e - e_neighbour;
            }
        } else if (is_charged_hadron(part_pid)) {
            ch.insert(R, g_e - e_neighbour, calo_part.pcal.e, calo_part.pcal.m2u, calo_part.pcal.m2v);
        } else if (is_neutral_hadron(part_pid)) {
            nh.insert(R, g_e - e_neighbour, calo_part.pcal.e, calo_part.pcal.m2u, calo_part.pcal.m2v);
        }
        // Any other known pid (e.g. 111, -11) contributes nothing but is not
        // skipped earlier -- it just falls off the end of the chain, as before.
    }

    // --- emit, GROUPED BY VARIABLE ------------------------------------------
    std::vector<float> feats;
    feats.reserve(kNumFeatures);

    const auto push = [&feats](double v) { feats.push_back(static_cast<float>(v)); };
    const auto push_all = [&feats](const double* a, std::size_t n) {
        for (std::size_t i = 0; i < n; ++i) feats.push_back(static_cast<float>(a[i]));
    };

    push(g_e);          // [0]
    push(g_epcal);      // [1]
    push(g_theta_rad);  // [2] RADIANS
    push(g_m2u);        // [3]
    push(g_m2v);        // [4]
    push(R_e);          // [5]
    push(dE_e);         // [6]

    push_all(gam.R, kNumPhotonNeighbours);      // [ 7.. 9]
    push_all(gam.dE, kNumPhotonNeighbours);     // [10..12]
    push_all(gam.Epcal, kNumPhotonNeighbours);  // [13..15]
    push_all(gam.m2u, kNumPhotonNeighbours);    // [16..18]
    push_all(gam.m2v, kNumPhotonNeighbours);    // [19..21]

    push_all(ch.R, kNumChargedNeighbours);      // [22..23]
    push_all(ch.dE, kNumChargedNeighbours);     // [24..25]
    push_all(ch.Epcal, kNumChargedNeighbours);  // [26..27]
    push_all(ch.m2u, kNumChargedNeighbours);    // [28..29]
    push_all(ch.m2v, kNumChargedNeighbours);    // [30..31]

    push_all(nh.R, kNumNeutralNeighbours);      // [32..33]
    push_all(nh.dE, kNumNeutralNeighbours);     // [34..35]
    push_all(nh.Epcal, kNumNeutralNeighbours);  // [36..37]
    push_all(nh.m2u, kNumNeutralNeighbours);    // [38..39]
    push_all(nh.m2v, kNumNeutralNeighbours);    // [40..41]

    push(num_photons_0_1);   // [42]
    push(num_photons_0_2);   // [43]
    push(num_photons_0_35);  // [44]

    // Cheap insurance: the models index features[0..44] unconditionally, so a
    // short vector is a heap overread, not a wrong answer.
    if (feats.size() != kNumFeatures) {
        throw std::logic_error(
            "pi0::photonid::build_features: built " + std::to_string(feats.size()) +
            " features, expected " + std::to_string(kNumFeatures));
    }
    return feats;
}

}  // namespace pi0::photonid
