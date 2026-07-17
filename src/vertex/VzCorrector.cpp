#include "vertex/VzCorrector.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace pi0::vertex {

namespace {

/// Pull the next non-blank, non-comment line. False at EOF.
bool next_data_line(std::istream& in, std::string& out) {
    while (std::getline(in, out)) {
        const auto start = out.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        if (out[start] == '#') continue;
        out = out.substr(start);
        return true;
    }
    return false;
}

/// Parse the next whitespace-delimited token as a double via std::strtod.
///
/// NOT `iss >> d`: the params file uses the bare token "nan" as the
/// theta-domain sentinel for empty cells, and some libstdc++ versions'
/// num_get rejects it, which would turn every empty cell into a parse error.
/// strtod handles "nan" per C99. Kept from the old loader for that reason.
bool read_double(std::istream& in, double& out) {
    std::string tok;
    if (!(in >> tok)) return false;
    char* endp = nullptr;
    out = std::strtod(tok.c_str(), &endp);
    return endp != tok.c_str() && *endp == '\0';
}

/// Horner, lowest-order coefficient first.
double poly_eval(const std::array<double, VzCorrector::kNPoly>& c, double x) {
    double y = c[VzCorrector::kNPoly - 1];
    for (int k = VzCorrector::kNPoly - 2; k >= 0; --k) y = y * x + c[k];
    return y;
}

template <typename T>
T clamp_to(T v, T lo, T hi) {
    return v < lo ? lo : (hi < v ? hi : v);
}

}  // namespace

std::string to_string(Target t) {
    switch (t) {
        case Target::LD2: return "LD2";
        case Target::CxC: return "CxC";
        case Target::Cu:  return "Cu";
        case Target::Sn:  return "Sn";
    }
    return "?";
}

std::string to_string(Polarity p) {
    return p == Polarity::Outbending ? "outbending" : "inbending";
}

std::optional<std::string> VzCorrector::variant_key(Target target, Polarity polarity) {
    std::string t;
    switch (target) {
        // Cu and Sn are the two foils of the SAME CuSn configuration, so they
        // share one fit. The target only picks which peak window applies.
        case Target::Cu:
        case Target::Sn:  t = "cusn"; break;
        case Target::CxC: t = "cxc";  break;
        case Target::LD2: return std::nullopt;  // uncorrected by construction
    }
    return t + "_" + to_string(polarity);
}

// ---------------------------------------------------------------------------
// Loading
// ---------------------------------------------------------------------------

VzCorrector VzCorrector::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("VzCorrector: cannot open '" + path + "'");

    VzCorrector vc;
    std::string line;

    // ---- global header: n_variants n_phi n_sectors n_poly ----------------
    if (!next_data_line(in, line))
        throw std::runtime_error("VzCorrector: missing global header in '" + path + "'");

    int n_variants = 0, n_phi_chk = 0, n_sec_chk = 0, n_poly_chk = 0;
    {
        std::istringstream iss(line);
        iss >> n_variants >> n_phi_chk >> n_sec_chk >> n_poly_chk;
        if (!iss)
            throw std::runtime_error("VzCorrector: bad global header in '" + path + "'");
        if (n_phi_chk != kNPhi || n_sec_chk != kNSectors || n_poly_chk != kNPoly)
            throw std::runtime_error(
                "VzCorrector: schema mismatch in '" + path + "': file says n_phi=" +
                std::to_string(n_phi_chk) + " n_sectors=" + std::to_string(n_sec_chk) +
                " n_poly=" + std::to_string(n_poly_chk) + ", this build expects " +
                std::to_string(kNPhi) + "/" + std::to_string(kNSectors) + "/" +
                std::to_string(kNPoly));
        if (n_variants <= 0)
            throw std::runtime_error("VzCorrector: n_variants=" + std::to_string(n_variants) +
                                     " in '" + path + "'");
    }

    // ---- one block per variant -------------------------------------------
    for (int v = 0; v < n_variants; ++v) {
        if (!next_data_line(in, line))
            throw std::runtime_error("VzCorrector: truncated before variant " + std::to_string(v) +
                                     " of " + std::to_string(n_variants) + " in '" + path + "'");
        std::string key;
        Params p;
        {
            std::istringstream iss(line);
            iss >> key;
            const bool ok = read_double(iss, p.mu0_lo)
                         && read_double(iss, p.sigma0_lo)
                         && read_double(iss, p.mu0_hi)
                         && read_double(iss, p.sigma0_hi)
                         && read_double(iss, p.p_lo)
                         && read_double(iss, p.p_binw)
                         && static_cast<bool>(iss >> p.n_p);
            if (!ok)
                throw std::runtime_error("VzCorrector: bad variant header for '" + key + "' in '" +
                                         path + "'");
            if (p.n_p <= 0)
                throw std::runtime_error("VzCorrector: variant '" + key + "' has n_p=" +
                                         std::to_string(p.n_p));
            if (!(p.p_binw > 0.0))
                throw std::runtime_error("VzCorrector: variant '" + key + "' has p_binw=" +
                                         std::to_string(p.p_binw) + ", must be > 0");
            if (vc.m_variants.count(key))
                throw std::runtime_error("VzCorrector: duplicate variant '" + key + "' in '" +
                                         path + "'");
        }

        const std::size_t n_cells = p.n_cells();
        p.mu_lo.assign(n_cells, {});
        p.mu_hi.assign(n_cells, {});
        p.sig_lo.assign(n_cells, {});
        p.sig_hi.assign(n_cells, {});
        p.theta_dom.assign(n_cells, {0.0, 0.0});

        // The old loader trusted (sector, ip, iphi) straight off the row and
        // indexed with it -- no bounds check, so a malformed file wrote past
        // the end of the vectors. We validate, and we require each cell to be
        // written exactly once, so a duplicated or missing row is an error
        // instead of a silently stale cell.
        std::vector<bool> seen(n_cells, false);

        for (std::size_t i = 0; i < n_cells; ++i) {
            if (!next_data_line(in, line))
                throw std::runtime_error("VzCorrector: truncated cells in variant '" + key +
                                         "' (" + std::to_string(i) + " of " +
                                         std::to_string(n_cells) + " read)");
            std::istringstream iss(line);
            int sec = 0, ip = 0, iphi = 0;
            if (!(iss >> sec >> ip >> iphi))
                throw std::runtime_error("VzCorrector: bad cell key in variant '" + key +
                                         "', row " + std::to_string(i));
            if (sec < 1 || sec > kNSectors || ip < 0 || ip >= p.n_p || iphi < 0 || iphi >= kNPhi)
                throw std::runtime_error("VzCorrector: cell (sector=" + std::to_string(sec) +
                                         ", ip=" + std::to_string(ip) + ", iphi=" +
                                         std::to_string(iphi) + ") out of range in variant '" +
                                         key + "', row " + std::to_string(i));

            const std::size_t idx = p.cell_index(sec, ip, iphi);
            if (seen[idx])
                throw std::runtime_error("VzCorrector: duplicate cell (sector=" +
                                         std::to_string(sec) + ", ip=" + std::to_string(ip) +
                                         ", iphi=" + std::to_string(iphi) + ") in variant '" +
                                         key + "'");
            seen[idx] = true;

            bool ok = true;
            auto read4 = [&](std::array<double, kNPoly>& dst) {
                for (int k = 0; k < kNPoly; ++k) ok = ok && read_double(iss, dst[k]);
            };
            read4(p.mu_lo[idx]);
            read4(p.mu_hi[idx]);
            read4(p.sig_lo[idx]);
            read4(p.sig_hi[idx]);
            ok = ok && read_double(iss, p.theta_dom[idx][0])
                    && read_double(iss, p.theta_dom[idx][1]);
            if (!ok)
                throw std::runtime_error("VzCorrector: bad row in variant '" + key + "' at cell (" +
                                         std::to_string(sec) + ", " + std::to_string(ip) + ", " +
                                         std::to_string(iphi) + ")");
        }
        vc.m_variants.emplace(std::move(key), std::move(p));
    }

    return vc;
}

// ---------------------------------------------------------------------------
// Variant selection
// ---------------------------------------------------------------------------

void VzCorrector::set_variant(Target target, Polarity polarity) {
    const auto key = variant_key(target, polarity);
    if (!key)
        throw std::runtime_error("VzCorrector: target " + to_string(target) +
                                 " has no correction variant (LD2 is uncorrected by construction; "
                                 "use pass_window(), which routes it to the raw range)");
    const auto it = m_variants.find(*key);
    if (it == m_variants.end())
        throw std::runtime_error("VzCorrector: variant '" + *key + "' not loaded");
    m_active_key = *key;
    m_active = &it->second;
}

bool VzCorrector::has_variant(Target target, Polarity polarity) const {
    const auto key = variant_key(target, polarity);
    return key && m_variants.count(*key) > 0;
}

const VzCorrector::Params& VzCorrector::params_for(Target target, Polarity polarity) const {
    const auto key = variant_key(target, polarity);
    if (!key)
        throw std::runtime_error("VzCorrector: target " + to_string(target) + " has no variant");
    const auto it = m_variants.find(*key);
    if (it == m_variants.end())
        throw std::runtime_error("VzCorrector: variant '" + *key + "' not loaded");
    return it->second;
}

std::vector<std::string> VzCorrector::variant_keys() const {
    std::vector<std::string> keys;
    keys.reserve(m_variants.size());
    for (const auto& kv : m_variants) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    return keys;
}

const VzCorrector::Params& VzCorrector::active() const {
    if (!m_active)
        throw std::runtime_error("VzCorrector: no active variant. Call set_variant() first.");
    return *m_active;
}

// ---------------------------------------------------------------------------
// Binning
// ---------------------------------------------------------------------------

int VzCorrector::phi_bin_index(double phi_deg, int sector) {
    // DEGREES throughout -- see the header. The +540/-180 dance wraps
    // (phi_deg - centre) into [-180, 180) without assuming any input range,
    // reproducing Python's `%` (which returns a non-negative remainder) since
    // C's fmod keeps the dividend's sign.
    const int s = clamp_to(sector, 1, kNSectors);
    const double center_deg = kSectorPhiCenterDeg[static_cast<std::size_t>(s - 1)];
    double r = std::fmod(phi_deg - center_deg + 540.0, 360.0);
    if (r < 0.0) r += 360.0;
    double phi_local_deg = r - 180.0;

    // Clip into [-30, +30). The 1e-9 keeps phi_local == +30 exactly (the
    // upper edge) from binning to kNPhi instead of kNPhi-1.
    if (phi_local_deg < -kPhiHalfWidthDeg) {
        phi_local_deg = -kPhiHalfWidthDeg;
    } else if (phi_local_deg > kPhiHalfWidthDeg - 1e-9) {
        phi_local_deg = kPhiHalfWidthDeg - 1e-9;
    }
    const int iphi = static_cast<int>((phi_local_deg + kPhiHalfWidthDeg) / kPhiBinWidthDeg);
    return clamp_to(iphi, 0, kNPhi - 1);
}

// ---------------------------------------------------------------------------
// The correction
// ---------------------------------------------------------------------------

double VzCorrector::correct(double vz, double p, double theta_deg, double phi_deg, int sector) const {
    const Params& pr = active();

    const int s_clamped = clamp_to(sector, 1, kNSectors);
    // Clamping the BIN INDEX, not p -- same thing for the shipped binning
    // (p_lo=2, p_binw=0.5, n_p=16 => p clamped to [2, 10] GeV/c), and it is
    // what the old code did.
    const int ip = clamp_to(static_cast<int>((p - pr.p_lo) / pr.p_binw), 0, pr.n_p - 1);
    const int iphi = phi_bin_index(phi_deg, sector);
    const std::size_t idx = pr.cell_index(s_clamped, ip, iphi);

    // Empty cells carry a NaN domain; fall back to the global fit range.
    double tmin_deg = pr.theta_dom[idx][0];
    double tmax_deg = pr.theta_dom[idx][1];
    if (!std::isfinite(tmin_deg)) tmin_deg = kThetaLoFallbackDeg;
    if (!std::isfinite(tmax_deg)) tmax_deg = kThetaHiFallbackDeg;
    const double th_c_deg = clamp_to(theta_deg, tmin_deg, tmax_deg);

    const double mu_lo  = poly_eval(pr.mu_lo[idx], th_c_deg);
    const double mu_hi  = poly_eval(pr.mu_hi[idx], th_c_deg);
    const double sig_lo = std::max(poly_eval(pr.sig_lo[idx], th_c_deg), kSigmaFloor);
    const double sig_hi = std::max(poly_eval(pr.sig_hi[idx], th_c_deg), kSigmaFloor);

    // z-score against each peak, then standardise to that peak's reference.
    const double z_lo  = (vz - mu_lo) / sig_lo;
    const double z_hi  = (vz - mu_hi) / sig_hi;
    const double vz_lo = pr.mu0_lo + pr.sigma0_lo * z_lo;
    const double vz_hi = pr.mu0_hi + pr.sigma0_hi * z_hi;

    // Posterior P(lo | vz) under a two-Gaussian mixture with equal priors,
    // written as a sigmoid of the log-likelihood ratio.
    const double log_ratio = 0.5 * (z_hi * z_hi - z_lo * z_lo) + std::log(sig_hi / sig_lo);
    const double w_lo = 1.0 / (1.0 + std::exp(-log_ratio));
    return vz_hi + w_lo * (vz_lo - vz_hi);
}

// ---------------------------------------------------------------------------
// Windows
// ---------------------------------------------------------------------------

const VzTargetCuts& VertexCuts::for_target(Target t) const {
    switch (t) {
        case Target::LD2: return ld2;
        case Target::CxC: return cxc;
        case Target::Cu:  return cu;
        case Target::Sn:  return sn;
    }
    throw std::runtime_error("VertexCuts::for_target: unknown target");
}

bool VzCorrector::pass_window(double vz_cm, const VzTargetCuts& tc) {
    switch (tc.rule) {
        case VzRule::RawWindow:
            // Exclusive at both ends, matching the old service's
            // `vz > min && vz < max`. vz_cm is the RAW vz here.
            return vz_cm > tc.vz_min_cm && vz_cm < tc.vz_max_cm;

        case VzRule::CorrectedPeaks:
            // Inside ANY listed peak. WHICH peaks are listed is what makes this
            // one rule serve Cu (the -7.9 foil only), Sn (the -2.9 foil only)
            // and CxC (either, both foils being carbon) -- see VzTargetCuts.
            for (const VzPeak& pk : tc.peaks) {
                if (std::fabs(vz_cm - pk.mu_cm) < tc.n_sigma * pk.sigma_cm) return true;
            }
            return false;
    }
    return false;
}

bool VzCorrector::pass_window(double vz_cm, Target target, const VertexCuts& v) {
    return pass_window(vz_cm, v.for_target(target));
}

}  // namespace pi0::vertex
