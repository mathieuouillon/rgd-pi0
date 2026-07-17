// Unit tests for src/vertex/VzCorrector -- the RG-D vertex-z correction and
// the per-target acceptance window.
//
// Three of these tests exist because of specific, named failure modes that the
// superseded analysis actually shipped or that this port could plausibly
// reintroduce:
//
//   1. THE NAMING TRAP. 'lo'/'hi' name the parameter set, not the z ordering:
//      mu_hi = -7.861 is the MORE NEGATIVE (upstream) foil, so Cu is the -7.9
//      peak and Sn is the -2.9 peak. The old VzCorrectorService.hpp doc-comment
//      asserted the opposite. Nothing about the *code* changes if you swap the
//      Cu/Sn routing -- it still compiles, still runs, still fills histograms --
//      it just silently analyses the wrong foil. So it is pinned here, from both
//      directions (each target accepts its peak AND rejects the other's).
//      The window itself now comes from config/cuts.json rather than from a
//      table in this library, so which peaks a target selects is DATA; these
//      tests pin the semantics, and tests/test_selection.cpp pins the shipped
//      file's numbers.
//
//   2. DEGREES. phi_bin_index takes degrees. The old analysis shipped a sector
//      helper documented as radians and fed degrees. The test below feeds
//      sector centres in degrees and demands the centre bin back, and
//      separately feeds the radian-valued equivalent and demands that it does
//      NOT land in the centre bin -- so a radian implementation cannot pass.
//
//   3. THE UNCHECKED CELL INDEX. The old loader read (sector, ip, iphi) off the
//      row and indexed the coefficient vectors with it WITHOUT a bounds check;
//      a typo'd sector in the file was an out-of-bounds write. That is now an
//      exception, and the test writes a malformed file to prove it.
//
// Reference values are computed independently OF the code under test -- either
// worked out in closed form here (the identity cell, the synthetic sigma-floor
// cell) or transcribed as literals from an offline evaluation of the raw params
// file. They are never produced by calling VzCorrector.
//
// NOTE (TAP): meson runs this with Catch2's TAP reporter, whose parser rejects
// the line WARN(...) emits. Use INFO(...) + CHECK(...), never WARN.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "vertex/VzCorrector.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using pi0::vertex::Polarity;
using pi0::vertex::Target;
using pi0::vertex::VertexCuts;
using pi0::vertex::VzCorrector;
using pi0::vertex::VzPeak;
using pi0::vertex::VzRule;
using pi0::vertex::VzTargetCuts;

namespace {

// Path to data/Vz/vz_corrector_params.txt, injected by meson (see
// tests/meson.build). No fallback on purpose: a test that silently skips
// because it could not find its data file is worse than one that fails.
#ifndef PI0_VZ_PARAMS_FILE
#error "PI0_VZ_PARAMS_FILE must be defined by the build (see tests/meson.build)"
#endif

const std::string& params_path() {
    static const std::string p = PI0_VZ_PARAMS_FILE;
    return p;
}

// Load once: the file is 460K x 4 variants and every test wants the same
// read-only view of it.
const VzCorrector& shipped() {
    static const VzCorrector vc = VzCorrector::load(params_path());
    return vc;
}

// ---------------------------------------------------------------------------
// The OUTBENDING vertex cut, TRANSCRIBED from config/cuts.json's /vertex block.
//
// Transcribed rather than loaded, because this binary links pi0_vertex ONLY --
// no pi0_config, no nlohmann, no ROOT (see tests/meson.build: keeping the link
// line bare is what keeps "src/vertex is plain C++17" checkable). So the split
// is:
//
//   * HERE: the acceptance SEMANTICS against known inputs -- that a target
//     accepts its own peaks, rejects the other's, honours n_sigma, and that LD2
//     is a flat raw range.
//   * tests/test_selection.cpp: that the SHIPPED cuts.json actually carries
//     these numbers, since that binary does load the real file.
//
// Neither half alone pins the Cu/Sn trap; together they do. If you change a
// number here you must change it there too, and cuts.json is what wins.
// ---------------------------------------------------------------------------
const VertexCuts& outbending_cuts() {
    static const VertexCuts v = [] {
        VertexCuts c;

        // LD2: raw v_z, correction disabled, flat window. Not a peak rule.
        c.ld2.correction_enabled = false;
        c.ld2.rule = VzRule::RawWindow;
        c.ld2.vz_min_cm = -15.0;
        c.ld2.vz_max_cm = 5.0;

        // Cu: corrected, the -7.9 (UPSTREAM) foil ONLY.
        c.cu.correction_enabled = true;
        c.cu.rule = VzRule::CorrectedPeaks;
        c.cu.n_sigma = 3.0;
        c.cu.peaks = {VzPeak{-7.861, 0.415}};

        // Sn: corrected, the -2.9 (DOWNSTREAM) foil ONLY.
        c.sn.correction_enabled = true;
        c.sn.rule = VzRule::CorrectedPeaks;
        c.sn.n_sigma = 3.0;
        c.sn.peaks = {VzPeak{-2.916, 0.370}};

        // CxC: corrected, EITHER foil -- both are carbon. Its peaks are close
        // to but NOT identical with the CuSn ones (-7.887/0.395 vs
        // -7.861/0.415): separate measurements, not to be deduplicated.
        c.cxc.correction_enabled = true;
        c.cxc.rule = VzRule::CorrectedPeaks;
        c.cxc.n_sigma = 3.0;
        c.cxc.peaks = {VzPeak{-7.887, 0.395}, VzPeak{-2.906, 0.373}};

        return c;
    }();
    return v;
}

// RAII temp file, for the malformed-input tests.
class TempFile {
   public:
    explicit TempFile(const std::string& contents) {
        static int counter = 0;
        m_path = (std::filesystem::temp_directory_path() /
                  ("pi0_vz_test_" + std::to_string(counter++) + ".txt"))
                     .string();
        std::ofstream out(m_path);
        out << contents;
    }
    ~TempFile() {
        std::error_code ec;
        std::filesystem::remove(m_path, ec);
    }
    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;

    [[nodiscard]] const std::string& path() const { return m_path; }

   private:
    std::string m_path;
};

}  // namespace

// ===========================================================================
// Loading
// ===========================================================================

TEST_CASE("params file loads all four variants", "[vertex][load]") {
    const VzCorrector& vc = shipped();

    const std::vector<std::string> expected = {
        "cusn_inbending", "cusn_outbending", "cxc_inbending", "cxc_outbending",
    };
    // variant_keys() sorts, so this compares as-is.
    CHECK(vc.variant_keys() == expected);

    // This project is outbending-only, but the inbending variants must stay
    // loadable -- that is the whole reason the file carries four.
    CHECK(vc.has_variant(Target::Cu, Polarity::Outbending));
    CHECK(vc.has_variant(Target::Cu, Polarity::Inbending));
    CHECK(vc.has_variant(Target::CxC, Polarity::Outbending));
    CHECK(vc.has_variant(Target::CxC, Polarity::Inbending));

    // Cu and Sn are two foils of ONE configuration and share the "cusn" fit.
    CHECK(*VzCorrector::variant_key(Target::Cu, Polarity::Outbending) == "cusn_outbending");
    CHECK(*VzCorrector::variant_key(Target::Sn, Polarity::Outbending) == "cusn_outbending");
    CHECK(*VzCorrector::variant_key(Target::CxC, Polarity::Outbending) == "cxc_outbending");
}

TEST_CASE("each variant yields 576 cells", "[vertex][load]") {
    const VzCorrector& vc = shipped();

    // 6 sectors x 16 p-bins x 6 phi-bins = 576.
    for (const auto target : {Target::Cu, Target::CxC}) {
        for (const auto pol : {Polarity::Outbending, Polarity::Inbending}) {
            const auto& p = vc.params_for(target, pol);
            INFO("variant " << *VzCorrector::variant_key(target, pol));

            CHECK(p.n_p == 16);
            CHECK(p.n_cells() == 576u);
            CHECK(6 * p.n_p * VzCorrector::kNPhi == 576);

            // Every coefficient vector must actually be that long -- n_cells()
            // is arithmetic on n_p and would report 576 even if the loader had
            // filled nothing.
            CHECK(p.mu_lo.size() == 576u);
            CHECK(p.mu_hi.size() == 576u);
            CHECK(p.sig_lo.size() == 576u);
            CHECK(p.sig_hi.size() == 576u);
            CHECK(p.theta_dom.size() == 576u);
        }
    }
}

TEST_CASE("p-binning spans 2 to 10 GeV", "[vertex][load]") {
    const auto& p = shipped().params_for(Target::Cu, Polarity::Outbending);
    CHECK_THAT(p.p_lo, WithinAbs(2.0, 1e-12));
    CHECK_THAT(p.p_binw, WithinAbs(0.5, 1e-12));
    // p_lo + n_p * p_binw = 2 + 16*0.5 = 10. This is the [2, 10] clamp domain.
    CHECK_THAT(p.p_lo + p.n_p * p.p_binw, WithinAbs(10.0, 1e-12));
}

TEST_CASE("cell_index matches the exporter's row order", "[vertex][load]") {
    const auto& p = shipped().params_for(Target::Cu, Polarity::Outbending);
    // (sector-1) * n_p * kNPhi + ip * kNPhi + iphi
    CHECK(p.cell_index(1, 0, 0) == 0u);
    CHECK(p.cell_index(1, 0, 5) == 5u);
    CHECK(p.cell_index(1, 1, 0) == 6u);
    CHECK(p.cell_index(2, 0, 0) == 96u);   // 1 * 16 * 6
    CHECK(p.cell_index(6, 15, 5) == 575u); // last cell
}

// ---------------------------------------------------------------------------
// Fail loudly
// ---------------------------------------------------------------------------

TEST_CASE("loader rejects a missing file", "[vertex][load][malformed]") {
    CHECK_THROWS_AS(VzCorrector::load("/nonexistent/vz_params.txt"), std::runtime_error);
}

TEST_CASE("loader rejects a schema mismatch", "[vertex][load][malformed]") {
    // n_phi = 5, but this build is compiled for 6.
    TempFile f("# bad schema\n1 5 6 4\n");
    CHECK_THROWS_AS(VzCorrector::load(f.path()), std::runtime_error);
}

TEST_CASE("loader rejects a truncated file", "[vertex][load][malformed]") {
    SECTION("header promises a variant that is not there") {
        TempFile f("1 6 6 4\n");
        CHECK_THROWS_AS(VzCorrector::load(f.path()), std::runtime_error);
    }
    SECTION("variant promises cells that are not there") {
        // n_p = 1 => 36 cells promised, one supplied.
        TempFile f(
            "1 6 6 4\n"
            "cusn_outbending -2.9 0.37 -7.9 0.41 2 0.5 1\n"
            "1 0 0  0 0 0 0  -10 0 0 0  1 0 0 0  1 0 0 0  nan nan\n");
        CHECK_THROWS_AS(VzCorrector::load(f.path()), std::runtime_error);
    }
}

TEST_CASE("loader rejects a non-positive n_p", "[vertex][load][malformed]") {
    TempFile f("1 6 6 4\ncusn_outbending -2.9 0.37 -7.9 0.41 2 0.5 0\n");
    CHECK_THROWS_AS(VzCorrector::load(f.path()), std::runtime_error);
}

TEST_CASE("loader rejects an out-of-range cell index", "[vertex][load][malformed]") {
    // THE MEMORY-SAFETY FIX. sector=9 does not exist. The old loader passed it
    // straight to cell_index() and wrote past the end of the coefficient
    // vectors; here it must throw.
    std::string s = "1 6 6 4\ncusn_outbending -2.9 0.37 -7.9 0.41 2 0.5 1\n";
    s += "9 0 0  0 0 0 0  -10 0 0 0  1 0 0 0  1 0 0 0  nan nan\n";
    for (int i = 1; i < 36; ++i) s += "1 0 1  0 0 0 0  -10 0 0 0  1 0 0 0  1 0 0 0  nan nan\n";
    TempFile f(s);
    CHECK_THROWS_AS(VzCorrector::load(f.path()), std::runtime_error);
}

TEST_CASE("loader rejects a duplicated cell", "[vertex][load][malformed]") {
    // 36 rows, all naming the same cell: the count is right but 35 cells would
    // silently keep their zero-initialised coefficients.
    std::string s = "1 6 6 4\ncusn_outbending -2.9 0.37 -7.9 0.41 2 0.5 1\n";
    for (int i = 0; i < 36; ++i) s += "1 0 0  0 0 0 0  -10 0 0 0  1 0 0 0  1 0 0 0  nan nan\n";
    TempFile f(s);
    CHECK_THROWS_AS(VzCorrector::load(f.path()), std::runtime_error);
}

TEST_CASE("loader rejects a short row", "[vertex][load][malformed]") {
    std::string s = "1 6 6 4\ncusn_outbending -2.9 0.37 -7.9 0.41 2 0.5 1\n";
    s += "1 0 0  0 0 0 0  -10 0 0 0  1 0 0 0  1 0 0\n";  // one coefficient short
    for (int i = 1; i < 36; ++i) s += "1 0 1  0 0 0 0  -10 0 0 0  1 0 0 0  1 0 0 0  nan nan\n";
    TempFile f(s);
    CHECK_THROWS_AS(VzCorrector::load(f.path()), std::runtime_error);
}

TEST_CASE("loader rejects a duplicate variant key", "[vertex][load][malformed]") {
    std::string block = "cusn_outbending -2.9 0.37 -7.9 0.41 2 0.5 1\n";
    for (int sec = 1; sec <= 6; ++sec)
        for (int iphi = 0; iphi < 6; ++iphi)
            block += std::to_string(sec) + " 0 " + std::to_string(iphi) +
                     "  0 0 0 0  -10 0 0 0  1 0 0 0  1 0 0 0  nan nan\n";
    TempFile f("2 6 6 4\n" + block + block);
    CHECK_THROWS_AS(VzCorrector::load(f.path()), std::runtime_error);
}

// ===========================================================================
// phi binning -- DEGREES
// ===========================================================================

TEST_CASE("phi_bin_index maps every sector centre to the centre bin", "[vertex][phi][degrees]") {
    // The sector centres, in DEGREES, straight from the geometry.
    const double centres_deg[6] = {0.0, 60.0, 120.0, 180.0, -120.0, -60.0};

    // A track at its sector's centre has phi_local = 0, which sits in bin
    // floor((0 + 30) / 10) = 3. Worked out here, not read off the code.
    for (int sector = 1; sector <= 6; ++sector) {
        INFO("sector " << sector << " centre " << centres_deg[sector - 1] << " deg");
        CHECK(VzCorrector::phi_bin_index(centres_deg[sector - 1], sector) == 3);
    }
}

TEST_CASE("phi_bin_index is fed DEGREES, not radians", "[vertex][phi][degrees]") {
    // Sector 2's centre is 60 deg == 1.0472 rad. Feeding the DEGREE value is
    // what lands in the centre bin.
    CHECK(VzCorrector::phi_bin_index(60.0, 2) == 3);

    // Feeding the RADIAN value of the same physical angle must NOT land in the
    // centre bin. If someone "fixes" this function to take radians, or a caller
    // starts passing radians, this fails. 1.0472 rad is read as 1.0472 deg,
    // which is 58.95 deg away from the sector 2 centre -- outside the +/-30 deg
    // sector entirely, so it clips to the -30 edge, i.e. bin 0.
    const double sixty_deg_in_rad = 60.0 * M_PI / 180.0;
    CHECK(VzCorrector::phi_bin_index(sixty_deg_in_rad, 2) == 0);
    CHECK(VzCorrector::phi_bin_index(sixty_deg_in_rad, 2) != 3);
}

TEST_CASE("phi_bin_index bins across the sector", "[vertex][phi]") {
    // Sector 1 has centre 0, so phi_local == phi_deg and the bin edges fall at
    // -30, -20, -10, 0, +10, +20, +30 degrees.
    CHECK(VzCorrector::phi_bin_index(-30.0, 1) == 0);
    CHECK(VzCorrector::phi_bin_index(-25.0, 1) == 0);
    CHECK(VzCorrector::phi_bin_index(-20.0, 1) == 1);
    CHECK(VzCorrector::phi_bin_index(-10.0, 1) == 2);
    CHECK(VzCorrector::phi_bin_index(0.0, 1) == 3);
    CHECK(VzCorrector::phi_bin_index(10.0, 1) == 4);
    CHECK(VzCorrector::phi_bin_index(20.0, 1) == 5);
    CHECK(VzCorrector::phi_bin_index(29.999, 1) == 5);

    // The upper edge is closed onto the last bin, not spilled into bin 6.
    CHECK(VzCorrector::phi_bin_index(30.0, 1) == 5);
}

TEST_CASE("phi_bin_index clips outside the sector and wraps at +/-180", "[vertex][phi]") {
    SECTION("outside the sector clips to the edge bins") {
        CHECK(VzCorrector::phi_bin_index(-40.0, 1) == 0);
        CHECK(VzCorrector::phi_bin_index(40.0, 1) == 5);
    }
    SECTION("phi is wrapped, so 370 deg == 10 deg") {
        CHECK(VzCorrector::phi_bin_index(370.0, 1) == VzCorrector::phi_bin_index(10.0, 1));
        CHECK(VzCorrector::phi_bin_index(-350.0, 1) == VzCorrector::phi_bin_index(10.0, 1));
    }
    SECTION("sector 4 straddles the +/-180 discontinuity") {
        // Its centre is 180; -180 is the same angle and must give the same bin.
        CHECK(VzCorrector::phi_bin_index(180.0, 4) == 3);
        CHECK(VzCorrector::phi_bin_index(-180.0, 4) == 3);
    }
    SECTION("sector 5 centre -120 == 240") {
        CHECK(VzCorrector::phi_bin_index(-120.0, 5) == 3);
        CHECK(VzCorrector::phi_bin_index(240.0, 5) == 3);
    }
    SECTION("out-of-range sector clamps") {
        CHECK(VzCorrector::phi_bin_index(0.0, 0) == VzCorrector::phi_bin_index(0.0, 1));
        CHECK(VzCorrector::phi_bin_index(0.0, 7) == VzCorrector::phi_bin_index(0.0, 6));
    }
}

// ===========================================================================
// correct()
// ===========================================================================

TEST_CASE("correct throws without an active variant", "[vertex][correct]") {
    VzCorrector vc;  // nothing loaded, nothing selected
    CHECK_THROWS_AS(vc.correct(-5.0, 4.0, 15.0, 0.0, 1), std::runtime_error);
}

TEST_CASE("correct is the identity on a degenerate cell", "[vertex][correct]") {
    // Cell (sector=1, ip=0, iphi=0) of cusn_outbending is a fallback cell: its
    // cubics are the constants (mu0_lo, mu0_hi, sigma0_lo, sigma0_hi) -- the
    // very reference values the correction standardises TO. So
    //     vz_lo = mu0_lo + sigma0_lo * (vz - mu0_lo)/sigma0_lo = vz
    // and likewise vz_hi = vz, hence the blend returns vz for ANY weight.
    // The correction must therefore be exactly the identity there, which
    // exercises the cell lookup, the Horner evaluation, the standardisation
    // and the blend against an answer known in closed form.
    VzCorrector vc = VzCorrector::load(params_path());
    vc.set_variant(Target::Cu, Polarity::Outbending);
    CHECK(vc.active_variant() == "cusn_outbending");

    // p = 2.0 -> bin 0; phi = -25 deg -> bin 0 in sector 1.
    for (const double vz : {-9.0, -7.861, -5.0, -2.916, 0.0, 3.5}) {
        INFO("vz = " << vz);
        CHECK_THAT(vc.correct(vz, 2.0, 15.0, -25.0, 1), WithinAbs(vz, 1e-12));
    }
}

TEST_CASE("correct reproduces an offline reference value", "[vertex][correct]") {
    // Cell (sector=1, ip=0, iphi=1) of cusn_outbending, theta_dom = [6.5, 25.5].
    // The literal below was evaluated offline, directly from the raw text of
    // data/Vz/vz_corrector_params.txt, by an independent implementation of the
    // published algorithm -- not by calling this class.
    VzCorrector vc = VzCorrector::load(params_path());
    vc.set_variant(Target::Cu, Polarity::Outbending);

    // p = 2.0 -> bin 0; phi = -15 deg -> bin 1 in sector 1; theta = 15 deg is
    // inside [6.5, 25.5] so it is not clipped.
    CHECK_THAT(vc.correct(-5.0, 2.0, 15.0, -15.0, 1),
               WithinRel(-3.8195917565454924, 1e-12));
}

TEST_CASE("correct clips theta to the cell's fitted domain", "[vertex][correct]") {
    VzCorrector vc = VzCorrector::load(params_path());
    vc.set_variant(Target::Cu, Polarity::Outbending);

    // Cell (1, 0, 1) has theta_dom = [6.5, 25.5] -- verified against the file.
    const auto& p = vc.params_for(Target::Cu, Polarity::Outbending);
    const auto idx = p.cell_index(1, 0, 1);
    REQUIRE_THAT(p.theta_dom[idx][0], WithinAbs(6.5, 1e-12));
    REQUIRE_THAT(p.theta_dom[idx][1], WithinAbs(25.5, 1e-12));

    // Below the domain -> evaluated at tmin. A cubic extrapolated to theta = 3
    // would give something else entirely.
    CHECK_THAT(vc.correct(-5.0, 2.0, 3.0, -15.0, 1),
               WithinRel(vc.correct(-5.0, 2.0, 6.5, -15.0, 1), 1e-12));
    // Above the domain -> evaluated at tmax.
    CHECK_THAT(vc.correct(-5.0, 2.0, 40.0, -15.0, 1),
               WithinRel(vc.correct(-5.0, 2.0, 25.5, -15.0, 1), 1e-12));

    // And the clipped answers are genuinely different from the interior one,
    // so the checks above are not vacuously comparing a constant to itself.
    CHECK(vc.correct(-5.0, 2.0, 3.0, -15.0, 1) != vc.correct(-5.0, 2.0, 15.0, -15.0, 1));
}

TEST_CASE("correct falls back to [6, 30] deg on a NaN theta domain", "[vertex][correct]") {
    // Cells with no statistics carry `nan nan` as their domain sentinel, and
    // the clip must fall back to the global fit range [6, 30] deg.
    //
    // WHY THIS TEST IS FUSSY ABOUT WHICH CELL IT USES
    // ----------------------------------------------
    // A missing fallback does NOT produce a NaN, which is the obvious thing to
    // check for and is useless here. `clamp(theta, NaN, NaN)` is written
    // `theta < lo ? lo : (hi < theta ? hi : theta)`, and EVERY comparison
    // against NaN is false -- so it silently returns theta, unclamped. Nothing
    // poisons, nothing throws; theta just quietly escapes the domain the cubic
    // was fitted over and the cubic gets extrapolated. So `isfinite` cannot see
    // this bug, and neither can any cell whose cubics are constant in theta
    // (611 of the 2304 cells have a NaN domain, and 606 of those are constant,
    // which makes theta irrelevant and the fallback unobservable there).
    //
    // Exactly 5 cells have BOTH a NaN domain and a theta-dependent cubic, and
    // all 5 are in the inbending variants -- verified by scanning the params
    // file. This uses one of them, which is also the only reason the inbending
    // variants earn their keep in this test file.
    VzCorrector vc = VzCorrector::load(params_path());
    vc.set_variant(Target::Cu, Polarity::Inbending);

    const auto& p = vc.params_for(Target::Cu, Polarity::Inbending);
    const auto idx = p.cell_index(2, 9, 5);
    REQUIRE(std::isnan(p.theta_dom[idx][0]));
    REQUIRE(std::isnan(p.theta_dom[idx][1]));

    // Cell (sector=2, ip=9, iphi=5): p-bin 9 is p in [6.5, 7.0), and phi = 85
    // deg is 25 deg from sector 2's centre (60 deg), i.e. phi-bin 5.
    const auto at = [&](double theta_deg) { return vc.correct(-5.0, 6.5, theta_deg, 85.0, 2); };

    // theta = 1 deg is below the fallback floor of 6 deg, so it must be clipped
    // to 6 and give exactly the theta = 6 answer. WITHOUT the fallback the
    // cubic is extrapolated to theta = 1 and returns -7.4018 instead: a 0.12 cm
    // shift, finite and plausible-looking, which is precisely why it needs a
    // test that compares values rather than one that checks for NaN.
    CHECK_THAT(at(1.0), WithinRel(at(6.0), 1e-12));

    // The literal is an offline evaluation of this cell straight from the raw
    // params text, with the fallback applied -- not a value this class produced.
    CHECK_THAT(at(1.0), WithinRel(-7.2860105783307754, 1e-12));

    // Pin the mistake too, so the test above cannot pass by both sides being
    // equally wrong: the un-clipped extrapolation is a genuinely different,
    // perfectly finite number.
    CHECK(std::fabs(-7.2860105783307754 - (-7.4018242459030379)) > 0.1);

    // And nothing anywhere in the theta range may go non-finite.
    for (const double theta_deg : {1.0, 6.0, 15.0, 30.0, 89.0}) {
        INFO("theta = " << theta_deg << " deg");
        CHECK(std::isfinite(at(theta_deg)));
    }
}

TEST_CASE("the NaN theta fallback clips at BOTH 6 and 30 deg", "[vertex][correct]") {
    // The test above pins the lower fallback (6 deg) against a real cell. The
    // UPPER one (30 deg) has no cell in the shipped file that can observe it:
    // of the 5 NaN-domain theta-dependent cells, none produces a different
    // answer above 30 deg. So it is pinned here on a synthetic variant instead
    // of left untested -- otherwise deleting `tmax_deg = kThetaHiFallbackDeg`
    // is a silent no-op on the whole suite.
    //
    // The cell is built so theta is maximally visible: mu_lo(theta) = theta,
    // everything else constant.
    std::string s = "1 6 6 4\n";
    // refs: mu0_lo = 0, sigma0_lo = 1, mu0_hi = -10, sigma0_hi = 1
    s += "cusn_outbending 0 1 -10 1 2 0.5 1\n";
    for (int sec = 1; sec <= 6; ++sec)
        for (int iphi = 0; iphi < 6; ++iphi)
            // mu_lo = theta | mu_hi = -10 | sig_lo = 1 | sig_hi = 1 | domain NaN
            s += std::to_string(sec) + " 0 " + std::to_string(iphi) +
                 "  0 1 0 0  -10 0 0 0  1 0 0 0  1 0 0 0  nan nan\n";
    TempFile f(s);

    VzCorrector vc = VzCorrector::load(f.path());
    vc.set_variant(Target::Cu, Polarity::Outbending);
    const auto at = [&](double theta_deg) { return vc.correct(30.0, 4.0, theta_deg, 0.0, 1); };

    // Worked out in closed form, e.g. at theta = 30 (the clip point):
    //   mu_lo = 30 -> z_lo = (30-30)/1 = 0    -> vz_lo = 0 + 1*0  = 0
    //   mu_hi = -10 -> z_hi = (30+10)/1 = 40  -> vz_hi = -10 + 40 = 30
    //   log_ratio = 0.5*(40^2 - 0^2) + log(1) = 800 -> w_lo = 1
    //   result = vz_hi + 1*(vz_lo - vz_hi) = 30 + (0 - 30) = 0
    CHECK_THAT(at(30.0), WithinAbs(0.0, 1e-9));
    CHECK_THAT(at(20.0), WithinAbs(10.0, 1e-9));  // interior, untouched
    CHECK_THAT(at(6.0), WithinAbs(24.0, 1e-9));

    // ABOVE the domain -> clipped to 30. Without the tmax fallback the cubic is
    // extrapolated to theta = 100 and the answer is 30, not 0.
    CHECK_THAT(at(100.0), WithinAbs(0.0, 1e-9));
    CHECK_THAT(at(100.0), WithinRel(at(30.0), 1e-12));

    // BELOW the domain -> clipped to 6. Without the tmin fallback: 29, not 24.
    CHECK_THAT(at(1.0), WithinAbs(24.0, 1e-9));
    CHECK_THAT(at(1.0), WithinRel(at(6.0), 1e-12));

    // The clip is doing real work: the interior answer differs from both edges.
    CHECK(at(20.0) != at(30.0));
    CHECK(at(20.0) != at(6.0));
}

TEST_CASE("correct clamps p to the [2, 10] GeV binning", "[vertex][correct]") {
    VzCorrector vc = VzCorrector::load(params_path());
    vc.set_variant(Target::Cu, Polarity::Outbending);

    // Below 2 GeV -> p-bin 0, same as p = 2.0 exactly.
    CHECK_THAT(vc.correct(-5.0, 0.5, 15.0, -15.0, 1),
               WithinRel(vc.correct(-5.0, 2.0, 15.0, -15.0, 1), 1e-12));
    CHECK_THAT(vc.correct(-5.0, 1.9, 15.0, -15.0, 1),
               WithinRel(vc.correct(-5.0, 2.0, 15.0, -15.0, 1), 1e-12));

    // At/above 10 GeV -> the last p-bin (15), same as p = 9.75.
    CHECK_THAT(vc.correct(-5.0, 10.0, 15.0, -15.0, 1),
               WithinRel(vc.correct(-5.0, 9.75, 15.0, -15.0, 1), 1e-12));
    CHECK_THAT(vc.correct(-5.0, 50.0, 15.0, -15.0, 1),
               WithinRel(vc.correct(-5.0, 9.75, 15.0, -15.0, 1), 1e-12));
}

TEST_CASE("correct clamps the sector", "[vertex][correct]") {
    VzCorrector vc = VzCorrector::load(params_path());
    vc.set_variant(Target::Cu, Polarity::Outbending);
    CHECK_THAT(vc.correct(-5.0, 4.0, 15.0, 0.0, 0), WithinRel(vc.correct(-5.0, 4.0, 15.0, 0.0, 1), 1e-12));
    CHECK_THAT(vc.correct(-5.0, 4.0, 15.0, 0.0, 9), WithinRel(vc.correct(-5.0, 4.0, 15.0, 0.0, 6), 1e-12));
}

TEST_CASE("correct floors sigma at 0.01", "[vertex][correct]") {
    // A synthetic single-p-bin variant whose sig_lo cubic is identically ZERO.
    // Without the floor, z_lo = (vz - mu_lo)/0 = -inf, the blend weight goes to
    // 0, and the result is 0 * -inf = NaN. With the floor the answer is finite
    // and, worked out in closed form below, exactly -2.
    std::string s = "1 6 6 4\n";
    // refs: mu0_lo = 0, sigma0_lo = 1, mu0_hi = -10, sigma0_hi = 1
    s += "cusn_outbending 0 1 -10 1 2 0.5 1\n";
    for (int sec = 1; sec <= 6; ++sec)
        for (int iphi = 0; iphi < 6; ++iphi)
            // mu_lo = 0 | mu_hi = -10 | sig_lo = 0 (!) | sig_hi = 1
            s += std::to_string(sec) + " 0 " + std::to_string(iphi) +
                 "  0 0 0 0  -10 0 0 0  0 0 0 0  1 0 0 0  nan nan\n";
    TempFile f(s);

    VzCorrector vc = VzCorrector::load(f.path());
    vc.set_variant(Target::Cu, Polarity::Outbending);

    // vz = -0.02, so with sig_lo floored to 0.01:
    //   z_lo  = (-0.02 - 0) / 0.01 = -2      -> vz_lo = 0 + 1*(-2)     = -2
    //   z_hi  = (-0.02 + 10) / 1   = 9.98    -> vz_hi = -10 + 1*(9.98) = -0.02
    //   log_ratio = 0.5*(9.98^2 - (-2)^2) + log(1/0.01) = 52.405...
    //   w_lo  = 1/(1 + exp(-52.4)) = 1 - O(1e-23)
    //   result = vz_hi + w_lo*(vz_lo - vz_hi) = -0.02 + 1*(-1.98) = -2
    const double got = vc.correct(-0.02, 4.0, 15.0, 0.0, 1);
    CHECK(std::isfinite(got));
    CHECK_THAT(got, WithinAbs(-2.0, 1e-12));
}

// ===========================================================================
// pass_window -- targets, peaks, and THE NAMING TRAP
// ===========================================================================

TEST_CASE("the cut's sigma OVERRIDES the params file's reference peak", "[vertex][window]") {
    // THE SPLIT THAT MUST SURVIVE: the CORRECTION POLYNOMIALS come from
    // data/Vz/vz_corrector_params.txt; the CUT WINDOW comes from
    // config/cuts.json. They disagree, on purpose. The file records
    // sigma0_hi = 0.3851 for cusn_outbending while the code that produced the
    // published result used 0.415 -- an 8% difference, open question Q2 in the
    // note -- and cuts.json carries 0.415.
    //
    // Pin the disagreement, so that "reconciling" the window to the params
    // file's own reference peaks is a visible change rather than a tidy-up.
    const auto& p = shipped().params_for(Target::Cu, Polarity::Outbending);
    CHECK_THAT(p.sigma0_hi, WithinAbs(0.3851194323, 1e-9));

    const VzPeak& cu_peak = outbending_cuts().cu.peaks.at(0);
    CHECK_THAT(cu_peak.sigma_cm, WithinAbs(0.415, 1e-12));
    CHECK(std::fabs(p.sigma0_hi - cu_peak.sigma_cm) > 0.02);
}

TEST_CASE("Cu's peak is the MORE NEGATIVE one", "[vertex][window][naming-trap]") {
    // ***** THE NAMING TRAP, PINNED (1 of 2) *****
    //
    // The fit's 'lo'/'hi' name its PARAMETER SET, not the z ordering: mu_hi =
    // -7.861 is the MORE NEGATIVE, further UPSTREAM foil. The old service's
    // doc-comment said the reverse. The cut no longer speaks lo/hi at all --
    // targets carry their own peaks -- so the surviving statement of the fact is
    // this one: the foil called Cu sits upstream of the foil called Sn.
    const VertexCuts& v = outbending_cuts();
    REQUIRE(v.cu.peaks.size() == 1u);
    REQUIRE(v.sn.peaks.size() == 1u);
    CHECK(v.cu.peaks.at(0).mu_cm < v.sn.peaks.at(0).mu_cm);

    // CxC occupies the same two z slots, so it has one peak near each foil.
    REQUIRE(v.cxc.peaks.size() == 2u);
    const double cxc_upstream = std::min(v.cxc.peaks.at(0).mu_cm, v.cxc.peaks.at(1).mu_cm);
    const double cxc_downstream = std::max(v.cxc.peaks.at(0).mu_cm, v.cxc.peaks.at(1).mu_cm);
    CHECK_THAT(cxc_upstream, WithinAbs(v.cu.peaks.at(0).mu_cm, 0.1));
    CHECK_THAT(cxc_downstream, WithinAbs(v.sn.peaks.at(0).mu_cm, 0.1));
    // ...but they are separate measurements, NOT the CuSn numbers reused.
    CHECK(cxc_upstream != v.cu.peaks.at(0).mu_cm);
    CHECK(cxc_downstream != v.sn.peaks.at(0).mu_cm);
}

TEST_CASE("pass_window accepts each peak centre and rejects 4 sigma away", "[vertex][window]") {
    const VertexCuts& v = outbending_cuts();

    SECTION("Sn: the -2.9 foil") {
        const VzPeak pk = v.sn.peaks.at(0);
        CHECK(VzCorrector::pass_window(pk.mu_cm, Target::Sn, v));
        // 4 sigma out, both sides. n_sigma is 3, so these must fail.
        CHECK_FALSE(VzCorrector::pass_window(pk.mu_cm + 4.0 * pk.sigma_cm, Target::Sn, v));
        CHECK_FALSE(VzCorrector::pass_window(pk.mu_cm - 4.0 * pk.sigma_cm, Target::Sn, v));
        // Just inside / just outside the 3 sigma edge.
        CHECK(VzCorrector::pass_window(pk.mu_cm + 2.99 * pk.sigma_cm, Target::Sn, v));
        CHECK_FALSE(VzCorrector::pass_window(pk.mu_cm + 3.01 * pk.sigma_cm, Target::Sn, v));
    }

    SECTION("Cu: the -7.9 foil") {
        const VzPeak pk = v.cu.peaks.at(0);
        CHECK(VzCorrector::pass_window(pk.mu_cm, Target::Cu, v));
        CHECK_FALSE(VzCorrector::pass_window(pk.mu_cm + 4.0 * pk.sigma_cm, Target::Cu, v));
        CHECK_FALSE(VzCorrector::pass_window(pk.mu_cm - 4.0 * pk.sigma_cm, Target::Cu, v));
        CHECK(VzCorrector::pass_window(pk.mu_cm + 2.99 * pk.sigma_cm, Target::Cu, v));
        CHECK_FALSE(VzCorrector::pass_window(pk.mu_cm + 3.01 * pk.sigma_cm, Target::Cu, v));
    }

    SECTION("CxC: either peak, since both foils are carbon") {
        for (const VzPeak& pk : v.cxc.peaks) {
            INFO("cxc peak at " << pk.mu_cm);
            CHECK(VzCorrector::pass_window(pk.mu_cm, Target::CxC, v));
            // 4 sigma from this peak, and nowhere near the other either.
            CHECK_FALSE(VzCorrector::pass_window(pk.mu_cm + 4.0 * pk.sigma_cm, Target::CxC, v));
            CHECK_FALSE(VzCorrector::pass_window(pk.mu_cm - 4.0 * pk.sigma_cm, Target::CxC, v));
        }
        // The gap between the peaks is empty.
        CHECK_FALSE(VzCorrector::pass_window(-5.4, Target::CxC, v));
    }
}

TEST_CASE("Cu and Sn select DIFFERENT peaks", "[vertex][window][naming-trap]") {
    // ***** THE NAMING TRAP, PINNED (2 of 2) *****
    //
    // Sn is the -2.9 (downstream) foil. Cu is the -7.9 (upstream) foil. The old
    // service's doc-comment claimed the reverse mapping to upstream/downstream;
    // config/cuts.json agrees with what is asserted here:
    //     /vertex/targets/cu/peaks = [{mu_cm: -7.861, ...}]
    //     /vertex/targets/sn/peaks = [{mu_cm: -2.916, ...}]
    //
    // The literals are deliberate: this test must fail if the fixture above is
    // "corrected" as well as if the routing is, and for Cu/Sn the vertex window
    // IS the target assignment -- swapping it silently analyses the wrong foil.
    //
    // Each target is checked from BOTH directions -- accepts its own peak AND
    // rejects the other's -- because a swapped routing still accepts *a* peak
    // and would otherwise look fine.
    const VertexCuts& v = outbending_cuts();
    const double sn_peak = -2.916;
    const double cu_peak = -7.861;

    SECTION("Sn accepts the -2.9 foil and rejects the -7.9 foil") {
        CHECK(VzCorrector::pass_window(sn_peak, Target::Sn, v));
        CHECK_FALSE(VzCorrector::pass_window(cu_peak, Target::Sn, v));
    }
    SECTION("Cu accepts the -7.9 foil and rejects the -2.9 foil") {
        CHECK(VzCorrector::pass_window(cu_peak, Target::Cu, v));
        CHECK_FALSE(VzCorrector::pass_window(sn_peak, Target::Cu, v));
    }
    SECTION("the two targets disagree about both peaks") {
        // The direct statement of "different peaks": at each foil centre,
        // exactly one of the two targets accepts.
        CHECK(VzCorrector::pass_window(sn_peak, Target::Sn, v) !=
              VzCorrector::pass_window(sn_peak, Target::Cu, v));
        CHECK(VzCorrector::pass_window(cu_peak, Target::Sn, v) !=
              VzCorrector::pass_window(cu_peak, Target::Cu, v));
    }
    SECTION("Cu is the upstream foil: its peak is the more negative one") {
        CHECK(cu_peak < sn_peak);
    }
}

TEST_CASE("LD2 is uncorrected and uses the raw (-15, 5) window", "[vertex][window][ld2]") {
    const VertexCuts& v = outbending_cuts();

    SECTION("LD2 has no correction variant at all") {
        CHECK_FALSE(VzCorrector::variant_key(Target::LD2, Polarity::Outbending).has_value());
        CHECK_FALSE(shipped().has_variant(Target::LD2, Polarity::Outbending));

        // Selecting it is an error rather than a silent no-op.
        VzCorrector vc = VzCorrector::load(params_path());
        CHECK_THROWS_AS(vc.set_variant(Target::LD2, Polarity::Outbending), std::runtime_error);

        // ...and the config says so too, rather than the skim inferring it from
        // the target's identity.
        CHECK(v.ld2.correction_enabled == false);
        CHECK(v.ld2.rule == VzRule::RawWindow);
        CHECK(v.cu.correction_enabled);
        CHECK(v.sn.correction_enabled);
        CHECK(v.cxc.correction_enabled);
    }

    SECTION("the window is the flat raw range, exclusive at both ends") {
        CHECK(VzCorrector::pass_window(0.0, Target::LD2, v));
        CHECK(VzCorrector::pass_window(-14.99, Target::LD2, v));
        CHECK(VzCorrector::pass_window(4.99, Target::LD2, v));

        CHECK_FALSE(VzCorrector::pass_window(-15.01, Target::LD2, v));
        CHECK_FALSE(VzCorrector::pass_window(5.01, Target::LD2, v));
        // Exclusive: the edges themselves are out.
        CHECK_FALSE(VzCorrector::pass_window(-15.0, Target::LD2, v));
        CHECK_FALSE(VzCorrector::pass_window(5.0, Target::LD2, v));
    }

    SECTION("the LD2 window ignores the foil peaks entirely") {
        // Both foil centres sit inside (-15, 5), so they pass -- LD2 is a long
        // cryotarget, not two foils. This is the observable difference between
        // "uncorrected, flat window" and "peak window".
        CHECK(VzCorrector::pass_window(-7.861, Target::LD2, v));
        CHECK(VzCorrector::pass_window(-2.916, Target::LD2, v));
        // And the gap between the foils, which Cu/Sn/CxC all reject, passes.
        CHECK(VzCorrector::pass_window(-5.4, Target::LD2, v));
        CHECK_FALSE(VzCorrector::pass_window(-5.4, Target::CxC, v));
    }

    // NOTE: the old "polarity does not change the LD2 window" section is gone
    // along with the polarity parameter itself. The window comes from
    // cuts.json's /vertex, which has no polarity dimension at all -- it is the
    // outbending set, and Cuts::load() refuses an inbending config outright
    // rather than apply it. tests/test_selection.cpp pins that refusal.
}

TEST_CASE("pass_window's two forms agree", "[vertex][window]") {
    // The target-routing form must be exactly for_target() fed to the
    // per-target form. Nothing may be decided in the routing.
    const VertexCuts& v = outbending_cuts();
    for (const auto t : {Target::LD2, Target::CxC, Target::Cu, Target::Sn}) {
        for (const double vz : {-9.0, -7.861, -5.4, -2.916, 0.0, 4.0, 6.0}) {
            INFO("target " << pi0::vertex::to_string(t) << ", vz = " << vz);
            CHECK(VzCorrector::pass_window(vz, t, v) ==
                  VzCorrector::pass_window(vz, v.for_target(t)));
        }
    }
}

TEST_CASE("for_target hands back the right block", "[vertex][window]") {
    // A mis-routed for_target() would silently cut every target against another
    // one's window -- and for Cu/Sn that is the target assignment itself.
    const VertexCuts& v = outbending_cuts();
    CHECK(&v.for_target(Target::LD2) == &v.ld2);
    CHECK(&v.for_target(Target::CxC) == &v.cxc);
    CHECK(&v.for_target(Target::Cu) == &v.cu);
    CHECK(&v.for_target(Target::Sn) == &v.sn);
}

TEST_CASE("a wider n_sigma admits what 3 sigma rejects", "[vertex][window]") {
    // Guards against n_sigma being ignored.
    const VzTargetCuts sn = outbending_cuts().sn;
    const VzPeak pk = sn.peaks.at(0);
    const double four_sigma_out = pk.mu_cm + 4.0 * pk.sigma_cm;

    CHECK_FALSE(VzCorrector::pass_window(four_sigma_out, sn));

    VzTargetCuts wider = sn;
    wider.n_sigma = 5.0;
    CHECK(VzCorrector::pass_window(four_sigma_out, wider));
}

TEST_CASE("an empty peak list accepts nothing", "[vertex][window]") {
    // The degenerate case, stated rather than left to inference: "inside ANY
    // listed peak" over an empty list is false, not true. Cuts::load() refuses
    // to build one of these (tests/test_selection.cpp), so this pins the
    // mechanism's own answer should one ever be constructed by hand.
    VzTargetCuts empty;
    empty.correction_enabled = true;
    empty.rule = VzRule::CorrectedPeaks;
    empty.n_sigma = 3.0;
    for (const double vz : {-9.0, -7.861, -5.4, -2.916, 0.0}) {
        CHECK_FALSE(VzCorrector::pass_window(vz, empty));
    }
}

TEST_CASE("corrected foil events land in their window", "[vertex][correct][window]") {
    // End to end: an event sitting exactly on a cell's fitted foil peak must
    // come out of correct() inside that foil's window. Uses the identity cell
    // (1, 0, 0), where correct() is the identity, so the reference peaks ARE
    // the fitted peaks and the two must agree by construction.
    //
    // This is the one test that crosses the params-file/cuts.json seam: mu0_hi
    // comes from the params file, the window from the (transcribed) config. It
    // passes because the two agree on the PEAK CENTRES -- they disagree only on
    // the hi peak's WIDTH, which is the subject of the override test above.
    const VertexCuts& v = outbending_cuts();
    VzCorrector vc = VzCorrector::load(params_path());
    vc.set_variant(Target::Cu, Polarity::Outbending);
    const auto& p = vc.params_for(Target::Cu, Polarity::Outbending);

    const double vz_cu = vc.correct(p.mu0_hi, 2.0, 15.0, -25.0, 1);
    const double vz_sn = vc.correct(p.mu0_lo, 2.0, 15.0, -25.0, 1);

    CHECK(VzCorrector::pass_window(vz_cu, Target::Cu, v));
    CHECK(VzCorrector::pass_window(vz_sn, Target::Sn, v));
    // ...and not in the other's.
    CHECK_FALSE(VzCorrector::pass_window(vz_cu, Target::Sn, v));
    CHECK_FALSE(VzCorrector::pass_window(vz_sn, Target::Cu, v));
}
