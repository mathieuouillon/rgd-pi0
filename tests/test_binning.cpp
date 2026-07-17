// Unit tests for the frozen factorized binning (src/core/Binning).
//
// These tests exist to make four specific historical failures impossible. All
// four shipped in the analysis this code replaces:
//
//   1. A BINNING NOBODY CAN REPRODUCE. The kd-tree was built from an unseeded,
//      thread-timing-dependent reservoir sample: two passes over the same data
//      gave different edges. Nothing about that is testable, which is the
//      point -- the grids here are loaded from JSON, so their geometry is a
//      value a test can assert on, and provenance_hash() is the thing that
//      makes an output traceable back to it.
//   2. AN INDEX FORMULA WRITTEN DOWN NOWHERE. The old leaf formula had to be
//      REVERSE-ENGINEERED out of the output files. The round-trip test below is
//      the executable form of the formula documented in core/Binning.hpp: it
//      decodes every bin back to (i_q2, i_xb, i_z, i_pt2), so the header's
//      comment and the code cannot drift apart.
//   3. THE TOP-EDGE DROP. The old find_1d_bin() returned -1 for a value exactly
//      equal to the top edge, silently discarding every pT2 that sat on it.
//      Pinned below, by name.
//   4. SILENT DEFAULTS AND SILENT MIS-BINNING. Non-monotonic edges do not
//      crash: the binary search over them is merely meaningless, so events get
//      filed into the wrong bins and every downstream number is quietly wrong.
//      load() must refuse, and that refusal is tested.
//
// Method note, inherited from test_selection.cpp: reference values are stated
// as literals worked out from the SOURCE OF TRUTH (config/binning/*.json) and
// never by calling the function under test. Where an index is checked, the
// formula is re-derived locally from the per-axis bin counts rather than by
// asking Binning -- so a mistake in the formula has to be made twice,
// identically, to pass.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <nlohmann/json.hpp>

#include "core/Binning.hpp"

using Catch::Matchers::ContainsSubstring;

namespace {

// PI0_GRID_A_JSON / PI0_GRID_B_JSON are the absolute paths to the shipped
// config/binning/*.json, injected by tests/meson.build. The tests read the REAL
// files: a test that parses a hand-written copy of the config proves only that
// the copy parses.
const char* grid_a_path() { return PI0_GRID_A_JSON; }
const char* grid_b_path() { return PI0_GRID_B_JSON; }

const pi0::Binning& shipped_binning() {
    static const pi0::Binning b = pi0::Binning::load(grid_a_path(), grid_b_path());
    return b;
}

nlohmann::json parsed(const char* path) {
    std::ifstream in(path);
    REQUIRE(in.good());
    nlohmann::json doc;
    in >> doc;
    return doc;
}

/// Write `doc` to a scratch file and return the path. Used to build a grid that
/// is the real one MINUS or PLUS one defect, so the loader is tested against a
/// document realistic in every respect except the fault under test.
std::string write_temp_json(const nlohmann::json& doc, const std::string& tag) {
    const std::string path = std::string(PI0_TEST_TMPDIR) + "/grid_" + tag + ".json";
    std::ofstream out(path);
    REQUIRE(out.good());
    out << doc.dump(2);
    out.close();
    return path;
}

/// Load with a doctored Grid A and the shipped Grid B. Returns the temp path so
/// the caller can remove it.
struct TempGrid {
    std::string path;
    explicit TempGrid(const nlohmann::json& doc, const std::string& tag)
        : path(write_temp_json(doc, tag)) {}
    ~TempGrid() { std::remove(path.c_str()); }
    TempGrid(const TempGrid&) = delete;
    TempGrid& operator=(const TempGrid&) = delete;
};

/// A hand-built axis, for the semantics tests. Deliberately NOT the shipped
/// edges: find()'s contract is about half-open intervals and boundaries, and it
/// is clearest to read against small round numbers with a known answer.
pi0::Grid1D toy_axis() { return pi0::Grid1D{"toy", {0.0, 1.0, 2.0, 4.0}}; }  // 3 bins

}  // namespace

// ===========================================================================
// Grid1D::find -- half-open semantics
// ===========================================================================

TEST_CASE("Grid1D::nbins is one fewer than the edge count", "[binning]") {
    CHECK(toy_axis().nbins() == 3);
    CHECK(pi0::Grid1D{"two", {0.0, 1.0}}.nbins() == 1);

    // A default-constructed axis is inert rather than undefined: no edges, no
    // bins, and find() must not read edges.front() off an empty vector.
    const pi0::Grid1D empty{};
    CHECK(empty.nbins() == 0);
    CHECK(empty.find(0.0) == -1);
    CHECK(pi0::Grid1D{"one", {1.0}}.nbins() == 0);
    CHECK(pi0::Grid1D{"one", {1.0}}.find(1.0) == -1);
}

TEST_CASE("Grid1D::find bins are half-open [lo, hi)", "[binning]") {
    const pi0::Grid1D g = toy_axis();  // edges 0, 1, 2, 4

    // A value on a bin's LOWER edge belongs to that bin.
    CHECK(g.find(0.0) == 0);
    CHECK(g.find(1.0) == 1);
    CHECK(g.find(2.0) == 2);

    // A value strictly inside a bin belongs to it.
    CHECK(g.find(0.5) == 0);
    CHECK(g.find(1.5) == 1);
    CHECK(g.find(3.0) == 2);

    // A value on an interior bin's UPPER edge belongs to the NEXT bin, not to
    // it. This is the half-open rule stated from the other side; 1.0 is bin 0's
    // upper edge and bin 1's lower edge, and it must land in 1.
    CHECK(g.find(1.0) != 0);
    CHECK(g.find(2.0) != 1);

    // Just below an edge stays in the lower bin.
    CHECK(g.find(std::nextafter(1.0, 0.0)) == 0);
    CHECK(g.find(std::nextafter(2.0, 0.0)) == 1);
}

TEST_CASE("Grid1D::find puts a value ON THE TOP EDGE in the LAST bin", "[binning]") {
    // PINS THE OLD BUG. The superseded find_1d_bin() treated the top edge as
    // out of range and returned -1, so every pT2 sitting exactly on it was
    // silently DROPPED -- not mis-binned, dropped, with no counter anywhere
    // recording the loss. The top edge is the kinematic limit; discarding the
    // events that reach it is not a defensible reading of "outside the grid".
    const pi0::Grid1D g = toy_axis();  // edges 0, 1, 2, 4

    CHECK(g.find(4.0) == 2);  // == nbins() - 1, NOT -1
    CHECK(g.find(4.0) != -1);

    // The rule is the LAST bin specifically, not merely "some bin".
    CHECK(g.find(4.0) == g.nbins() - 1);

    // And it is the only place the half-open rule bends: one ulp above the top
    // edge is out, one ulp below is still the last bin.
    CHECK(g.find(std::nextafter(4.0, 5.0)) == -1);
    CHECK(g.find(std::nextafter(4.0, 0.0)) == 2);

    // Same rule on the shipped pT2 axis, which is the axis the bug actually
    // bit. Its top edge is 1.5 (config/binning/grid_B_z_pt2.json).
    const pi0::Grid1D& pt2 = shipped_binning().B.y;
    REQUIRE(pt2.name == "pt2");
    CHECK(pt2.find(1.5) == pt2.nbins() - 1);
    CHECK(pt2.find(1.5) == 4);
}

TEST_CASE("Grid1D::find returns -1 outside the axis", "[binning]") {
    const pi0::Grid1D g = toy_axis();  // edges 0, 1, 2, 4

    // Below range.
    CHECK(g.find(-0.001) == -1);
    CHECK(g.find(-1e9) == -1);
    CHECK(g.find(std::nextafter(0.0, -1.0)) == -1);

    // Above range.
    CHECK(g.find(4.001) == -1);
    CHECK(g.find(1e9) == -1);

    // There are NO flow bins. -1 is the only out-of-range answer; it is never a
    // valid bin index, so a caller cannot mistake one for the other.
    CHECK(g.find(-1e9) == g.find(1e9));

    // NaN has no ordering against the edges, so it must take the -1 branch
    // rather than reach the binary search and land in an arbitrary bin.
    CHECK(g.find(std::numeric_limits<double>::quiet_NaN()) == -1);
    CHECK(g.find(std::numeric_limits<double>::infinity()) == -1);
    CHECK(g.find(-std::numeric_limits<double>::infinity()) == -1);
}

// ===========================================================================
// Grid2D
// ===========================================================================

TEST_CASE("Grid2D::find is row-major with y fast", "[binning]") {
    // x: 2 bins (edges 0,1,2), y: 3 bins (edges 0,1,2,4).
    const pi0::Grid2D g{pi0::Grid1D{"x", {0.0, 1.0, 2.0}}, toy_axis()};
    REQUIRE(g.x.nbins() == 2);
    REQUIRE(g.y.nbins() == 3);
    CHECK(g.ncells() == 6);

    // Formula re-derived locally, NOT by asking Grid2D: cell = ix * n_y + iy.
    CHECK(g.find(0.5, 0.5) == 0 * 3 + 0);
    CHECK(g.find(0.5, 1.5) == 0 * 3 + 1);
    CHECK(g.find(0.5, 3.0) == 0 * 3 + 2);
    CHECK(g.find(1.5, 0.5) == 1 * 3 + 0);
    CHECK(g.find(1.5, 3.0) == 1 * 3 + 2);

    // Every cell index is in [0, ncells()).
    CHECK(g.find(1.5, 3.0) == g.ncells() - 1);
}

TEST_CASE("Grid2D::find returns -1 if EITHER axis misses", "[binning]") {
    const pi0::Grid2D g{pi0::Grid1D{"x", {0.0, 1.0, 2.0}}, toy_axis()};

    CHECK(g.find(-1.0, 0.5) == -1);  // x below
    CHECK(g.find(99.0, 0.5) == -1);  // x above
    CHECK(g.find(0.5, -1.0) == -1);  // y below
    CHECK(g.find(0.5, 99.0) == -1);  // y above
    CHECK(g.find(-1.0, 99.0) == -1); // both

    // There is no partial hit: a good x with a bad y is not "x's cell".
    REQUIRE(g.find(0.5, 0.5) == 0);
    CHECK(g.find(0.5, 99.0) != 0);
}

// ===========================================================================
// The shipped grids
// ===========================================================================

TEST_CASE("the shipped grids have the geometry the note specifies", "[binning]") {
    const pi0::Binning& b = shipped_binning();

    // Axis names and order, transcribed from config/binning/*.json. Order is
    // load-bearing: it is the slow-then-fast order of the flat index formula.
    CHECK(b.A.x.name == "q2");
    CHECK(b.A.y.name == "xb");
    CHECK(b.B.x.name == "z");
    CHECK(b.B.y.name == "pt2");

    // Grid A = 8 x 7 = 56, Grid B = 5 x 5 = 25, 4D = 1400, 3D = 280.
    CHECK(b.A.x.nbins() == 8);
    CHECK(b.A.y.nbins() == 7);
    CHECK(b.A.ncells() == 56);
    CHECK(b.B.x.nbins() == 5);
    CHECK(b.B.y.nbins() == 5);
    CHECK(b.B.ncells() == 25);
    CHECK(b.n4d() == 1400);
    CHECK(b.n3d() == 280);

    // Edges, transcribed from the JSON (the source of truth). If these ever
    // disagree with it, IT WINS and this test is what changes -- but a
    // disagreement means someone moved a bin edge, and that must be a
    // deliberate act, which is exactly why these assertions force it.
    CHECK(b.A.x.edges == std::vector<double>{1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 11.0});
    CHECK(b.A.y.edges == std::vector<double>{0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7});
    CHECK(b.B.x.edges == std::vector<double>{0.0, 0.2, 0.4, 0.6, 0.8, 1.0});
    CHECK(b.B.y.edges == std::vector<double>{0.0, 0.1, 0.2, 0.4, 0.7, 1.5});
}

TEST_CASE("the shipped grids are still PLACEHOLDERS", "[binning]") {
    // This test is a TRIPWIRE, not a specification. It fails the day someone
    // replaces the placeholders with make_grid output -- at which point delete
    // it, and delete the placeholder warnings in the JSON with it. Until then
    // it stops the placeholder status from being forgotten, which is exactly
    // how the superseded analysis ended up quoting numbers nobody could source.
    for (const char* path : {grid_a_path(), grid_b_path()}) {
        const nlohmann::json doc = parsed(path);
        INFO("grid file: " << path);
        CHECK(doc.at("provenance").at("source").get<std::string>() == "placeholder");
        CHECK(doc.at("provenance").at("n_events").get<long long>() == 0);
    }
}

// ===========================================================================
// find_4d / find_3d -- composition and the flat index
// ===========================================================================

TEST_CASE("Binning::find_4d composes A and B with the B cell fast", "[binning]") {
    const pi0::Binning& b = shipped_binning();

    // Q2 = 1.2 -> i_q2 = 0; xB = 0.12 -> i_xb = 0; z = 0.1 -> i_z = 0;
    // pT2 = 0.05 -> i_pt2 = 0. The very first bin.
    CHECK(b.find_4d(1.2, 0.12, 0.1, 0.05) == 0);

    // Formula re-derived locally from the per-axis indices, exactly as the
    // downstream Python must:
    //     bin4d = ((i_q2 * n_xb) + i_xb) * (n_z * n_pt2) + (i_z * n_pt2) + i_pt2
    const int n_xb = b.A.y.nbins();
    const int n_z = b.B.x.nbins();
    const int n_pt2 = b.B.y.nbins();

    // Per-axis indices, read off the edge arrays in config/binning/*.json:
    //   Q2  = 3.5  -> i_q2  = 4, box [3.0, 4.0)
    //   xB  = 0.27 -> i_xb  = 3, box [0.25, 0.3)   <- edges are 0.1, 0.15, 0.2,
    //                                                 0.25, ... so [0.25,0.3)
    //                                                 is the FOURTH bin, index 3
    //   z   = 0.45 -> i_z   = 2, box [0.4, 0.6)
    //   pT2 = 0.3  -> i_pt2 = 2, box [0.2, 0.4)
    const int expect = ((4 * n_xb) + 3) * (n_z * n_pt2) + (2 * n_pt2) + 2;
    CHECK(expect == 787);  // spelled out, so the formula and the arithmetic are both pinned
    CHECK(b.find_4d(3.5, 0.27, 0.45, 0.3) == expect);

    // The composition is exactly A's cell times B's cell count plus B's cell --
    // asserted against the grids themselves, not just against the literal.
    CHECK(b.find_4d(3.5, 0.27, 0.45, 0.3) ==
          b.A.find(3.5, 0.27) * b.B.ncells() + b.B.find(0.45, 0.3));

    // The last bin: everything on its top edge.
    CHECK(b.find_4d(11.0, 0.7, 1.0, 1.5) == b.n4d() - 1);

    // Every in-range bin index is in [0, n4d()).
    CHECK(b.find_4d(3.5, 0.27, 0.45, 0.3) >= 0);
    CHECK(b.find_4d(3.5, 0.27, 0.45, 0.3) < b.n4d());
}

TEST_CASE("Binning::find_4d returns -1 if EITHER grid misses", "[binning]") {
    const pi0::Binning& b = shipped_binning();

    // A good 4-tuple, so that each case below spoils exactly one coordinate.
    REQUIRE(b.find_4d(3.5, 0.27, 0.45, 0.3) >= 0);

    CHECK(b.find_4d(0.5, 0.27, 0.45, 0.3) == -1);   // Q2 below Grid A
    CHECK(b.find_4d(50.0, 0.27, 0.45, 0.3) == -1);  // Q2 above Grid A
    CHECK(b.find_4d(3.5, 0.05, 0.45, 0.3) == -1);   // xB below Grid A
    CHECK(b.find_4d(3.5, 0.95, 0.45, 0.3) == -1);   // xB above Grid A
    CHECK(b.find_4d(3.5, 0.27, -0.1, 0.3) == -1);   // z below Grid B
    CHECK(b.find_4d(3.5, 0.27, 1.5, 0.3) == -1);    // z above Grid B
    CHECK(b.find_4d(3.5, 0.27, 0.45, -0.1) == -1);  // pT2 below Grid B
    CHECK(b.find_4d(3.5, 0.27, 0.45, 9.9) == -1);   // pT2 above Grid B

    // A miss in A alone is enough even when B hits, and vice versa: there is no
    // partial 4D bin.
    REQUIRE(b.A.find(0.5, 0.27) == -1);
    REQUIRE(b.B.find(0.45, 0.3) >= 0);
    CHECK(b.find_4d(0.5, 0.27, 0.45, 0.3) == -1);

    REQUIRE(b.A.find(3.5, 0.27) >= 0);
    REQUIRE(b.B.find(0.45, 9.9) == -1);
    CHECK(b.find_4d(3.5, 0.27, 0.45, 9.9) == -1);

    CHECK(b.find_4d(std::numeric_limits<double>::quiet_NaN(), 0.27, 0.45, 0.3) == -1);
}

TEST_CASE("the 4D flat index round-trips over EVERY bin", "[binning]") {
    // THE EXECUTABLE FORM OF THE HEADER'S FORMULA. The superseded analysis's
    // leaf index was documented nowhere and had to be reverse-engineered from
    // output files by chaining box centres. This test decodes every one of the
    // 1400 bins back to (i_q2, i_xb, i_z, i_pt2) using the decoder written in
    // core/Binning.hpp -- so if the code and that comment ever disagree, this
    // fails, and the downstream Python that implements the same decoder is
    // safe.
    const pi0::Binning& b = shipped_binning();

    const int n_q2 = b.A.x.nbins();
    const int n_xb = b.A.y.nbins();
    const int n_z = b.B.x.nbins();
    const int n_pt2 = b.B.y.nbins();

    std::vector<bool> seen(static_cast<std::size_t>(b.n4d()), false);

    for (int i_q2 = 0; i_q2 < n_q2; ++i_q2) {
        for (int i_xb = 0; i_xb < n_xb; ++i_xb) {
            for (int i_z = 0; i_z < n_z; ++i_z) {
                for (int i_pt2 = 0; i_pt2 < n_pt2; ++i_pt2) {
                    // A point strictly inside the (i_q2, i_xb, i_z, i_pt2) box:
                    // its midpoint. Used ONLY to name a bin -- this is not an
                    // abscissa, and the box midpoint is never what a result is
                    // reported at (see core/Binning.hpp's preamble).
                    const double q2 = 0.5 * (b.A.x.edges[i_q2] + b.A.x.edges[i_q2 + 1]);
                    const double xb = 0.5 * (b.A.y.edges[i_xb] + b.A.y.edges[i_xb + 1]);
                    const double z = 0.5 * (b.B.x.edges[i_z] + b.B.x.edges[i_z + 1]);
                    const double pt2 = 0.5 * (b.B.y.edges[i_pt2] + b.B.y.edges[i_pt2 + 1]);

                    const int bin = b.find_4d(q2, xb, z, pt2);
                    INFO("i_q2=" << i_q2 << " i_xb=" << i_xb << " i_z=" << i_z
                                 << " i_pt2=" << i_pt2 << " -> bin " << bin);
                    REQUIRE(bin >= 0);
                    REQUIRE(bin < b.n4d());

                    // DECODE, per the header:
                    //   b_cell = bin % (n_z * n_pt2)   a_cell = bin / (n_z * n_pt2)
                    //   i_pt2  = b_cell % n_pt2        i_z    = b_cell / n_pt2
                    //   i_xb   = a_cell % n_xb         i_q2   = a_cell / n_xb
                    const int b_cell = bin % (n_z * n_pt2);
                    const int a_cell = bin / (n_z * n_pt2);
                    CHECK(b_cell % n_pt2 == i_pt2);
                    CHECK(b_cell / n_pt2 == i_z);
                    CHECK(a_cell % n_xb == i_xb);
                    CHECK(a_cell / n_xb == i_q2);

                    // The decode agrees with what the grids themselves say.
                    CHECK(a_cell == b.A.find(q2, xb));
                    CHECK(b_cell == b.B.find(z, pt2));

                    // No two boxes share an index.
                    REQUIRE_FALSE(seen[static_cast<std::size_t>(bin)]);
                    seen[static_cast<std::size_t>(bin)] = true;
                }
            }
        }
    }

    // The map is onto: all 1400 indices are reachable. A product grid has no
    // kinematically-empty holes in its INDEX space (it may well have empty
    // bins in the data -- that is a different statement).
    for (std::size_t i = 0; i < seen.size(); ++i) {
        INFO("bin " << i << " is unreachable");
        CHECK(seen[i]);
    }
}

TEST_CASE("Binning::find_3d is the A cell times the z bin", "[binning]") {
    const pi0::Binning& b = shipped_binning();
    const int n_xb = b.A.y.nbins();
    const int n_z = b.B.x.nbins();

    CHECK(b.find_3d(1.2, 0.12, 0.1) == 0);

    // Formula re-derived locally: bin3d = ((i_q2 * n_xb) + i_xb) * n_z + i_z.
    // Q2 = 3.5 -> i_q2 = 4; xB = 0.27 -> i_xb = 3 (box [0.25, 0.3));
    // z = 0.45 -> i_z = 2.
    CHECK(((4 * n_xb) + 3) * n_z + 2 == 157);  // pin the arithmetic too
    CHECK(b.find_3d(3.5, 0.27, 0.45) == ((4 * n_xb) + 3) * n_z + 2);
    CHECK(b.find_3d(11.0, 0.7, 1.0) == b.n3d() - 1);

    // -1 if the A cell or the z bin misses. pT2 is deliberately not an
    // argument: it is the OBSERVABLE of the pT-broadening measurement, not a
    // binning dimension, so a pT2 outside Grid B must not remove an event from
    // its 3D bin.
    CHECK(b.find_3d(0.5, 0.27, 0.45) == -1);
    CHECK(b.find_3d(3.5, 0.95, 0.45) == -1);
    CHECK(b.find_3d(3.5, 0.27, 1.5) == -1);

    // The 3D round-trip decodes with the header's decoder.
    for (int i_q2 = 0; i_q2 < b.A.x.nbins(); ++i_q2) {
        for (int i_xb = 0; i_xb < n_xb; ++i_xb) {
            for (int i_z = 0; i_z < n_z; ++i_z) {
                const double q2 = 0.5 * (b.A.x.edges[i_q2] + b.A.x.edges[i_q2 + 1]);
                const double xb = 0.5 * (b.A.y.edges[i_xb] + b.A.y.edges[i_xb + 1]);
                const double z = 0.5 * (b.B.x.edges[i_z] + b.B.x.edges[i_z + 1]);

                const int bin = b.find_3d(q2, xb, z);
                INFO("i_q2=" << i_q2 << " i_xb=" << i_xb << " i_z=" << i_z);
                REQUIRE(bin >= 0);
                REQUIRE(bin < b.n3d());
                const int a_cell = bin / n_z;
                CHECK(bin % n_z == i_z);
                CHECK(a_cell % n_xb == i_xb);
                CHECK(a_cell / n_xb == i_q2);
            }
        }
    }
}

TEST_CASE("find_3d is not derived from find_4d", "[binning]") {
    const pi0::Binning& b = shipped_binning();
    const int n_pt2 = b.B.y.nbins();

    // They agree whenever both are in range...
    CHECK(b.find_3d(3.5, 0.27, 0.45) == b.find_4d(3.5, 0.27, 0.45, 0.3) / n_pt2);

    // ...but the agreement is a coincidence of the index layout, not a
    // contract, and it breaks exactly where it matters: a pT2 outside Grid B
    // kills the 4D bin and must NOT kill the 3D one, because pT2 is the
    // pT-broadening measurement's observable. A caller who computes bin3d as
    // bin4d / n_pt2 silently drops the high-pT2 tail -- which is the tail the
    // measurement is about.
    REQUIRE(b.find_4d(3.5, 0.27, 0.45, 9.9) == -1);
    CHECK(b.find_3d(3.5, 0.27, 0.45) >= 0);
}

// ===========================================================================
// load() -- fails loudly
// ===========================================================================

TEST_CASE("Binning::load throws on a missing file", "[binning]") {
    CHECK_THROWS_AS(pi0::Binning::load("/nonexistent/grid_A.json", grid_b_path()),
                    std::runtime_error);
    CHECK_THROWS_AS(pi0::Binning::load(grid_a_path(), "/nonexistent/grid_B.json"),
                    std::runtime_error);

    // The message must name the file, or it is not actionable.
    CHECK_THROWS_WITH(pi0::Binning::load("/nonexistent/grid_A.json", grid_b_path()),
                      ContainsSubstring("/nonexistent/grid_A.json"));
}

TEST_CASE("Binning::load throws on a missing key", "[binning]") {
    SECTION("axes") {
        nlohmann::json doc = parsed(grid_a_path());
        REQUIRE(doc.contains("axes"));  // not vacuous: fails if the key is renamed
        doc.erase("axes");
        const TempGrid t(doc, "no_axes");
        CHECK_THROWS_WITH(pi0::Binning::load(t.path, grid_b_path()),
                          ContainsSubstring("axes") && ContainsSubstring("MISSING"));
    }

    SECTION("an axis's edges") {
        nlohmann::json doc = parsed(grid_a_path());
        REQUIRE(doc["axes"][0].contains("edges"));
        doc["axes"][0].erase("edges");
        const TempGrid t(doc, "no_edges");
        CHECK_THROWS_WITH(pi0::Binning::load(t.path, grid_b_path()),
                          ContainsSubstring("axes[0].edges") && ContainsSubstring("MISSING"));
    }

    SECTION("an axis's name") {
        nlohmann::json doc = parsed(grid_b_path());
        REQUIRE(doc["axes"][1].contains("name"));
        doc["axes"][1].erase("name");
        const TempGrid t(doc, "no_name");
        CHECK_THROWS_WITH(pi0::Binning::load(grid_a_path(), t.path),
                          ContainsSubstring("axes[1].name") && ContainsSubstring("MISSING"));
    }

    SECTION("provenance") {
        // /provenance is required, so that it cannot decay into another
        // declared-set-and-never-read key -- the failure mode config/README.md
        // is a monument to.
        nlohmann::json doc = parsed(grid_a_path());
        REQUIRE(doc.contains("provenance"));
        doc.erase("provenance");
        const TempGrid t(doc, "no_prov");
        CHECK_THROWS_WITH(pi0::Binning::load(t.path, grid_b_path()),
                          ContainsSubstring("provenance") && ContainsSubstring("MISSING"));
    }

    SECTION("provenance.source") {
        nlohmann::json doc = parsed(grid_a_path());
        REQUIRE(doc["provenance"].contains("source"));
        doc["provenance"].erase("source");
        const TempGrid t(doc, "no_prov_source");
        CHECK_THROWS_WITH(pi0::Binning::load(t.path, grid_b_path()),
                          ContainsSubstring("provenance.source"));
    }

    SECTION("provenance.n_events") {
        nlohmann::json doc = parsed(grid_b_path());
        REQUIRE(doc["provenance"].contains("n_events"));
        doc["provenance"].erase("n_events");
        const TempGrid t(doc, "no_prov_nevents");
        CHECK_THROWS_WITH(pi0::Binning::load(grid_a_path(), t.path),
                          ContainsSubstring("provenance.n_events"));
    }
}

TEST_CASE("Binning::load throws on NON-MONOTONIC edges", "[binning]") {
    // THE FAILURE THAT DOES NOT ANNOUNCE ITSELF. Out-of-order edges do not
    // crash: upper_bound over an unsorted range is merely meaningless. The
    // events get filed into the wrong bins, every histogram fills, every fit
    // converges, and every number downstream is quietly wrong. If this is not
    // caught at load it is not caught at all.
    SECTION("a descending pair") {
        nlohmann::json doc = parsed(grid_a_path());
        // Q2 edges 1.0, 1.5, 2.0, ... -> put 2.0 before 1.5.
        doc["axes"][0]["edges"] = {1.0, 2.0, 1.5, 2.5, 3.0, 4.0, 5.0, 7.0, 11.0};
        const TempGrid t(doc, "nonmono");
        CHECK_THROWS_AS(pi0::Binning::load(t.path, grid_b_path()), std::runtime_error);
        CHECK_THROWS_WITH(pi0::Binning::load(t.path, grid_b_path()),
                          ContainsSubstring("strictly increasing"));
    }

    SECTION("a fully reversed axis") {
        nlohmann::json doc = parsed(grid_b_path());
        doc["axes"][1]["edges"] = {1.5, 0.7, 0.4, 0.2, 0.1, 0.0};
        const TempGrid t(doc, "reversed");
        CHECK_THROWS_WITH(pi0::Binning::load(grid_a_path(), t.path),
                          ContainsSubstring("strictly increasing"));
    }

    SECTION("EQUAL adjacent edges") {
        // Rejected too. A zero-width bin can never be filled -- find() is
        // half-open, so [e, e) is empty -- and it would surface downstream as a
        // permanently empty bin with no explanation.
        nlohmann::json doc = parsed(grid_b_path());
        doc["axes"][0]["edges"] = {0.0, 0.2, 0.2, 0.6, 0.8, 1.0};
        const TempGrid t(doc, "dup_edge");
        CHECK_THROWS_WITH(pi0::Binning::load(grid_a_path(), t.path),
                          ContainsSubstring("strictly increasing"));
    }

    SECTION("an edge literal too large for a double") {
        // JSON has no NaN or infinity literal, so a non-finite edge cannot
        // actually be spelled in a grid file -- read_axis()'s std::isfinite
        // guard is defence in depth for a Grid1D built in code, not a path
        // reachable from disk, and this test does NOT claim to drive it.
        //
        // What IS reachable from disk is an edge written as 1e999. nlohmann
        // rejects it with out_of_range (406) -- a json::exception but NOT a
        // json::parse_error. That distinction is the bug this section caught:
        // load() originally caught only parse_error, so 406 escaped as a
        // json::exception, and json is PRIVATE to Binning.cpp. Binning.hpp
        // promises std::runtime_error and does not include nlohmann, so a
        // caller following the header had no type to catch and the process
        // would die unhandled while trying to report a bad grid.
        std::string raw = parsed(grid_a_path()).dump(2);
        const std::string needle = "0.25,";
        const std::size_t at = raw.find(needle);
        REQUIRE(at != std::string::npos);
        raw.replace(at, needle.size() - 1, "1e999");
        const std::string path = std::string(PI0_TEST_TMPDIR) + "/grid_overflow.json";
        {
            std::ofstream out(path);
            REQUIRE(out.good());
            out << raw;
        }
        // std::runtime_error specifically -- that is the header's contract, and
        // the whole point is that no nlohmann type crosses this boundary.
        CHECK_THROWS_AS(pi0::Binning::load(path, grid_b_path()), std::runtime_error);
        CHECK_THROWS_WITH(pi0::Binning::load(path, grid_b_path()), ContainsSubstring(path));
        std::remove(path.c_str());
    }
}

TEST_CASE("Binning::load throws on FEWER THAN 2 edges", "[binning]") {
    SECTION("one edge") {
        nlohmann::json doc = parsed(grid_a_path());
        doc["axes"][0]["edges"] = {1.0};
        const TempGrid t(doc, "one_edge");
        CHECK_THROWS_WITH(pi0::Binning::load(t.path, grid_b_path()),
                          ContainsSubstring("at least 2"));
    }

    SECTION("no edges") {
        nlohmann::json doc = parsed(grid_b_path());
        doc["axes"][1]["edges"] = nlohmann::json::array();
        const TempGrid t(doc, "no_edge");
        CHECK_THROWS_WITH(pi0::Binning::load(grid_a_path(), t.path),
                          ContainsSubstring("at least 2"));
    }
}

TEST_CASE("Binning::load throws on a malformed document", "[binning]") {
    SECTION("not valid JSON") {
        const std::string path = std::string(PI0_TEST_TMPDIR) + "/grid_notjson.json";
        {
            std::ofstream out(path);
            REQUIRE(out.good());
            out << "{ this is not json";
        }
        CHECK_THROWS_WITH(pi0::Binning::load(path, grid_b_path()),
                          ContainsSubstring("not valid JSON"));
        std::remove(path.c_str());
    }

    SECTION("valid JSON but not an object") {
        const std::string path = std::string(PI0_TEST_TMPDIR) + "/grid_array.json";
        {
            std::ofstream out(path);
            REQUIRE(out.good());
            out << "[1, 2, 3]";
        }
        CHECK_THROWS_AS(pi0::Binning::load(path, grid_b_path()), std::runtime_error);
        std::remove(path.c_str());
    }

    SECTION("the wrong number of axes") {
        nlohmann::json doc = parsed(grid_a_path());
        doc["axes"].erase(1);
        const TempGrid t(doc, "one_axis");
        CHECK_THROWS_WITH(pi0::Binning::load(t.path, grid_b_path()),
                          ContainsSubstring("exactly 2"));
    }

    SECTION("edges that are not numbers") {
        nlohmann::json doc = parsed(grid_b_path());
        doc["axes"][0]["edges"] = {0.0, "0.2", 0.4, 0.6, 0.8, 1.0};
        const TempGrid t(doc, "string_edge");
        CHECK_THROWS_WITH(pi0::Binning::load(grid_a_path(), t.path),
                          ContainsSubstring("is not a number"));
    }

    SECTION("SWAPPED axes") {
        // Axis order is load-bearing: it is the slow-then-fast order of the
        // flat index formula the downstream Python decodes against. A swapped
        // pair binds and indexes perfectly cleanly -- and labels every result
        // with the wrong variable. Refuse it.
        nlohmann::json doc = parsed(grid_b_path());
        std::swap(doc["axes"][0], doc["axes"][1]);
        const TempGrid t(doc, "swapped");
        CHECK_THROWS_WITH(pi0::Binning::load(grid_a_path(), t.path),
                          ContainsSubstring("must be 'z'"));
    }

    SECTION("grid A's file handed in as grid B") {
        // The same guard catches the obvious argument-order slip.
        CHECK_THROWS_AS(pi0::Binning::load(grid_b_path(), grid_a_path()), std::runtime_error);
    }
}

// ===========================================================================
// provenance_hash
// ===========================================================================

TEST_CASE("provenance_hash is stable across loads", "[binning]") {
    const pi0::Binning b1 = pi0::Binning::load(grid_a_path(), grid_b_path());
    const pi0::Binning b2 = pi0::Binning::load(grid_a_path(), grid_b_path());

    CHECK(b1.provenance_hash() == b2.provenance_hash());
    CHECK(b1.provenance_hash().size() == 16);
    CHECK(b1.provenance_hash().find_first_not_of("0123456789abcdef") == std::string::npos);
}

TEST_CASE("provenance_hash ignores comments and provenance", "[binning]") {
    // The hash covers the axis NAMES and the EDGES, and nothing else. Rewording
    // a comment must not invalidate an existing result -- if it did, nobody
    // would ever improve a comment.
    const std::string base = shipped_binning().provenance_hash();

    nlohmann::json doc = parsed(grid_a_path());
    doc["_comment"] = "a completely different comment";
    doc["axes"][0]["_comment"] = "also different";
    doc["provenance"]["source"] = "some_scan_v3";
    doc["provenance"]["n_events"] = 123456789;
    doc["_generated"] = "2099-01-01";
    const TempGrid t(doc, "recomment");

    const pi0::Binning b = pi0::Binning::load(t.path, grid_b_path());
    CHECK(b.provenance_hash() == base);
}

TEST_CASE("provenance_hash CHANGES when an edge changes", "[binning]") {
    const std::string base = shipped_binning().provenance_hash();

    SECTION("a Grid A edge moves") {
        nlohmann::json doc = parsed(grid_a_path());
        doc["axes"][0]["edges"] = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 11.5};  // 11.0 -> 11.5
        const TempGrid t(doc, "edge_moved");
        const pi0::Binning b = pi0::Binning::load(t.path, grid_b_path());
        CHECK(b.provenance_hash() != base);
    }

    SECTION("a Grid B edge moves by ONE ULP") {
        // The hash serialises at 17 significant digits, so it separates doubles
        // that differ in the last bit. Two grids that are not bit-identical are
        // not the same binning, and a result stamped with one must not appear
        // to have come from the other.
        nlohmann::json doc = parsed(grid_b_path());
        doc["axes"][1]["edges"] = {0.0, 0.1, 0.2, 0.4, std::nextafter(0.7, 1.0), 1.5};
        const TempGrid t(doc, "edge_ulp");
        const pi0::Binning b = pi0::Binning::load(grid_a_path(), t.path);
        CHECK(b.provenance_hash() != base);
    }

    SECTION("an edge is added") {
        nlohmann::json doc = parsed(grid_b_path());
        doc["axes"][0]["edges"] = {0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0};  // 5 -> 6 bins
        const TempGrid t(doc, "edge_added");
        const pi0::Binning b = pi0::Binning::load(grid_a_path(), t.path);
        CHECK(b.provenance_hash() != base);
        CHECK(b.n4d() != shipped_binning().n4d());
    }

    SECTION("an axis is renamed") {
        // Hand-built rather than loaded: load() enforces the axis names, so a
        // renamed axis cannot come off disk. The hash must separate them anyway
        // -- provenance_hash() is what a downstream reader trusts to say "these
        // are the same grids", and it should not depend on the loader's
        // validation to be true.
        pi0::Binning b = shipped_binning();
        REQUIRE(b.provenance_hash() == base);
        b.A.x.name = "Q2";  // same edges, different label
        CHECK(b.provenance_hash() != base);
    }
}

TEST_CASE("provenance_hash separates grid A from grid B", "[binning]") {
    // Two grids with the SAME multiset of axes in a different arrangement must
    // hash differently: the A/B roles are part of the geometry, not a naming
    // convention.
    const pi0::Binning& ship = shipped_binning();

    pi0::Binning swapped;
    swapped.A = ship.B;
    swapped.B = ship.A;
    CHECK(swapped.provenance_hash() != ship.provenance_hash());

    // And a Binning assembled by hand from the shipped grids hashes the same as
    // the shipped one -- the hash is a function of the geometry alone, not of
    // how the object was built.
    pi0::Binning rebuilt;
    rebuilt.A = ship.A;
    rebuilt.B = ship.B;
    CHECK(rebuilt.provenance_hash() == ship.provenance_hash());
}
