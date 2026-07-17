// Unit tests for the frozen event-mixing donor pool (src/stageB_bin).
//
// These tests exist to make five specific properties of the SUPERSEDED mixing
// scheme impossible to reintroduce. All five shipped:
//
//   1. A POOL THAT DEPENDED ON THREAD SCHEDULING. The old rolling FIFO's
//      contents at event k depended on how the record reader happened to
//      interleave, so the background was not reproducible -- it needed 32 shards
//      and a shared_mutex to be merely SAFE, which was never the same property.
//      Here: the same offers and the same seed must give BIT-IDENTICAL donors.
//   2. AN UNSEEDED SAMPLER. The old kd-tree binning drew its reservoir with no
//      seed at all, so two passes over one dataset gave different edges and the
//      production's binning is permanently unrecoverable. Mixing now samples too,
//      so the seed has to be real: a different seed must give a different pool,
//      or the seed is decorative and the defect is back in a new place.
//   3. A POOL THAT WAS COLD EARLY AND WARM LATE. Depth is uniform by
//      construction here, so the test is that a bin offered fewer photons than
//      the depth keeps ALL of them -- no silent loss at the shallow end.
//   4. SELF-MIXING PREVENTED ONLY BY A MIX-THEN-INSERT ORDERING. donors() throws
//      before freeze(), so there is no phase in which a pool is both readable
//      and still growing.
//   5. SILENTLY POOLING PHOTONS THAT FAILED THEIR OWN EVENT'S e-gamma CUT. That
//      contract is the CALLER'S and no test can prove it from in here -- see
//      DonorPool.hpp. What is testable is the tripwire: offering a photon count
//      that contradicts the bin's own multiplicity class must throw.
//
// Method note, and it matters for #2: every test that asserts randomness fixes
// its seeds as literals. A seeded sampler is deterministic, so a "uniform-ish"
// assertion over fixed seeds is a claim about ONE known sequence and cannot
// flake -- it either passes forever or it caught a real change. There is no
// tolerance tuned to make a flaky test quiet, because there is no flaky test.
//
// The pool grid is read from the REAL config/cuts.json rather than hand-built,
// for the same reason test_selection reads it: a hand-written copy proves only
// that the copy parses. The 224-bin claim in the note (sec:mixing, eq:pool-key)
// is a claim about the SHIPPED file.

#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "config/Cuts.hpp"
#include "stageB_bin/DonorPool.hpp"
#include "stageB_bin/PoolGrid.hpp"

#ifndef PI0_CUTS_JSON
#error "PI0_CUTS_JSON must be defined by the build (absolute path to config/cuts.json)."
#endif

using namespace pi0;

namespace {

/// The shipped pool grid, from the real config/cuts.json.
const PoolGridCuts& shipped_grid() {
    static const Cuts c = Cuts::load(PI0_CUTS_JSON);
    return c.mixing.pool_grid;
}

const MixingCuts& shipped_mixing() {
    static const Cuts c = Cuts::load(PI0_CUTS_JSON);
    return c.mixing;
}

/// A small hand-built grid for the tests that are about the RESERVOIR rather
/// than about the shipped configuration: 2 Q^2 bins x 2 x_B bins x 2 classes = 8.
/// Built in code precisely to prove the claim in meson.build that the donor pool
/// needs no JSON to be exercised.
PoolGridCuts toy_grid() {
    PoolGridCuts g;
    g.q2 = Grid1D{"q2", {1.0, 2.0, 3.0}};
    g.xb = Grid1D{"xb", {0.1, 0.2, 0.3}};
    g.n_photons_classes = {
        MultClass{"1", 1, 1},
        MultClass{">=2", 2, std::nullopt},
    };
    return g;
}

/// A donor photon tagged by an integer, so a reservoir's contents can be read
/// back as "which offers survived". px carries the tag; the rest is filler that
/// keeps the struct honest about being a photon.
DonorPhoton tagged(int tag) {
    DonorPhoton p;
    p.px = static_cast<float>(tag);
    p.py = 0.0f;
    p.pz = 1.0f;
    p.e = 1.0f;
    p.inv_p = 1.0f;
    return p;
}

int tag_of(const DonorPhoton& p) { return static_cast<int>(p.px); }

/// Flatten a multi-line string onto one line for INFO().
///
/// NOT cosmetic. meson's TAP parser rejects any stdout line it cannot parse and
/// reports a FULLY PASSING run as an ERROR -- the same trap tests/meson.build
/// already documents for Catch2's WARN(). report_underfilled() writes a
/// deliberately multi-line report, so INFO-ing it raw turns a green suite red
/// without a single failing assertion.
std::string one_line(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (const char c : s) out += (c == '\n') ? ' ' : c;
    return out;
}

/// Offer `n` singleton events to one bin of a toy pool, in tag order 0..n-1.
/// One photon per offer, so the bin's class must be the "1" class (index 0).
void offer_singletons(DonorPool& pool, int bin, int n) {
    for (int i = 0; i < n; ++i) {
        pool.offer(bin, {tagged(i)});
    }
}

/// Bit-exact comparison of two donor lists. `==` on float is deliberate: the
/// claim under test is BIT-identical reproduction, not numerical agreement, so
/// an epsilon here would defeat the entire point of the test.
bool identical(const std::vector<DonorPhoton>& a, const std::vector<DonorPhoton>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].px != b[i].px || a[i].py != b[i].py || a[i].pz != b[i].pz || a[i].e != b[i].e ||
            a[i].inv_p != b[i].inv_p) {
            return false;
        }
    }
    return true;
}

}  // namespace

// ---------------------------------------------------------------------------
// The shipped configuration
// ---------------------------------------------------------------------------

TEST_CASE("the shipped pool grid is 8 x 7 x 4 = 224 bins", "[donorpool][config]") {
    const PoolGridCuts& g = shipped_grid();

    // Transcribed from config/cuts.json, which is the source of truth. If these
    // ever disagree with it, IT WINS and this file is what changes -- but a
    // disagreement means someone moved the pool grid, and the note's eq:pool-key
    // says 224. That should be a deliberate act, which is why this is here.
    CHECK(g.q2.nbins() == 8);
    CHECK(g.xb.nbins() == 7);
    CHECK(g.n_photons_classes.size() == 4u);
    CHECK(g.n_bins() == 224u);

    CHECK(shipped_mixing().donors_per_bin == 200u);
    CHECK(shipped_mixing().seed_mode == "file_hash");

    // The >=4 class must be OPEN-ENDED. A finite max here would silently drop
    // every event above it out of the pool entirely.
    REQUIRE(g.n_photons_classes.size() == 4u);
    CHECK(g.n_photons_classes[3].label == ">=4");
    CHECK(g.n_photons_classes[3].min == 4);
    CHECK_FALSE(g.n_photons_classes[3].max.has_value());
}

TEST_CASE("Cuts::load reads the mixing block rather than defaulting it",
          "[donorpool][config]") {
    // The block was declared and never read for the whole of this project's life
    // until DonorPool consumed it. A zero here would mean it went back to being
    // parsed into nothing.
    CHECK(shipped_mixing().donors_per_bin > 0u);
    CHECK_FALSE(shipped_mixing().pool_grid.q2.edges.empty());
    CHECK_FALSE(shipped_mixing().pool_grid.xb.edges.empty());
}

// ---------------------------------------------------------------------------
// pool_bin: the grid mapping
// ---------------------------------------------------------------------------

TEST_CASE("pool_bin maps the shipped grid correctly", "[donorpool][binning]") {
    const DonorPool pool(shipped_grid(), 1u, 10u);
    REQUIRE(pool.n_bins() == 224u);

    const int n_xb = 7, n_mult = 4;
    const auto expect = [&](int i_q2, int i_xb, int i_mult) {
        return (i_q2 * n_xb + i_xb) * n_mult + i_mult;
    };

    // Q^2 edges [1, 1.5, 2, 2.5, 3, 4, 5, 7, 11]; x_B edges [.1,.15,.2,.25,.3,.4,.5,.7].
    // The first cell of the grid, 1 photon -> bin 0.
    CHECK(pool.pool_bin(1.0, 0.1, 1) == expect(0, 0, 0));
    CHECK(pool.pool_bin(1.2, 0.12, 1) == 0);

    // Multiplicity is the FASTEST axis.
    CHECK(pool.pool_bin(1.2, 0.12, 2) == expect(0, 0, 1));
    CHECK(pool.pool_bin(1.2, 0.12, 3) == expect(0, 0, 2));
    CHECK(pool.pool_bin(1.2, 0.12, 4) == expect(0, 0, 3));
    // ... and >=4 absorbs everything above it, rather than falling off the grid.
    CHECK(pool.pool_bin(1.2, 0.12, 9) == expect(0, 0, 3));
    CHECK(pool.pool_bin(1.2, 0.12, 500) == expect(0, 0, 3));

    // x_B is the middle axis, Q^2 the slowest.
    CHECK(pool.pool_bin(1.2, 0.22, 1) == expect(0, 2, 0));
    CHECK(pool.pool_bin(2.2, 0.12, 1) == expect(2, 0, 0));
    CHECK(pool.pool_bin(2.2, 0.22, 3) == expect(2, 2, 2));

    // The last cell, and the top edge of BOTH axes lands in the last bin rather
    // than off the grid -- the one place the half-open rule is broken, and it is
    // Grid1D's documented behaviour (the old find_1d_bin() returned -1 there and
    // so dropped every value sitting exactly on the kinematic limit).
    CHECK(pool.pool_bin(11.0, 0.7, 4) == expect(7, 6, 3));
    CHECK(pool.pool_bin(11.0, 0.7, 4) == 223);
    CHECK(pool.pool_bin(10.9, 0.69, 1) == expect(7, 6, 0));

    // Bins are half-open at the bottom: a value ON an interior edge belongs to
    // the bin ABOVE it.
    CHECK(pool.pool_bin(1.5, 0.12, 1) == expect(1, 0, 0));
    CHECK(pool.pool_bin(1.2, 0.15, 1) == expect(0, 1, 0));
}

TEST_CASE("pool_bin returns -1 off the grid and for zero photons",
          "[donorpool][binning]") {
    const DonorPool pool(shipped_grid(), 1u, 10u);

    SECTION("Q^2 out of range") {
        CHECK(pool.pool_bin(0.99, 0.2, 1) == -1);   // below dis.q2_min
        CHECK(pool.pool_bin(11.01, 0.2, 1) == -1);  // above the top edge
    }
    SECTION("x_B out of range -- there is NO x_B cut, so this genuinely happens") {
        CHECK(pool.pool_bin(2.0, 0.09, 1) == -1);
        CHECK(pool.pool_bin(2.0, 0.71, 1) == -1);
    }
    SECTION("zero photons is in no class, whatever the kinematics") {
        // An event with no photons produces no pairs. It must not land in the
        // "1" class by an off-by-one, which is exactly what a `n >= min` test
        // with min = 0 would do.
        CHECK(pool.pool_bin(2.0, 0.2, 0) == -1);
        CHECK(pool.pool_bin(1.0, 0.1, 0) == -1);
        CHECK(pool.pool_bin(11.0, 0.7, 0) == -1);
    }
    SECTION("NaN is off the grid rather than in an arbitrary bin") {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK(pool.pool_bin(nan, 0.2, 1) == -1);
        CHECK(pool.pool_bin(2.0, nan, 1) == -1);
    }
}

// ---------------------------------------------------------------------------
// Determinism -- the property the whole design exists for
// ---------------------------------------------------------------------------

TEST_CASE("the same offers and seed give a BIT-IDENTICAL pool", "[donorpool][determinism]") {
    // 500 offers into a depth-20 reservoir: the reservoir overflows 25x over, so
    // essentially every retained donor got there through a random replacement.
    // If anything in the sampler consulted the clock, the address of an object,
    // or an unseeded engine, these two pools would differ.
    const auto build = [](std::uint64_t seed) {
        DonorPool p(toy_grid(), seed, 20u);
        offer_singletons(p, 0, 500);
        offer_singletons(p, 4, 500);
        p.freeze();
        return p;
    };

    const DonorPool a = build(0xDEADBEEFULL);
    const DonorPool b = build(0xDEADBEEFULL);

    REQUIRE(a.n_bins() == b.n_bins());
    for (std::size_t bin = 0; bin < a.n_bins(); ++bin) {
        const int ib = static_cast<int>(bin);
        INFO("pool bin " << ib);
        CHECK(identical(a.donors(ib), b.donors(ib)));
    }

    // And the pools are non-trivial: a pool of empty bins would pass the above
    // while proving nothing at all.
    CHECK(a.donors(0).size() == 20u);
    CHECK(a.donors(4).size() == 20u);
    CHECK(a.n_filled() == 2u);
}

TEST_CASE("a different seed gives a different pool", "[donorpool][determinism]") {
    // Without this, a sampler that ignored its seed entirely would pass the
    // determinism test above perfectly. That is not hypothetical: the defect the
    // note records is an UNSEEDED reservoir, which is bit-identical across two
    // runs of a single-threaded pass and differs only under load.
    const auto build = [](std::uint64_t seed) {
        DonorPool p(toy_grid(), seed, 20u);
        offer_singletons(p, 0, 500);
        p.freeze();
        return p;
    };

    const DonorPool a = build(1u);
    const DonorPool b = build(2u);

    REQUIRE(a.donors(0).size() == 20u);
    REQUIRE(b.donors(0).size() == 20u);
    CHECK_FALSE(identical(a.donors(0), b.donors(0)));

    // Adjacent seeds too: mt19937_64 seeded from nearby integers produces
    // correlated early output, which is why the seed goes through SplitMix64
    // before it reaches an engine. Seeds 1 and 2 differing is the cheap check
    // that the decorrelation is actually wired in.
    const DonorPool c = build(3u);
    CHECK_FALSE(identical(b.donors(0), c.donors(0)));
}

TEST_CASE("bins do not share an RNG stream", "[donorpool][determinism]") {
    // Per-bin engines mean a bin's reservoir depends only on the sequence offered
    // to THAT bin. So filling bin 0 alone, and filling bin 0 with bin 4's traffic
    // interleaved between every offer, must give bin 0 the same donors. With one
    // shared engine this fails: bin 4's draws would advance bin 0's stream.
    DonorPool alone(toy_grid(), 7u, 20u);
    offer_singletons(alone, 0, 300);
    alone.freeze();

    DonorPool interleaved(toy_grid(), 7u, 20u);
    for (int i = 0; i < 300; ++i) {
        interleaved.offer(0, {tagged(i)});
        interleaved.offer(4, {tagged(1000 + i)});  // bin 4: a different Q^2 cell
    }
    interleaved.freeze();

    CHECK(identical(alone.donors(0), interleaved.donors(0)));
}

// ---------------------------------------------------------------------------
// The reservoir itself
// ---------------------------------------------------------------------------

TEST_CASE("a bin offered fewer photons than the depth keeps ALL of them",
          "[donorpool][reservoir]") {
    // The shallow end. The old FIFO's cold-start meant early events mixed against
    // a nearly empty pool; here a thin bin is thin but must never be LOSSY.
    DonorPool pool(toy_grid(), 42u, 200u);
    offer_singletons(pool, 0, 7);
    pool.freeze();

    const auto& d = pool.donors(0);
    REQUIRE(d.size() == 7u);

    // Every tag present, exactly once, and in offer order -- while the reservoir
    // is filling there is no replacement to make and so nothing to permute.
    for (int i = 0; i < 7; ++i) {
        CHECK(tag_of(d[static_cast<std::size_t>(i)]) == i);
    }
    CHECK(pool.n_offered(0) == 7u);
}

TEST_CASE("a bin offered exactly the depth keeps all of them", "[donorpool][reservoir]") {
    // The boundary: `reservoir.size() < depth` must admit the depth'th photon and
    // take no draw for it. An off-by-one here would silently drop one donor per
    // bin or take a spurious draw that shifts every later replacement.
    DonorPool pool(toy_grid(), 42u, 20u);
    offer_singletons(pool, 0, 20);
    pool.freeze();

    REQUIRE(pool.donors(0).size() == 20u);
    for (int i = 0; i < 20; ++i) {
        CHECK(tag_of(pool.donors(0)[static_cast<std::size_t>(i)]) == i);
    }
}

TEST_CASE("the reservoir never exceeds its depth", "[donorpool][reservoir]") {
    DonorPool pool(toy_grid(), 42u, 20u);
    offer_singletons(pool, 0, 5000);
    pool.freeze();
    CHECK(pool.donors(0).size() == 20u);
    CHECK(pool.n_offered(0) == 5000u);
}

TEST_CASE("multi-photon events contribute every photon", "[donorpool][reservoir]") {
    // Bin 1 is the ">=2" class of the toy grid's first cell. A 3-photon event
    // offers three donors, not one: the reservoir's unit is the photon.
    DonorPool pool(toy_grid(), 42u, 200u);
    pool.offer(1, {tagged(0), tagged(1), tagged(2)});
    pool.offer(1, {tagged(3), tagged(4)});
    pool.freeze();

    CHECK(pool.donors(1).size() == 5u);
    CHECK(pool.n_offered(1) == 5u);
}

TEST_CASE("reservoir sampling is uniform over many offers", "[donorpool][reservoir]") {
    // THE TEST THAT CATCHES A BROKEN SAMPLER. Every classic reservoir bug is a
    // bias: keeping the first k (never replacing), keeping the last k (always
    // replacing), or drawing over the wrong range so early offers are favoured.
    // All of them survive the determinism tests above untouched.
    //
    // 400 pools, seeds 1..400 as literals, each offering 100 tagged photons into
    // a depth-10 reservoir. Expected selections per tag = 400 * 10 / 100 = 40.
    // Fixed seeds -> one fixed known sequence -> a chi^2 that cannot flake.
    constexpr int kTrials = 400;
    constexpr int kOffers = 100;
    constexpr std::size_t kDepth = 10;
    constexpr double kExpected = kTrials * static_cast<double>(kDepth) / kOffers;

    std::vector<int> hits(kOffers, 0);
    for (int t = 1; t <= kTrials; ++t) {
        DonorPool pool(toy_grid(), static_cast<std::uint64_t>(t), kDepth);
        offer_singletons(pool, 0, kOffers);
        pool.freeze();

        const auto& d = pool.donors(0);
        REQUIRE(d.size() == kDepth);

        // A reservoir must never hold a duplicate: each offer is placed at most
        // once, so a repeat means an index was written twice.
        std::set<int> seen;
        for (const DonorPhoton& p : d) {
            const int tag = tag_of(p);
            REQUIRE(tag >= 0);
            REQUIRE(tag < kOffers);
            CHECK(seen.insert(tag).second);
            hits[static_cast<std::size_t>(tag)]++;
        }
    }

    double chi2 = 0.0;
    int zero_hits = 0;
    for (const int h : hits) {
        const double d = h - kExpected;
        chi2 += d * d / kExpected;
        if (h == 0) ++zero_hits;
    }

    INFO("chi2 = " << chi2 << " over " << (kOffers - 1) << " dof, expected ~" << (kOffers - 1));

    // A LOOSE bound, on purpose. 99 dof has mean 99 and sd ~14, so a healthy
    // sampler lands near 99 and this passes with enormous margin. The bound is
    // not trying to measure uniformity to a p-value -- it is trying to catch a
    // sampler that is WRONG, and every way of being wrong here is catastrophic:
    // "keep the first 10" pins 90 tags at zero hits and sends chi2 past 35000.
    CHECK(chi2 < 300.0);

    // The blunt instrument, which no biased sampler survives: every tag must get
    // picked at least once across 400 trials.
    CHECK(zero_hits == 0);

    // And the late offers are not starved -- the specific signature of drawing
    // over a stale count. Compare the first and last deciles.
    int first_ten = 0, last_ten = 0;
    for (int i = 0; i < 10; ++i) first_ten += hits[static_cast<std::size_t>(i)];
    for (int i = kOffers - 10; i < kOffers; ++i) last_ten += hits[static_cast<std::size_t>(i)];
    INFO("first decile " << first_ten << " vs last decile " << last_ten);
    CHECK(first_ten > 250);
    CHECK(last_ten > 250);
    CHECK(first_ten < 550);
    CHECK(last_ten < 550);
}

// ---------------------------------------------------------------------------
// The frozen lifecycle
// ---------------------------------------------------------------------------

TEST_CASE("donors() throws before freeze()", "[donorpool][lifecycle]") {
    DonorPool pool(toy_grid(), 1u, 10u);
    offer_singletons(pool, 0, 5);

    CHECK_FALSE(pool.frozen());
    CHECK_THROWS_AS(pool.donors(0), std::logic_error);
    CHECK_THROWS_WITH(pool.donors(0), Catch::Matchers::ContainsSubstring("not frozen"));

    pool.freeze();
    CHECK(pool.frozen());
    CHECK_NOTHROW(pool.donors(0));
}

TEST_CASE("offer() throws after freeze()", "[donorpool][lifecycle]") {
    // The invariant the old code left to a mix-then-insert ORDERING. Here it is
    // not an ordering: after freeze() there is no insert to reorder.
    DonorPool pool(toy_grid(), 1u, 10u);
    pool.freeze();
    CHECK_THROWS_AS(pool.offer(0, {tagged(1)}), std::logic_error);
}

TEST_CASE("freeze() is idempotent and does not disturb the pool",
          "[donorpool][lifecycle]") {
    DonorPool pool(toy_grid(), 1u, 10u);
    offer_singletons(pool, 0, 100);
    pool.freeze();

    const std::vector<DonorPhoton> after_first = pool.donors(0);
    CHECK_NOTHROW(pool.freeze());
    CHECK(identical(pool.donors(0), after_first));
}

TEST_CASE("out-of-range bins throw rather than being clamped",
          "[donorpool][lifecycle]") {
    DonorPool pool(toy_grid(), 1u, 10u);
    const int n = static_cast<int>(pool.n_bins());
    REQUIRE(n == 8);

    // -1 is what pool_bin() returns for an off-grid event. It is not a bin, and
    // a caller that forgets to test for it must be told, not quietly given bin 0
    // or the last bin.
    CHECK_THROWS_AS(pool.offer(-1, {tagged(1)}), std::out_of_range);
    CHECK_THROWS_AS(pool.offer(n, {tagged(1)}), std::out_of_range);
    CHECK_THROWS_AS(pool.n_offered(-1), std::out_of_range);

    pool.freeze();
    CHECK_THROWS_AS(pool.donors(-1), std::out_of_range);
    CHECK_THROWS_AS(pool.donors(n), std::out_of_range);
}

// ---------------------------------------------------------------------------
// The e-gamma tripwire
// ---------------------------------------------------------------------------

TEST_CASE("offer() rejects a photon count contradicting the bin's class",
          "[donorpool][egamma]") {
    // THE TRIPWIRE, not a proof. The e-gamma filtering contract belongs to the
    // caller and cannot be checked from in here -- a photon that failed the cut
    // looks exactly like one that passed. What IS checkable is the disagreement
    // it produces when someone filters the photons but computes the multiplicity
    // from the unfiltered list: the bin then says "1 photon" while two arrive.
    DonorPool pool(toy_grid(), 1u, 10u);

    // Bin 0 is the "1" class of the toy grid's cell (0,0); bin 1 is that cell's
    // ">=2" class. Q^2 = 2.5 is the NEXT Q^2 cell (the toy edges are {1,2,3},
    // not the shipped grid's), so it is bin (1*2 + 0)*2 + 0 = 4.
    REQUIRE(pool.pool_bin(1.2, 0.12, 1) == 0);
    REQUIRE(pool.pool_bin(1.2, 0.12, 2) == 1);
    REQUIRE(pool.pool_bin(2.5, 0.15, 1) == 4);

    SECTION("too many photons for the class") {
        CHECK_THROWS_AS(pool.offer(0, {tagged(1), tagged(2)}), std::invalid_argument);
        CHECK_THROWS_WITH(pool.offer(0, {tagged(1), tagged(2)}),
                          Catch::Matchers::ContainsSubstring("e-gamma"));
    }
    SECTION("too few photons for the class") {
        CHECK_THROWS_AS(pool.offer(1, {tagged(1)}), std::invalid_argument);
    }
    SECTION("an empty event is in no class and cannot be offered anywhere") {
        CHECK_THROWS_AS(pool.offer(0, {}), std::invalid_argument);
        CHECK_THROWS_AS(pool.offer(1, {}), std::invalid_argument);
    }
    SECTION("the matching count is accepted") {
        CHECK_NOTHROW(pool.offer(0, {tagged(1)}));
        CHECK_NOTHROW(pool.offer(1, {tagged(1), tagged(2)}));
        CHECK_NOTHROW(pool.offer(1, {tagged(1), tagged(2), tagged(3)}));  // >=2 is open-ended
    }
}

// ---------------------------------------------------------------------------
// report_underfilled -- the old code's best QA feature
// ---------------------------------------------------------------------------

TEST_CASE("report_underfilled names bins below depth", "[donorpool][qa]") {
    DonorPool pool(shipped_grid(), 1u, 50u);

    // Fill one bin to depth, leave one thin, and leave the other 222 empty.
    const int full = pool.pool_bin(1.2, 0.12, 1);
    const int thin = pool.pool_bin(2.2, 0.22, 1);
    REQUIRE(full >= 0);
    REQUIRE(thin >= 0);
    REQUIRE(full != thin);

    offer_singletons(pool, full, 400);
    offer_singletons(pool, thin, 3);
    pool.freeze();

    std::ostringstream os;
    pool.report_underfilled(os);
    const std::string out = os.str();
    INFO(one_line(out));

    // The thin bin is named, with what it holds AND what it was offered -- the
    // second number is what makes it diagnosable rather than merely reported.
    CHECK_THAT(out, Catch::Matchers::ContainsSubstring("UNDERFILLED"));
    CHECK_THAT(out, Catch::Matchers::ContainsSubstring("bin " + std::to_string(thin)));
    CHECK_THAT(out, Catch::Matchers::ContainsSubstring("3/50"));

    // Empty bins are called out separately: no donors at all is a different
    // problem from a sparse bin, and must not be lost in a list of near-misses.
    CHECK_THAT(out, Catch::Matchers::ContainsSubstring("EMPTY"));
    CHECK_THAT(out, Catch::Matchers::ContainsSubstring("222 bin(s) received NO photons"));

    // The bin that reached depth is NOT named. A QA report that warns about
    // healthy bins is a report nobody reads.
    CHECK_THAT(out, !Catch::Matchers::ContainsSubstring("bin " + std::to_string(full) + " "));

    CHECK(pool.n_filled() == 2u);
    CHECK(pool.n_bins() == 224u);
}

TEST_CASE("report_underfilled is quiet when every bin reached depth",
          "[donorpool][qa]") {
    DonorPool pool(toy_grid(), 1u, 2u);
    for (int bin = 0; bin < 8; ++bin) {
        // Even bins are the "1" class, odd bins the ">=2" class.
        if (bin % 2 == 0) {
            offer_singletons(pool, bin, 50);
        } else {
            for (int i = 0; i < 25; ++i) pool.offer(bin, {tagged(i), tagged(i + 100)});
        }
    }
    pool.freeze();

    std::ostringstream os;
    pool.report_underfilled(os);
    const std::string out = os.str();
    INFO(one_line(out));

    CHECK(pool.n_filled() == 8u);
    CHECK_THAT(out, Catch::Matchers::ContainsSubstring("all 8 bins reached depth"));
    CHECK_THAT(out, !Catch::Matchers::ContainsSubstring("WARNING"));
}

TEST_CASE("report_underfilled throws before freeze()", "[donorpool][qa]") {
    DonorPool pool(toy_grid(), 1u, 10u);
    std::ostringstream os;
    CHECK_THROWS_AS(pool.report_underfilled(os), std::logic_error);
}

// ---------------------------------------------------------------------------
// Grid validation
// ---------------------------------------------------------------------------

TEST_CASE("the constructor rejects an unusable grid", "[donorpool][validation]") {
    SECTION("zero depth") {
        CHECK_THROWS_AS(DonorPool(toy_grid(), 1u, 0u), std::invalid_argument);
    }
    SECTION("no multiplicity classes") {
        PoolGridCuts g = toy_grid();
        g.n_photons_classes.clear();
        CHECK_THROWS_AS(DonorPool(g, 1u, 10u), std::invalid_argument);
    }
    SECTION("a class starting at zero photons") {
        PoolGridCuts g = toy_grid();
        g.n_photons_classes[0].min = 0;
        CHECK_THROWS_AS(DonorPool(g, 1u, 10u), std::invalid_argument);
    }
    SECTION("an open-ended class that is not last -- every later class is dead") {
        PoolGridCuts g = toy_grid();
        g.n_photons_classes = {
            MultClass{"open", 1, std::nullopt},
            MultClass{"unreachable", 5, 9},
        };
        CHECK_THROWS_AS(DonorPool(g, 1u, 10u), std::invalid_argument);
    }
    SECTION("overlapping classes") {
        PoolGridCuts g = toy_grid();
        g.n_photons_classes = {
            MultClass{"1-3", 1, 3},
            MultClass{"2+", 2, std::nullopt},
        };
        CHECK_THROWS_AS(DonorPool(g, 1u, 10u), std::invalid_argument);
    }
    SECTION("non-monotonic edges -- Grid1D does not check its own") {
        // upper_bound over unsorted edges does not crash; it is merely
        // meaningless, so photons would be filed into the wrong pool bins and
        // every mixed spectrum would be quietly wrong. Nothing else validates
        // this grid: Binning::load() never sees it.
        PoolGridCuts g = toy_grid();
        g.q2 = Grid1D{"q2", {3.0, 2.0, 1.0}};
        CHECK_THROWS_AS(DonorPool(g, 1u, 10u), std::invalid_argument);
    }
    SECTION("duplicate edges make a zero-width bin") {
        PoolGridCuts g = toy_grid();
        g.q2 = Grid1D{"q2", {1.0, 2.0, 2.0}};
        CHECK_THROWS_AS(DonorPool(g, 1u, 10u), std::invalid_argument);
    }
    SECTION("a single edge is not an axis") {
        PoolGridCuts g = toy_grid();
        g.q2 = Grid1D{"q2", {1.0}};
        CHECK_THROWS_AS(DonorPool(g, 1u, 10u), std::invalid_argument);
    }
}
