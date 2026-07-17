// Unit tests for pi0::NeumaierSum and the fixed-partition helpers.
//
// The property these guard is the one the whole multi-threaded design rests on:
// a partitioned, compensated reduction gives the SAME BITS regardless of how the
// stream is split -- so a run on 8 threads equals a run on 1, and equals no
// threads at all. If that ever stops holding, Stage B's abscissa stops being
// reproducible, which is the exact defect the stage exists to prevent.
//
// Reference values are computed independently here, never by calling the code
// under test in a way that would make the test pass by construction.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>
#include <random>
#include <vector>

#include "core/Summation.hpp"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinULP;
using pi0::NeumaierSum;
using pi0::partition_count;
using pi0::partition_range;

namespace {

/// A partitioned Neumaier reduction of `xs` into `p` fixed partitions, merged in
/// partition-index order -- exactly what Stage B does per cell. Factored here so
/// the tests can vary p and assert the result does not move.
double partitioned_sum(const std::vector<double>& xs, std::size_t p) {
    const std::size_t n = xs.size();
    std::vector<NeumaierSum> parts(p);
    for (std::size_t i = 0; i < p; ++i) {
        const auto r = partition_range(i, p, n);
        for (std::size_t k = r.begin; k < r.end; ++k) parts[i].add(xs[k]);
    }
    NeumaierSum total;
    for (const auto& part : parts) total.merge(part);
    return total.value();
}

}  // namespace

TEST_CASE("NeumaierSum recovers precision a naive double sum loses", "[summation]") {
    // The textbook catastrophe: a big value, then many small ones whose running
    // contribution falls below the big value's ULP. Naive summation drops them;
    // the compensated sum keeps them. 1e16 + 1 added 1e6 times: the exact answer
    // is 1e16 + 1e6, but 1e16 + 1.0 == 1e16 in double.
    const double big = 1e16;
    NeumaierSum s;
    s.add(big);
    double naive = big;
    for (int i = 0; i < 1'000'000; ++i) {
        s.add(1.0);
        naive += 1.0;
    }
    CHECK(naive == big);  // the naive sum lost every single 1.0
    CHECK_THAT(s.value(), WithinAbs(big + 1'000'000.0, 0.5));
}

TEST_CASE("merge carries the compensation, not just the running sum", "[summation]") {
    // Two partials, each individually lossy, merged. If merge() folded only the
    // running `sum` and dropped `c`, the recovered precision would be lost at the
    // seam. Build two partials that each need their compensation.
    NeumaierSum a, b;
    a.add(1e16);
    b.add(1e16);
    for (int i = 0; i < 500'000; ++i) {
        a.add(1.0);
        b.add(1.0);
    }
    NeumaierSum merged;
    merged.merge(a);
    merged.merge(b);
    CHECK_THAT(merged.value(), WithinAbs(2e16 + 1'000'000.0, 1.0));
}

TEST_CASE("a partitioned reduction is bit-identical for every partition count",
          "[summation][reproducibility]") {
    // THE test. This models exactly what changing --threads does: the same data,
    // split into a different number of fixed partitions, reduced in index order.
    // The result must not move by a single bit.
    std::mt19937_64 rng(0xC0FFEE);  // fixed seed: the data is the same every run
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    std::vector<double> xs(100'000);
    for (auto& x : xs) x = dist(rng);

    const double ref = partitioned_sum(xs, 1);
    for (std::size_t p : {2u, 3u, 4u, 7u, 8u, 16u, 64u, 997u}) {
        INFO("partition count = " << p);
        // Bit-for-bit: the design promises byte-identical output across thread
        // counts, so this is == and not an approximate match.
        CHECK(partitioned_sum(xs, p) == ref);
    }
}

TEST_CASE("partition_count depends only on n and the target", "[summation]") {
    CHECK(partition_count(0) == 1);       // never zero: a reduction needs >= 1 bucket
    CHECK(partition_count(1) == 1);
    CHECK(partition_count(200000, 200000) == 1);
    CHECK(partition_count(200001, 200000) == 2);
    CHECK(partition_count(1'000'000, 200000) == 5);
    CHECK(partition_count(5, 200000) == 1);
    // Never more partitions than items.
    CHECK(partition_count(3, 1) == 3);
}

TEST_CASE("partition_range tiles [0, n) with no gaps or overlaps", "[summation]") {
    for (std::size_t n : {0u, 1u, 5u, 100u, 1000u}) {
        for (std::size_t count : {1u, 2u, 3u, 7u, 13u}) {
            if (count > n && n != 0) continue;  // partition_count would not produce this
            const std::size_t c = (n == 0) ? 1 : count;
            std::size_t covered = 0;
            std::size_t prev_end = 0;
            for (std::size_t i = 0; i < c; ++i) {
                const auto r = partition_range(i, c, n);
                CHECK(r.begin == prev_end);   // contiguous with the previous
                CHECK(r.end >= r.begin);      // never inverted
                covered += r.end - r.begin;
                prev_end = r.end;
            }
            CHECK(prev_end == n);   // the last partition ends exactly at n
            CHECK(covered == n);    // every item covered exactly once
        }
    }
}

TEST_CASE("partition_range splits as evenly as integer division allows", "[summation]") {
    // 100 items into 7 partitions: 100 = 7*14 + 2, so the first two get 15 and
    // the rest 14. This is what keeps every core's slice about the same size.
    const std::size_t n = 100, count = 7;
    std::vector<std::size_t> lens;
    for (std::size_t i = 0; i < count; ++i) {
        const auto r = partition_range(i, count, n);
        lens.push_back(r.end - r.begin);
    }
    CHECK(lens == std::vector<std::size_t>{15, 15, 14, 14, 14, 14, 14});
}

TEST_CASE("an empty stream reduces to zero without special-casing", "[summation]") {
    std::vector<double> xs;
    CHECK(partitioned_sum(xs, partition_count(0)) == 0.0);
}
