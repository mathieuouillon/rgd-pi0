#pragma once

/// \file Summation.hpp
/// \brief Order-stable floating-point summation for the kinematic sums, so that
///        a multi-threaded run gives the SAME bits as a single-threaded one.
///
/// WHY THIS EXISTS. Stage B accumulates sum_q2, sum_xb, sum_z, sum_pt2 per (4D
/// bin, m_gg bin). Those sums ARE the count-weighted abscissae -- the whole
/// reason Stage B exists (note sec:binning-caveat): a bin is reported at
/// sum_X / n_same, not at the geometric box centre. A plain `double += v` over N
/// events is exact-in-intent but its rounding depends on the ORDER the values
/// arrive, so splitting the events across T threads and adding the partials back
/// would make sum_q2 depend on T. The provenance used to state, flatly, that the
/// stage is "single-threaded, deliberately" for exactly this reason.
///
/// The fix has two independent parts, and BOTH are needed:
///
///  1. NeumaierSum -- a compensated accumulator (the Neumaier refinement of
///     Kahan) that tracks the rounding error the naive sum throws away. It makes
///     the result far less sensitive to order, and its merge() combines two
///     partial sums while carrying both corrections.
///
///  2. A FIXED partition count, decided by the data and NOT by the thread count.
///     Compensation alone still leaves the last ulp depending on how the stream
///     was split; if the split itself is a fixed function of the data, the result
///     is bit-identical for ANY thread count -- including one thread, and
///     including no threads at all. That last point is the design's backbone:
///     the sums never see a thread boundary, so parallelism CANNOT change them,
///     and correctness can be verified with the threads turned off.
///
/// The one-time cost is that these numbers differ, in their last bits, from the
/// old naive sequential sum. That shift is to the MORE accurate value, it happens
/// once, and it is detectable -- the config hash and the Stage B code version are
/// both stamped into every output.
///
/// No dependency on ROOT or threads: pure arithmetic, unit-tested without a data
/// file, like the rest of src/core.

#include <cmath>
#include <cstddef>
#include <vector>

namespace pi0 {

/// A running sum that also tracks the low-order bits ordinary addition discards.
///
/// Neumaier's variant of Kahan summation: it handles the case where the new term
/// is larger in magnitude than the running total, which plain Kahan mishandles.
/// value() returns sum + correction; read it once at the end, not per add.
struct NeumaierSum {
    double sum = 0.0;
    double c = 0.0;  ///< the compensation: the sum of the bits that did not fit.

    /// Add one value.
    void add(double v) {
        const double t = sum + v;
        // The compensation captures whichever operand's low bits were lost, which
        // depends on which is larger -- this branch is the whole difference from
        // plain Kahan, and the reason a big value added to a small running sum
        // (or vice versa) stays accurate.
        if (std::abs(sum) >= std::abs(v)) {
            c += (sum - t) + v;
        } else {
            c += (v - t) + sum;
        }
        sum = t;
    }

    /// Fold another partial sum into this one, carrying ITS compensation too.
    ///
    /// This is what makes a partitioned reduction exact-to-Neumaier: merging is
    /// not `add(other.value())` (which would round other's correction back into a
    /// single term), but adding other's running sum AND its correction as two
    /// separate contributions.
    void merge(const NeumaierSum& other) {
        add(other.sum);
        add(other.c);
    }

    /// The compensated total. Call once, at the end.
    [[nodiscard]] double value() const { return sum + c; }
};

/// How many fixed partitions to split N items into for an order-stable reduction.
///
/// The count depends ONLY on N (and a target partition size), never on the thread
/// count or the machine, so the boundaries -- and therefore the per-partition
/// sums and their merge order -- are a pure function of the data. `target` is the
/// approximate number of items per partition; the returned count is at least 1
/// and never exceeds N (an empty tail partition would just be a no-op, but there
/// is no reason to create one).
///
/// The default target is deliberately coarse: partitions exist for reproducible
/// reduction first and parallelism second, so a few thousand items each is plenty
/// -- enough partitions to keep every core busy on a real file, few enough that
/// the merge is negligible.
[[nodiscard]] inline std::size_t partition_count(std::size_t n, std::size_t target = 200000) {
    if (n == 0) return 1;
    if (target == 0) target = 1;
    const std::size_t p = (n + target - 1) / target;
    return p < 1 ? 1 : (p > n ? n : p);
}

/// The half-open [begin, end) entry range of partition `i` of `count`, over
/// `n` items. Contiguous, gap-free, and covering [0, n): partition i gets the
/// items the earlier partitions did not, so a reduction over the partitions in
/// index order visits every item exactly once, in the original order WITHIN each
/// partition. The split is as even as integer division allows; the first
/// `n % count` partitions get one extra item.
struct PartitionRange {
    std::size_t begin;
    std::size_t end;
};

[[nodiscard]] inline PartitionRange partition_range(std::size_t i, std::size_t count, std::size_t n) {
    const std::size_t base = n / count;
    const std::size_t rem = n % count;
    // The first `rem` partitions are one longer; offsets accumulate that.
    const std::size_t begin = i * base + (i < rem ? i : rem);
    const std::size_t len = base + (i < rem ? 1 : 0);
    return {begin, begin + len};
}

}  // namespace pi0
