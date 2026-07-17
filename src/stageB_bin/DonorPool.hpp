#pragma once

/// \file DonorPool.hpp
/// \brief The FROZEN donor pool for event mixing: a seeded, deterministic
///        pre-pass builds it, and after freeze() mixing is a PURE FUNCTION of
///        (event, pool).
///
/// ---------------------------------------------------------------------------
/// WHY THIS IS NOT THE OLD ROLLING FIFO
/// ---------------------------------------------------------------------------
/// The superseded analysis kept a rolling FIFO of depth 50 events per pool bin,
/// mixed the current event against the pool, then inserted it. That works, and
/// four things were wrong with it:
///
///   * NOT REPRODUCIBLE. Under multithreading the pool contents at event k
///     depend on record scheduling, so the background depended on how the
///     threads happened to interleave. It needed 32 shards and a shared_mutex to
///     be merely SAFE -- and safety was never the same thing as reproducibility.
///   * TIME-DEPENDENT WITHIN A FILE. The pool is COLD for early events and WARM
///     for late ones, so the background an event sees depends on where in the
///     file it sits -- in an analysis whose largest unevaluated systematic is
///     run-to-run drift.
///   * SELF-MIXING PREVENTED ONLY BY ORDERING. "Mix then insert" is an invariant
///     a maintainer can break with a two-line reorder, and no test would notice.
///   * IT POOLED PHOTONS THAT FAILED THEIR OWN EVENT'S e-gamma CUT. All of the
///     current event's photons went in, so a photon rejected from the same-event
///     spectrum still contributed to the mixed one.
///
/// This design fixes all four, and the fixes are structural rather than
/// procedural:
///
///   * BIT-FOR-BIT REPRODUCIBLE. Same (file, seed, donors_per_bin) -> same pool,
///     always. The RNG is seeded and per-bin; nothing consults the clock, the
///     thread id, or global rand.
///   * UNIFORM IN TIME. Every bin is sampled across the WHOLE file before any
///     mixing happens, so there is no cold start and no warm tail.
///   * ZERO LOCKS. After freeze() the pool is const, so the mixing pass shares it
///     across threads with no shards and no mutex. MT becomes trivial because
///     there is nothing to synchronise, not because the locking is good.
///   * SELF-MIXING IS IMPOSSIBLE BY CONSTRUCTION. The pool is complete before the
///     first mixed pair exists, so "mix then insert" is not an ordering to
///     preserve -- there is no insert. (One residual, stated below.)
///   * MIXED STATISTICS IS A KNOB. `donors_per_bin` sets it directly, rather than
///     it emerging from how fast a FIFO warmed up.
///
/// ---------------------------------------------------------------------------
/// THE ONE CONTRACT THIS CLASS CANNOT ENFORCE -- READ THIS
/// ---------------------------------------------------------------------------
/// *** ONLY OFFER PHOTONS THAT PASSED THEIR OWN EVENT'S e-gamma CUT. ***
///
/// The caller filters. Stage A does NOT apply the e-gamma cut -- it ships
/// `g_e_gamma_deg` per photon and leaves the cut to the consumer -- so the
/// pre-pass MUST drop photons failing `pairing.e_gamma_min_angle_deg` BEFORE
/// calling offer(). Violating this silently reintroduces the old defect in full:
/// a photon rejected from the same-event spectrum would contribute to the mixed
/// one, distorting the background shape near the electron direction. It is
/// silent because the result is a perfectly ordinary-looking spectrum. No test
/// in this project can catch it for you, because from inside offer() a photon
/// that failed the cut is indistinguishable from one that passed.
///
/// offer() does check ONE consequence of the contract: that `event_photons.size()`
/// falls in the multiplicity class encoded in `pool_bin`. That is a tripwire, NOT
/// a proof -- it catches a caller who filtered the photons but computed the class
/// from the unfiltered count (or the reverse), and it catches nothing at all if
/// you get both wrong the same way.
///
/// THE ASYMMETRY IS CORRECT AND IS NOT THIS DEFECT. For a mixed pair only the
/// CURRENT event's photon is tested against the CURRENT event's electron; the
/// donor is never re-tested against it. That is right by construction: the donor
/// was already tested against its OWN event's electron, which is the physically
/// meaningful reference. Do not "fix" it.
///
/// ---------------------------------------------------------------------------
/// THREADING
/// ---------------------------------------------------------------------------
/// TWO PHASES, TWO RULES, AND THEY ARE OPPOSITE:
///
///   BUILD  (offer): NOT thread-safe, and NOT MERELY AS AN OMISSION. Reservoir
///          sampling is order-dependent by nature: the pool is a function of the
///          SEQUENCE of offers, so a parallel pre-pass would make it depend on
///          scheduling again -- the exact defect this class exists to remove. The
///          pre-pass MUST be sequential over the slim file. It is a single scan
///          of a few columns and is not the bottleneck; the mixing pass is.
///   MIXED  (donors, pool_bin, report_underfilled, n_bins, n_filled): const, and
///          safe to call concurrently from any number of threads with no
///          synchronisation whatsoever.
///
/// freeze() is the barrier between them. Publish the pool to the worker threads
/// after it, and the C++ memory model does the rest.
///
/// ---------------------------------------------------------------------------
/// A RESIDUAL, STATED PLAINLY
/// ---------------------------------------------------------------------------
/// An event sampled INTO the pool will later be mixed against the pool, and so
/// can be paired with its own photons -- reintroducing, for that event alone, a
/// real pi0 correlation into the "background". This is not the old self-mixing
/// bug (which hit every event); it hits at most `donors_per_bin` photons' worth
/// of events per bin out of the whole file, and its weight in the mixed spectrum
/// is O(1 / N_events_in_bin). It is bounded, it is stated, and it has not been
/// quantified. If it ever needs removing, the fix is to carry the donor's event
/// id and skip it at mixing time -- cheap, but a change to the mixing loop, not
/// to this class.

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <random>
#include <vector>

#include "stageB_bin/PoolGrid.hpp"

namespace pi0 {

/// One pooled photon, as flat floats.
///
/// FLAT FLOATS ARE THE ONE THING WORTH KEEPING from the old code's PhotonCache.
/// The mixed inner loop is pure arithmetic over millions of donors and is
/// memory-bound, so the layout is the algorithm: 20 bytes per photon, three per
/// cache line, no indirection, no virtual, no ROOT object.
///
/// PRECISION IS NOT A CORNER CUT. float carries ~7 significant digits; the m_gg
/// spectrum these feed is histogrammed in 1.5 MeV bins over [0, 0.3] GeV
/// (cuts.json, mgg_histogram), so the representation error is some five orders of
/// magnitude below the bin width. Stage A's doubles are converted by the caller
/// at offer() time.
struct DonorPhoton {
    float px{}, py{}, pz{};

    /// Energy, GeV. Photons are MASSLESS here: e == |p|, exactly as Stage A and
    /// core/Types.hpp define them.
    float e{};

    /// 1 / |p|, GeV^-1. PRECOMPUTED, and this is the whole reason the struct is
    /// not just a 4-vector: the mixing loop needs the opening angle for every
    /// (current, donor) pair, so it needs cos(theta) = (p1 . p2) / (|p1| |p2|).
    /// The donor's |p| is loop-invariant, so its reciprocal is computed once here
    /// rather than a division per pair -- and a division is the one arithmetic op
    /// that would dominate this loop.
    ///
    /// Since photons are massless this is also 1/e. Both names are kept because
    /// the loop reads it as a momentum reciprocal, not as an energy one.
    float inv_p{};
};

/// A frozen pool of donor photons, keyed on (i_Q2, i_xB, i_mult).
///
/// Lifecycle, and it is strictly one-way:
///     DonorPool pool(cuts.mixing.pool_grid, seed, cuts.mixing.donors_per_bin);
///     for (each event in the slim file, SEQUENTIALLY) {
///         const int b = pool.pool_bin(q2, xb, passing.size());
///         if (b >= 0) pool.offer(b, passing);      // `passing`: e-gamma cut APPLIED
///     }
///     pool.freeze();
///     pool.report_underfilled(std::clog);
///     // ... now const, shared, lock-free, and mixing is a pure function of it.
///
/// THIS CLASS TAKES NO `Binning`, DELIBERATELY -- see the constructor.
///
/// Every method throws rather than misbehaves: offer() after freeze(), donors()
/// before it, and an out-of-range bin anywhere are all std::logic_error or
/// std::out_of_range, never a silently wrong pool.
class DonorPool {
   public:
    /// \param pool_grid       cuts.json's `mixing.pool_grid`. Validated here.
    /// \param seed            Already derived by the caller (cuts.json's
    ///                        `mixing.seed_mode` says how -- "file_hash"). This
    ///                        class samples from a seed; deciding what to hash is
    ///                        not its job.
    /// \param donors_per_bin  Reservoir depth in PHOTONS. Must be > 0.
    ///
    /// \throws std::invalid_argument if `pool_grid` is not a usable grid or
    ///         `donors_per_bin` is 0.
    ///
    /// WHY THERE IS NO `const Binning&` PARAMETER HERE.
    /// -----------------------------------------------
    /// It is tempting -- the caller has one, and it would let this constructor
    /// assert that the analysis binning's Grid A agrees with the pool grid's Q^2
    /// and x_B edges. That assertion would be WRONG, and wrong in the most
    /// expensive way available.
    ///
    /// The two grids are not the same grid and are not meant to become one. The
    /// pool grid's edges are hand-authored configuration product edges, shared
    /// with the N_DIS(Q^2, x_B) normalisation denominator and required to stay
    /// identical TO IT. Grid A's edges are EQUAL-STATISTICS, fitted to data by
    /// the make_grid tool, and will diverge from both the day it first runs on
    /// real data. They coincide only because Grid A's shipped file is still a
    /// placeholder copied out of cuts.json; that file's own _coincidence_warning
    /// says it outright -- code built on the coincidence "will be silently wrong
    /// on that day and correct on every day before it, which is the worst
    /// possible failure schedule".
    ///
    /// So an equality check here would pass every test today and throw on the
    /// first production grid. And a `Binning` parameter that is accepted but
    /// never read would be a parameter declared, set, and never used -- the exact
    /// defect this project is a monument to ending, reproduced inside the class
    /// written to end the mixing half of it.
    ///
    /// A donor is selected by its POOL bin, never by its analysis leaf. The pool
    /// needs no binning, so it does not take one.
    DonorPool(const PoolGridCuts& pool_grid, std::uint64_t seed, std::size_t donors_per_bin);

    /// Offer one event's photons to `pool_bin` by seeded reservoir sampling.
    ///
    /// *** `event_photons` MUST ALREADY HAVE THE e-gamma CUT APPLIED. *** See the
    /// file header. This is the contract; it is the caller's to keep.
    ///
    /// The reservoir's unit is the PHOTON, not the event: `donors_per_bin` photons
    /// are retained per bin, and donors() hands back exactly that flat list. Whole
    /// donor events are not kept because mixing never pairs a donor with a donor
    /// -- only current-event photons with pool photons -- so the grouping carries
    /// no information the mixing loop can use, and dropping it is what makes the
    /// inner loop a flat scan.
    ///
    /// Sampling is Algorithm R with a PER-BIN engine seeded from (seed, bin), so a
    /// bin's reservoir depends only on the sequence of photons offered to THAT
    /// bin -- not on how bins interleave.
    ///
    /// \param pool_bin       from pool_bin(). Must be in `[0, n_bins())`.
    /// \param event_photons  the event's e-gamma-passing photons.
    ///
    /// \throws std::logic_error   if the pool is already frozen.
    /// \throws std::out_of_range  if `pool_bin` is not a bin.
    /// \throws std::invalid_argument if `event_photons.size()` does not fall in
    ///         the multiplicity class `pool_bin` encodes. THIS IS THE TRIPWIRE
    ///         described in the file header: it means the count used to pick the
    ///         bin and the photons actually handed over disagree, which is what
    ///         happens when someone filters one and not the other. An empty vector
    ///         is in no class, so it throws here rather than passing as a no-op.
    void offer(int pool_bin, const std::vector<DonorPhoton>& event_photons);

    /// End the build phase. After this the pool is immutable and donors() works.
    /// Idempotent; releases the reservoirs' surplus capacity.
    void freeze();

    /// \return the frozen donor list for `pool_bin`. Possibly empty -- an empty
    ///         bin is a real and reportable state, not an error (see
    ///         report_underfilled).
    /// \throws std::logic_error   if the pool is not frozen. Reading a reservoir
    ///         mid-build hands out a list that is still moving, and "mixing is a
    ///         pure function of (event, pool)" is only true of a pool that has
    ///         stopped changing.
    /// \throws std::out_of_range  if `pool_bin` is not a bin. In particular
    ///         pool_bin() returns -1 for out-of-grid, and -1 is not a bin: check
    ///         it, do not pass it.
    [[nodiscard]] const std::vector<DonorPhoton>& donors(int pool_bin) const;

    /// \param n_photons  the number of photons in the event that PASSED the
    ///                   e-gamma cut -- the same count whose photons are handed to
    ///                   offer(), and the same definition on the mixing pass. That
    ///                   is the count governing how many pairs the event actually
    ///                   makes in the same-event spectrum, which is what
    ///                   "comparable combinatorics" has to mean. Use ONE definition
    ///                   everywhere or the multiplicity matching is decorative.
    /// \return the bin in `[0, n_bins())`, or -1 when Q^2 or x_B is off the grid,
    ///         or when `n_photons` is in no class -- which INCLUDES 0 photons.
    [[nodiscard]] int pool_bin(double q2_gev2, double xb, std::size_t n_photons) const;

    /// Write a QA report naming every bin that never reached `donors_per_bin`.
    ///
    /// KEPT FROM THE OLD CODE, WHICH GOT THIS RIGHT. Its finalisation step warned
    /// for every bin that never reached depth, and that is the single most useful
    /// thing it printed: it is where you look for leaves whose background estimate
    /// is thin. A thin bin does not fail -- it produces a noisy background shape
    /// that looks fine and subtracts badly.
    ///
    /// Reports the count retained AND the count offered, because they say
    /// different things: a bin that saw 50000 photons and holds 200 is healthy,
    /// while a bin that holds 200 because it only ever saw 200 is at the edge of
    /// being thin and the next file may push it over. Bins that saw NOTHING are
    /// listed separately and first -- an empty bin means no mixed background at
    /// all there, which is a different problem from a sparse one.
    ///
    /// \throws std::logic_error if the pool is not frozen. Occupancy is a
    ///         finalisation report; taken mid-build it describes a moving target.
    void report_underfilled(std::ostream& os) const;

    /// \return the number of pool bins. 224 for the shipped grid.
    [[nodiscard]] std::size_t n_bins() const { return m_donors.size(); }

    /// \return the number of bins holding at least one donor. Compare with
    ///         n_bins(): the shortfall is bins the file never populated.
    [[nodiscard]] std::size_t n_filled() const;

    /// \return the reservoir depth this pool was built to. Not a cut; a knob.
    [[nodiscard]] std::size_t donors_per_bin() const { return m_donors_per_bin; }

    /// \return true once freeze() has been called.
    [[nodiscard]] bool frozen() const { return m_frozen; }

    /// \return how many photons were OFFERED to `pool_bin` over the whole
    ///         pre-pass -- the reservoir's denominator, not its contents. This is
    ///         what makes an underfilled bin diagnosable rather than merely
    ///         reportable. Available before freeze().
    /// \throws std::out_of_range if `pool_bin` is not a bin.
    [[nodiscard]] std::uint64_t n_offered(int pool_bin) const;

   private:
    /// Decompose a bin into (i_q2, i_xb, i_mult). Row-major, and the numbering is
    /// frozen: bin = (i_q2 * n_xb + i_xb) * n_mult + i_mult.
    ///
    /// This mirrors core/Binning.hpp's convention -- slowest axis first, and the
    /// formula written down in the header rather than nowhere. The old analysis's
    /// leaf formula was recorded in no file and had to be reverse-engineered out
    /// of output histograms; a pool bin is no more self-describing than a leaf.
    void decompose(int pool_bin, std::size_t& i_q2, std::size_t& i_xb, std::size_t& i_mult) const;

    void require_bin(int pool_bin) const;

    PoolGridCuts m_grid;
    std::size_t m_donors_per_bin{};
    bool m_frozen{false};

    std::vector<std::vector<DonorPhoton>> m_donors;  ///< [bin] -> reservoir
    std::vector<std::uint64_t> m_seen;               ///< [bin] -> photons offered
    std::vector<std::mt19937_64> m_engines;          ///< [bin] -> that bin's stream
};

}  // namespace pi0
