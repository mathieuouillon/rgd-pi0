#pragma once

/// \file PoolGrid.hpp
/// \brief The runtime image of cuts.json's `mixing` block: the event-mixing
///        donor-pool grid (Q^2, x_B, n_photons class) and its two knobs.
///
/// Defined HERE, next to the code that consumes it, rather than in
/// config/Cuts.hpp -- the same arrangement as `PairingCuts` in core/Pairing.hpp
/// and `VertexCuts` in vertex/VzCorrector.hpp. `Cuts` composes this block; this
/// header knows nothing about `Cuts`, so the mixing code is testable without
/// ever parsing a JSON file.
///
/// Cuts::load() is the only thing that fills a MixingCuts from the config, and
/// it fails loudly on a missing key like every other block. Until it did, the
/// entire `mixing` block sat in cuts.json DECLARED AND NEVER READ -- precisely
/// the defect (see that file's header) this project exists to end, sitting
/// inside the file written to end it.
///
/// THE POOL GRID IS NOT GRID A. READ THIS BEFORE UNIFYING ANYTHING.
/// ---------------------------------------------------------------
/// The Q^2 and x_B axes here are the CONFIGURATION PRODUCT EDGES: hand-authored
/// round numbers, shared with -- and required to stay identical to -- the DIS
/// normalisation denominator N_DIS(Q^2, x_B). They are NOT the analysis
/// binning's Grid A (src/core/Binning.hpp), which is EQUAL-STATISTICS edges
/// fitted to data by the make_grid tool.
///
/// Today the two happen to hold the same numbers, because Grid A's shipped file
/// is still a placeholder copied from this block. That coincidence is temporary
/// and NOTHING MAY BE BUILT ON IT -- see config/binning/grid_A_q2_xb.json's own
/// _coincidence_warning. Code that lets one grid stand in for the other is
/// correct every day until make_grid runs on real data and silently wrong every
/// day after, which is the worst failure schedule available. Keep them apart.
///
/// TRAP, inherited: the old [Pi0Analysis] `W_edges` and `nu_edges` keys were
/// declared and registered but NEVER READ -- the pool key dropped W and nu while
/// the config comment still described a five-dimensional key. There is no W and
/// no nu here. Do not add them back believing you are restoring something.
///
/// Units: GeV^2 for Q^2, dimensionless for x_B.

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "core/Binning.hpp"

namespace pi0 {

/// One photon-multiplicity class of the pool key, from
/// `mixing.pool_grid.n_photons_classes`. The shipped set is {1, 2, 3, >=4}.
///
/// Matching donors on multiplicity is what makes the mixed combinatorics
/// comparable to the same-event combinatorics: an event with two photons and an
/// event with nine do not produce backgrounds of the same shape.
struct MultClass {
    std::string label{};  ///< the config's own label, e.g. "1" or ">=4". For messages only.
    int min{};            ///< inclusive lower bound on the photon count.

    /// Inclusive upper bound, or EMPTY for an open-ended class (the config's
    /// `"max": null`).
    ///
    /// std::optional rather than a sentinel, deliberately. cuts.json's
    /// `extraction.fit_bounds` block already warns that a null must map to the
    /// platform infinity and "NOT to zero and NOT to a large finite sentinel";
    /// an int has no infinity, so the only way to keep that promise is to make
    /// "no upper bound" a state the type can represent and a reader cannot
    /// misread as a number.
    std::optional<int> max{};

    /// \return true when a `n_photons`-photon event belongs to this class.
    [[nodiscard]] bool contains(std::size_t n_photons) const;
};

/// The pool key's grid: (i_Q2, i_xB, i_mult), 8 x 7 x 4 = 224 bins as shipped.
///
/// The axes are `Grid1D` -- the project's one axis type, from core/Binning.hpp.
/// Reusing it rather than inventing a second edge-vector-with-a-lookup means the
/// pool key and the analysis binning share their bin-boundary semantics exactly:
/// half-open [lo, hi), the top edge lands in the LAST bin rather than nowhere,
/// out-of-range and NaN give -1, no flow bins. Two axis types would be two
/// chances to disagree about which bin an edge value belongs to, and that
/// disagreement would be invisible in every histogram it corrupted.
///
/// Sharing the TYPE is not sharing the EDGES. See this file's header.
struct PoolGridCuts {
    Grid1D q2{};                                 ///< GeV^2. 8 bins as shipped.
    Grid1D xb{};                                 ///< dimensionless. 7 bins as shipped.
    std::vector<MultClass> n_photons_classes{};  ///< 4 classes as shipped.

    /// \return the class index in `[0, n_photons_classes.size())`, or -1 when
    ///         `n_photons` belongs to no class -- which INCLUDES zero photons.
    ///         A zero-photon event produces no pairs and enters no class.
    [[nodiscard]] int mult_class(std::size_t n_photons) const;

    /// \return q2.nbins() * xb.nbins() * n_photons_classes.size(). 224 as shipped.
    [[nodiscard]] std::size_t n_bins() const;

    /// \throws std::invalid_argument if this grid could not key a pool.
    ///
    /// Checks the EDGES too -- fewer than 2, non-finite, or not strictly
    /// increasing -- and not because it is tidy. `Grid1D` is a plain aggregate
    /// whose own header says it is "NOT validated by this struct's own methods,
    /// which are hot-path code": Binning::load() validates the grids IT loads,
    /// and it never sees this one. Non-monotonic edges do not crash; upper_bound
    /// over an unsorted range is merely meaningless, so photons would be filed
    /// into the wrong pool bins and every mixed spectrum would be quietly wrong.
    /// That has to be caught here or it is not caught at all.
    ///
    /// Called by Cuts::load() and by DonorPool's constructor. Both, not either:
    /// load() catches a bad config, and the constructor catches a PoolGridCuts
    /// assembled in code that never went through load().
    void validate() const;
};

/// The whole `mixing` block. Composed into `Cuts` as `Cuts::mixing`.
struct MixingCuts {
    /// Donor PHOTONS retained per pool bin -- not donor events. See DonorPool.hpp
    /// on why the reservoir's unit is the photon.
    ///
    /// A DESIGN KNOB, and a starting value rather than a measurement. It is NOT
    /// the old rolling FIFO's depth of 50 and is not comparable to it: that 50
    /// counted events in a pool that was cold early in a file and warm late, so
    /// the old mixed statistics were an emergent property of warm-up. Here the
    /// depth is uniform across the file by construction, so it is a number you
    /// set.
    std::size_t donors_per_bin{};

    /// How the reservoir seed is derived. "file_hash" is the only implemented
    /// mode: the seed comes from the input file's identity, so a given slim file
    /// always yields the same pool regardless of thread scheduling, and
    /// (file, seed, donors_per_bin) fully determines every mixed pair.
    ///
    /// DonorPool does not read this. It takes an already-derived
    /// `std::uint64_t seed`, because deciding what to hash is the caller's job
    /// and sampling from a seed is the pool's. Cuts::load() rejects any other
    /// mode rather than let a config name a scheme nothing implements.
    std::string seed_mode{};

    PoolGridCuts pool_grid{};
};

}  // namespace pi0
