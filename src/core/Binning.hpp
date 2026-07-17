#pragma once

/// \file Binning.hpp
/// \brief The frozen, factorized kinematic binning: two rectangular grids whose
///        edges are loaded from version-controlled JSON and hashed into every
///        output. No ROOT, no HIPO.
///
/// WHY THIS TYPE EXISTS, AND WHY IT IS NOT A KD-TREE
/// ------------------------------------------------
/// The superseded analysis binned in 4D with an adaptive kd-tree grown from an
/// UNSEEDED, thread-timing-dependent reservoir sample. Two passes over the same
/// data produced DIFFERENT edges, and the serialised tree was never archived
/// with the results -- so the binning of the production is unrecoverable, and
/// the published geometry had to be REVERSE-ENGINEERED out of the output files
/// by chaining box centres (note sec:binning-reconstruction). Even the leaf
/// index formula was written down nowhere and had to be inferred.
///
/// This project replaces that with FACTORIZED RECTANGULAR EQUAL-STATISTICS
/// GRIDS:
///
///   Grid A = (Q^2, x_B)   default 8 x 7 = 56 cells
///   Grid B = (z,   p_T^2) default 5 x 5 = 25 cells, INDEPENDENT of the A cell
///   4D bin = A x B                          = 1400
///   p_T broadening = A x (z axis of B)      = 8 x 7 x 5 = 280
///
/// The edges are computed ONCE by the `make_grid` tool, written to
/// config/binning/*.json, committed, and hashed by provenance_hash() into every
/// output file. Two runs of this code against the same JSON produce the same
/// bins, by construction: nothing here samples, and nothing here depends on
/// thread scheduling.
///
/// What this deliberately sacrifices is stated plainly in the note
/// (sec:binning-future): because Grid B does not adapt to the Grid A cell, the
/// z and p_T^2 edges can no longer follow the strong p_T^2-z correlation the
/// nested tree captured. That is a known, accepted trade for reproducibility.
///
/// WHAT THIS TYPE DOES NOT SOLVE
/// -----------------------------
/// An equal-statistics rectangular grid has EXACTLY the same unbounded
/// outermost bin as the kd-tree did, and therefore exactly the same bin-centre
/// disease (note sec:binning-caveat): 52.6% of the old bins sat in the top box
/// of at least one dimension, 15.5% of them implied y > 0.85 -- violating the
/// DIS cut their own events passed -- and 6.9% implied nu > E_beam, which is
/// kinematically impossible. Reporting a bin at (lo + hi) / 2 is wrong here for
/// the same reason it was wrong there.
///
/// So: THIS HEADER DELIBERATELY OFFERS NO bin_centre() FUNCTION, and no
/// accessor that hands back a box midpoint. There is nothing to call, because
/// there is nothing anyone should call. The abscissa of a bin is the
/// COUNT-WEIGHTED MEAN over the events under the fitted peak, accumulated per
/// (4D bin, m_gg bin) by the Stage B accumulators. The grid's job is to answer
/// "which bin is this event in"; it is not the authority on where that bin is.
///
/// Units: Q^2 and p_T^2 in GeV^2, x_B and z dimensionless. No angles here, so
/// the project's degrees-at-every-boundary rule has nothing to bite on.

#include <string>
#include <vector>

namespace pi0 {

/// One axis of a grid: a name and a monotonically increasing edge array.
///
/// `edges.size() - 1` bins. A default-constructed Grid1D has no edges, zero
/// bins, and find() returns -1 for every value -- it is inert rather than
/// undefined, but it is not a usable grid. Only Binning::load() produces one
/// that is, and it validates.
struct Grid1D {
    /// Axis name, e.g. "q2". Carried into error messages and into
    /// provenance_hash(), so a mislabelled axis is a different hash.
    std::string name;

    /// Bin edges, strictly increasing, at least 2 of them. Validated by
    /// Binning::load(); NOT validated by this struct's own methods, which are
    /// hot-path code called once per pi0.
    std::vector<double> edges;

    /// Number of bins = edges.size() - 1, or 0 if there are fewer than 2 edges.
    [[nodiscard]] int nbins() const;

    /// Bin index of `v` in [0, nbins()), or -1 if `v` is outside the axis.
    ///
    /// SEMANTICS, stated exactly because the old code got the last one wrong:
    ///   * bins are HALF-OPEN [lo, hi): a value on a bin's lower edge belongs
    ///     to that bin, a value on its upper edge belongs to the NEXT bin;
    ///   * a value exactly on the TOP edge, edges.back(), lands in the LAST
    ///     bin. It does NOT return -1. This is the one place the half-open rule
    ///     is broken, and it is broken deliberately: the top edge is the
    ///     kinematic limit, so the alternative is silently discarding the
    ///     events that sit on it. The old find_1d_bin() returned -1 there and
    ///     so dropped every p_T^2 exactly equal to the top edge; tests/
    ///     test_binning.cpp pins this behaviour so it cannot come back.
    ///   * `v` below edges.front() -> -1; above edges.back() -> -1;
    ///   * NaN -> -1 (the comparisons below are written so that NaN fails them
    ///     rather than reaching upper_bound with a value that has no ordering).
    [[nodiscard]] int find(double v) const;
};

/// Two axes and a flat, row-major cell index over them.
struct Grid2D {
    Grid1D x;  ///< the SLOW axis: it multiplies y.nbins() in the flat index
    Grid1D y;  ///< the FAST axis: it is the flat index's remainder

    /// Flat, ROW-MAJOR cell index:
    ///
    ///     cell = ix * y.nbins() + iy
    ///
    /// with ix = x.find(vx) and iy = y.find(vy). Returns -1 if EITHER axis
    /// returns -1 -- there is no partial hit and no flow bin. A caller that
    /// wants to know which axis missed must ask the axes.
    [[nodiscard]] int find(double vx, double vy) const;

    /// x.nbins() * y.nbins().
    [[nodiscard]] int ncells() const;
};

/// The whole binning: Grid A over (Q^2, x_B), Grid B over (z, p_T^2).
///
/// THE FLAT INDEX FORMULAE. Read this before writing any decoder.
/// -------------------------------------------------------------
/// These are written down here, in the header, because the superseded
/// analysis's leaf formula was written down NOWHERE and had to be
/// reverse-engineered from output files. The downstream Python decodes these
/// indices; this comment is the contract it implements against.
///
/// Let
///     n_q2  = A.x.nbins()      n_xb  = A.y.nbins()
///     n_z   = B.x.nbins()      n_pt2 = B.y.nbins()
/// and let i_q2, i_xb, i_z, i_pt2 be the per-axis indices.
///
/// Grid A cell (row-major, Q^2 slow):
///     a = i_q2 * n_xb + i_xb                        in [0, n_q2 * n_xb)
/// Grid B cell (row-major, z slow):
///     b = i_z  * n_pt2 + i_pt2                      in [0, n_z * n_pt2)
///
/// The 4D bin is the A cell times the B cell, with the B CELL FAST:
///
///     bin4d = a * B.ncells() + b
///           = ((i_q2 * n_xb) + i_xb) * (n_z * n_pt2) + (i_z * n_pt2) + i_pt2
///
/// which for the defaults (8, 7, 5, 5) runs over [0, 1400). Note this is the
/// same lexicographic Q^2 -> x_B -> z -> p_T^2 ordering the old leaf formula
/// used, so the two are directly comparable -- the bins differ, the index
/// convention does not.
///
/// TO DECODE bin4d (integer division throughout):
///     b     = bin4d % (n_z * n_pt2)          a     = bin4d / (n_z * n_pt2)
///     i_pt2 = b     % n_pt2                  i_z   = b     / n_pt2
///     i_xb  = a     % n_xb                   i_q2  = a     / n_xb
///
/// The 3D p_T-broadening bin is the A cell times the z AXIS of B -- p_T^2 is
/// the observable there, so it is not a binning dimension:
///
///     bin3d = a * n_z + i_z
///           = ((i_q2 * n_xb) + i_xb) * n_z + i_z    in [0, 280) by default
///
/// TO DECODE bin3d:
///     i_z   = bin3d % n_z                    a     = bin3d / n_z
///     i_xb  = a     % n_xb                   i_q2  = a     / n_xb
///
/// bin3d is NOT bin4d / n_pt2. It happens to equal it whenever both are in
/// range, but bin4d == -1 carries no i_z, so do not derive one from the other.
struct Binning {
    Grid2D A;  ///< (Q^2, x_B). x = Q^2 (GeV^2), y = x_B.
    Grid2D B;  ///< (z, p_T^2). x = z, y = p_T^2 (GeV^2). Independent of A.

    /// Flat 4D bin index, or -1 if either grid misses. See the formula above.
    [[nodiscard]] int find_4d(double q2, double xb, double z, double pt2) const;

    /// A.ncells() * B.ncells(). 1400 for the default grids.
    [[nodiscard]] int n4d() const;

    /// Flat 3D bin index for p_T broadening (A cell x z bin), or -1 if the A
    /// cell or the z bin misses. p_T^2 is deliberately not an argument.
    [[nodiscard]] int find_3d(double q2, double xb, double z) const;

    /// A.ncells() * B.x.nbins(). 280 for the default grids.
    [[nodiscard]] int n3d() const;

    /// A 16-hex-digit fingerprint of the grid GEOMETRY, to stamp into every
    /// output so a result can always be matched back to the binning that made
    /// it. That is the provenance hole this whole rewrite exists to close: the
    /// old production's edges are gone, and no output carries anything that
    /// could identify them.
    ///
    /// It covers the axis names and every edge, at full double precision, in
    /// order. It covers NOTHING ELSE: the JSON's _comment and provenance blocks
    /// are documentation and do not enter, so re-wording a comment does not
    /// invalidate a result, while moving an edge by one ulp does.
    ///
    /// It is FNV-1a over a canonical decimal serialisation, not a cryptographic
    /// digest. It answers "are these the same edges", which is a
    /// collision-by-accident question, not "did someone forge these edges".
    ///
    /// Stable across runs, builds and machines: the serialisation is fixed
    /// (17 significant digits, C locale) and the arithmetic is exact 64-bit.
    [[nodiscard]] std::string provenance_hash() const;

    /// Load Grid A from `grid_a_json` and Grid B from `grid_b_json`.
    ///
    /// TWO FILES, NOT ONE, AND NOT config/cuts.json. The grids are a build
    /// product of the `make_grid` tool run over real data -- machine-written,
    /// regenerated per dataset, and carrying their own provenance and their own
    /// hash. cuts.json is hand-maintained physics that a human argues about.
    /// Mixing the two would mean a tool rewrites the file that holds the cut
    /// values, and would put the edges behind Cuts::load()'s
    /// checked-against-Constants.hpp contract, which does not apply to them.
    /// Grid A and Grid B are separate files for the same reason they are
    /// separate grids: they are independent, and the note's `make_grid` scans
    /// them independently.
    ///
    /// cuts.json's /mixing/pool_grid keeps its OWN Q^2 and x_B edges, and that
    /// is not an oversight to be tidied away. The pool grid and the DIS
    /// normalisation denominator N_DIS(Q^2, x_B) are histogrammed on the
    /// configuration PRODUCT edges and must stay identical to each other (see
    /// that block's _comment). Grid A is an EQUAL-STATISTICS grid fitted to
    /// data and will not equal them once make_grid has run. Today the
    /// placeholder files happen to hold the same numbers; that coincidence is
    /// temporary and nothing may be built on it.
    ///
    /// FAILS LOUDLY, with no defaulting and no partial success, on: a file that
    /// cannot be opened; a document that is not valid JSON or not an object; a
    /// missing key; a key of the wrong JSON type; an axes array that is not
    /// exactly two axes; fewer than 2 edges on an axis; edges that are not
    /// strictly increasing; a non-finite edge.
    ///
    /// The monotonicity check is not ceremony. Non-monotonic edges do not
    /// crash: upper_bound over an unsorted range is merely meaningless, so the
    /// events would be silently filed into the wrong bins and every downstream
    /// number would be quietly wrong. That failure has to be caught at load or
    /// it is not caught at all.
    ///
    /// \throws std::runtime_error naming the file and the full key path.
    [[nodiscard]] static Binning load(const std::string& grid_a_json,
                                      const std::string& grid_b_json);
};

}  // namespace pi0
