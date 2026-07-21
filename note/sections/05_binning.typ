#import "../template/lib.typ": *

= Kinematic binning <sec:binning>

The four-dimensional binning is a factorized equal-statistics product grid.
Two grids are fitted to the data by `make_grid` and frozen into
`config/binning/`: Grid A over $(Q^2, x_B)$ and Grid B over $(z, p_T^2)$. The
4D bin is the product $A times B$, with Grid B independent of the Grid A
cell. Each axis is cut at data quantiles so that the marginal counts are
equal along it, and the edges are global, quotable and hashed into every
output (@sec:provenance-gaps).

This is a deliberate trade against an adaptive binning. An adaptive
tree can follow the strong correlations in the data --- $Q^2$ with $x_B$, and
$p_T^2$ with $z$ --- and equalise its leaves; a product grid cannot, and
leaves its corner cells under- or un-populated. What the product grid buys in
return is reproducibility: global rectangular edges that are identical on
every run, version-controlled, hashed into every result, and quotable in every
dimension --- none of which an unseeded, per-cell-adaptive tree offered.

== The grids <sec:grids>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, left, center, center),
    table.header([*Grid / tree*], [*Axes*], [*Shape*], [*Bins*]),
    [Grid A (per event)], [$Q^2, x_B$], [$8 times 7$], [$56$ ($47$ populated)],
    [Grid B (per $pi^0$)], [$z, p_T^2$], [$5 times 5$], [$25$],
    [4D --- multiplicity, @BSA], [$A times B$], [$56 times 25$], [$1400$],
    [3D --- $p_T$ broadening], [$A times (z "of" B)$], [$56 times 5$], [$280$],
  ),
  caption: [The grids. The same/mixed $m_(gamma gamma)$ spectra and the
  per-bin $phi_h$ histogram share the 4D grid; the broadening uses
  $A times (z "axis of" B)$, with $p_T^2$ as the observable rather than a
  binning dimension. Nine of the $56$ Grid A cells are the empty low-$x_B$ /
  high-$Q^2$ corner (@fig:kinematics), so $47$ carry data.],
) <tab:trees>

Bins are indexed by a single row-major integer,

$ "bin4d" = (i_(Q^2) dot 7 + i_(x_B)) dot 25 + (i_z dot 5 + i_(p_T^2)) , $ <eq:leaf>

with the axis order load-bearing --- `Binning::load` checks the axis names
against the order it expects and refuses a swapped pair. Grid A and Grid B are
computed from a single pooled scan of all four targets, so the edges are
common to numerator and denominator and the ratio of @eq:multiplicity is
formed bin-for-bin.

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    table.header([*Axis*], [*Edges*]),
    [$Q^2$ (GeV#super[2])],
      [$1.000, 1.103, 1.226, 1.378, 1.571, 1.831, 2.220, 2.927, 10.02$],
    [$x_B$],
      [$0.0625, 0.097, 0.120, 0.144, 0.174, 0.213, 0.276, 0.756$],
    [$z$], [$0.0475, 0.139, 0.193, 0.261, 0.375, 1.000$],
    [$p_T^2$ (GeV#super[2])], [$0.000, 0.025, 0.059, 0.110, 0.210, 4.995$],
  ),
  caption: [The frozen edges, computed by `make_grid` over the $95$ pooled
  slims ($5.49 times 10^6$ events, $2.51 times 10^6$ $pi^0$;
  `provenance_hash 2acf618294a6c3b0`). Each axis is quantile-spaced so its
  marginal is flat: the Grid A row sums are equal to $approx 686 000$ per
  $Q^2$ row and $approx 784 000$ per $x_B$ column, the Grid B marginals to
  $approx 502 000$ each, to within rounding. The cells are deliberately not
  flat --- the two axes are correlated, and a factorized grid cannot follow
  that.],
) <tab:kdtree-geometry>

#wide-figure(
  "../figures/binning_grid_factorized.pdf",
  [The factorized equal-statistics grid, drawn over the pooled data density it
  partitions (Grid A per event, Grid B per $pi^0$; log colour scale). The grid
  lines are global and rectangular in every dimension --- the reproducibility
  this scheme is built for --- and the density makes the trade visible: the
  edges bunch in tight where the data piles up, while the $Q^2$--$x_B$
  correlation leaves the top-left and bottom-right corner cells empty, the
  adaptivity a product grid gives up in exchange.],
  <fig:binning-factorized>,
  width: 88%,
)

== The top-bin problem <sec:binning-caveat>

#important-box(title: "The widest bins must be plotted at their mean, not their midpoint")[
  An equal-statistics grid packs equal counts into each bin, so the
  outermost bin of each axis --- which runs to the kinematic limit --- is
  enormously wide while holding no more $pi^0$ than its narrow neighbours. The
  top $z$ bin spans $[0.375, 1.0]$ and the top $p_T^2$ bin $[0.21, 4.99]$
  GeV#super[2], against interior widths of a few hundredths.

  Its geometric midpoint is therefore not where the data sits. The top $z$
  box has a true, sideband-subtracted $chevron.l z chevron.r approx 0.53$, not
  its midpoint $0.69$; the top $p_T^2$ box has
  $chevron.l p_T^2 chevron.r approx 0.3$ GeV#super[2], not $approx 2.6$.
  Plotting a bin at its midpoint would understate the depth of the attenuation
  at high $z$ and manufacture an apparent high-$p_T$ (Cronin-like) rise that is
  purely an abscissa artefact.

  The values in a bin are correct integrals over it; only a midpoint
  abscissa would be wrong --- which is why this analysis never uses one.
]

=== How it is handled <sec:binning-caveat-fix>

Both sides of the chain carry the true bin positions, so the caveat above is a
description of the machinery, not an outstanding bug:

/ C++: Stage B accumulates the per-bin kinematic sums `sum_q2`, `sum_xb`,
  `sum_z`, `sum_pt2` in the `spectra` tree and `sum_pt2`, `sum_pt4` in
  `ptb3d`, filled over exactly the $pi^0$ entering each bin, with a
  compensated (Neumaier) accumulator so the sums are bit-identical across
  thread counts.

/ Python: every abscissa resolves to the sideband-subtracted mean
  $chevron.l X chevron.r = "sum"_X \/ n_"same"$ over the signal mass region,
  never the box centre. This is why the top $z$ point in @tab:R-vs-z and
  @fig:R-vs-z sits at $chevron.l z chevron.r = 0.53$, and why the broadening
  of @sec:results-ptb is a moment difference rather than a box-labelled point.

== Auxiliary grids <sec:aux-grids>

Two fixed grids in `config/cuts.json` accompany the fitted ones and should not
be confused with them.

/ Normalisation grid: the @DIS denominator $N^"DIS"(Q^2, x_B)$ of
  @eq:multiplicity is counted once per event into the `n_dis` tree, on Grid A
  itself --- the same $(Q^2, x_B)$ cells the numerator uses --- so no
  area-fraction reweighting is needed between them (@sec:normalisation).

/ Mixing pool grid: the mixed-event donor pool is keyed on a hand-maintained
  $(Q^2, x_B)$ product grid (`/mixing/pool_grid`) plus a photon-multiplicity
  class (@sec:mixing). It is kept identical to the $N^"DIS"$ denominator
  binning by construction, and is deliberately separate from the fitted Grid A
  so that a re-fit of the physics grid cannot silently move the background
  model.
