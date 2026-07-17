#import "../template/lib.typ": *

= Kinematic binning <sec:binning>

The four-dimensional binning is *adaptive*: rather than a fixed rectangular
grid, a kd-tree is grown from the data so that every leaf holds
approximately the same number of $pi^0$. This is what makes a genuinely
four-dimensional measurement possible --- a product grid fine enough to
resolve $z$ and $p_T^2$ at high $Q^2$ would be mostly empty.

== The stratified kd-tree <sec:kdtree>

The tree is built by `AdaptiveBinnerND<4>` in *stratified* mode. The
dimensions are split *in a fixed order*,

$ Q^2 arrow.r x_B arrow.r z arrow.r p_T^2 , $

and at each level the parent cell is cut into $k$ *equal-count* groups at
data quantiles. A leaf is emitted once all four dimensions have been
split. There is no minimum-occupancy requirement: the only guards are
$k >= 1$ per dimension and that the sample is at least as large as the
requested bin count.

#note-box(title: "Only $Q^2$ has global edges")[
  Because $Q^2$ is split first, at the root, its seven split values are
  common to the whole tree and may be quoted as global edges. *Every other
  dimension's edges are local to its parent cell*: the $x_B$ edges differ
  per $Q^2$ slice, the $z$ edges per $(Q^2, x_B)$ cell, and the $p_T^2$
  edges per $(Q^2, x_B, z)$ cell.

  Quoting a single global $x_B$, $z$ or $p_T^2$ edge array would therefore
  be *wrong*. This also explains why the analysis code groups results by
  each leaf's exact box rather than by a display grid. It is confirmed
  directly in the output: the $p_T^2$ binning has $1450$ distinct centres
  against only $290$ distinct $z$ centres --- one $p_T^2$ set per
  $(Q^2,x_B,z)$ group, exactly as the nesting implies.
]

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, left, center, center),
    table.header([*Tree*], [*Dimensions*], [*Splits*], [*Leaves*]),
    [Multiplicity], [$Q^2, x_B, z, p_T^2$], [$8 times 8 times 5 times 5$],
      [$1600$ ($1450$ filled)],
    [$p_T$ broadening], [$Q^2, x_B, z$], [$8 times 8 times 5$], [$320$],
    [@BSA], [$Q^2, x_B, z, p_T^2$], [free split],
      [$250$ ($237$ used)],
  ),
  caption: [The three kd-trees. The multiplicity tree is shared with the
  event-mixing histograms, so signal and background are binned identically.
  The $150$ missing multiplicity leaves are the kinematically empty
  low-$x_B$ corner: $58$ of the $64$ $(Q^2, x_B)$ cells are populated, and
  $58 times 5 times 5 = 1450$.],
) <tab:trees>

Leaves are indexed in lexicographic order,

$ "leaf" = ((i_(Q^2) dot 8 + i_(x_B)) dot 5 + i_z) dot 5 + i_(p_T^2) , $ <eq:leaf>

which is verifiable against every row of the output files. The tree is
built from a reservoir sample capped at $10^6$ $pi^0$, giving a target
occupancy of $approx 625$ per leaf.

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    table.header([*Quantity*], [*Value*]),
    [Global $Q^2$ edges (GeV#super[2])],
      [$1.0000$, $1.1001$, $1.2193$, $1.3652$, $1.5502$, $1.7981$,
       $2.1644$, $2.8245$, $9.5805$],
    [Global $x_B$ range], [$[0.0625, 0.7414]$],
    [Global $z$ range], [$[0.0494, 0.9999]$],
    [Global $p_T^2$ range (GeV#super[2])], [$[0, 3.6395]$],
  ),
  caption: [Global kd-tree geometry, reconstructed from the output files
  (@sec:binning-reconstruction). Only the $Q^2$ edges are true global bin
  boundaries; the other three rows are the sample's bounding box, which
  clamps the outermost leaf boxes.],
) <tab:kdtree-geometry>

#subfig2(
  (
    ("../figures/adaptive_binning_xB_Q2.pdf", [$(x_B, Q^2)$ plane],
      <fig:kdtree-a>),
    ("../figures/adaptive_binning_z_pT2.pdf", [$(z, p_T^2)$ plane],
      <fig:kdtree-b>),
  ),
  [The kd-tree rendered as exact two-dimensional tilings, sliced at
  count-weighted quantiles of the two dimensions not shown. The
  equal-occupancy construction is visible directly: cells are small where
  the data is dense and grow toward the edges of the acceptance. The wide
  outer cells in @fig:kdtree-b are the top $z$ and $p_T^2$ boxes discussed
  in @sec:binning-caveat --- they hold the same number of $pi^0$ as their
  narrow neighbours, but span a far larger kinematic range.],
  <fig:kdtree>,
)

== Reconstructing the binning <sec:binning-reconstruction>

The leaf-box tables (`pi0_adaptive_binning*.txt.boxes.csv`) are written by
the C++ pass onto `/work` and are *not archived with the results*
(@sec:provenance-gaps). The geometry in @tab:kdtree-geometry was therefore
recovered from the output files themselves, by chaining box centres: given
a known low edge $e_0$ and the successive centres $c_i$ of a cell's bins,

$ e_(i+1) = 2 c_i - e_i . $ <eq:chain>

Applied to all $58$ independent $x_B$ chains and all $290$ independent
$p_T^2$ chains, every chain closes on the same global top edge --- a
decisive validation that the reconstruction is correct and that
@eq:leaf is the right index decoding.

#warning-box(title: "Reproducibility: the tree is not deterministic")[
  The kd-tree is built from an *unseeded, thread-timing-dependent*
  reservoir sample. Two scan passes over the same data give *different
  trees*. The serialised `.txt` tree (17 significant digits) is the only
  reproducible record of a given binning; the companion `boxes.csv` is
  written at the iostream default of *6* significant digits and is
  display-only.

  Since neither file is archived with the results, *the exact binning of
  the current production is recoverable only from `/work`*. If those files
  are lost, the results cannot be re-binned or reproduced --- only
  re-derived approximately, as in @eq:chain. Archiving the `.txt` trees
  alongside the CSVs is the single cheapest provenance fix available.
]

== The bin-centre problem <sec:binning-caveat>

#warning-box(title: "Reported bin positions are not where the data is")[
  Every plot abscissa and every kinematic column in the result files is the
  *geometric midpoint* of the leaf box, $(e_"lo" + e_"hi") \/ 2$. For
  interior boxes this is harmless. For the *outermost* box in each
  dimension it is badly wrong, because that box extends to the global
  kinematic limit and is enormously wide.

  #figure(
    table(
      columns: (auto, auto, auto, auto),
      align: (left, center, center, center),
      table.header([*Dimension*], [*Top box*], [*Plotted centre*],
                   [*Typical interior width*]),
      [$z$], [$[0.369, 0.999]$], [$0.684$], [$0.05 - 0.10$],
      [$p_T^2$ (GeV#super[2])], [$[0.18, 3.64]$], [$approx 1.92$],
        [$0.03 - 0.07$],
      [$Q^2$ (GeV#super[2])], [$[2.824, 9.581]$], [$6.203$], [$0.1 - 0.7$],
      [$x_B$], [$[0.465, 0.741]$], [$0.603$], [$0.02 - 0.04$],
    ),
    caption: [The top box in each dimension, for a representative cell. The
    $z$ example is cell $(i_(Q^2), i_(x_B)) = (3,3)$, whose reconstructed
    edges are $0.049$, $0.148$, $0.199$, $0.261$, $0.369$, $0.999$ --- the
    last box is $6-12 times$ wider than its neighbours.],
  ) <tab:top-boxes>

  *$52.6%$ of all leaves sit in the top box of at least one dimension.*

  How badly wrong the centres are is best shown by a quantity that does not
  appear in the binning at all. The energy transfer follows from the two
  binned variables, $nu = Q^2 \/ (2 M_p x_B)$, so every leaf's centre implies
  a $nu$ --- and that implied $nu$ can be tested against the cuts the events
  themselves had to pass ($y < 0.85$, i.e. $nu < 8.95$ GeV; $W > 2$ GeV).
  Evaluating it leaf by leaf:

  #figure(
    table(
      columns: (1fr, auto),
      align: (left, right),
      table.header([*Box-centre kinematics of the 1450 leaves*], [*Leaves*]),
      [imply $y > 0.85$ --- violating the @DIS cut applied to their events],
        [225 ($15.5%$)],
      [imply $W < 2$ GeV --- violating the $W$ cut applied to their events],
        [175 ($12.1%$)],
      [imply $nu > E_"beam"$ --- *kinematically impossible*], [100 ($6.9%$)],
    ),
    caption: [Box centres tested against the selection their own events
    passed. The largest implied energy transfer is $nu = 21.7$ GeV, against a
    beam energy of $10.53$ GeV.],
  ) <tab:centre-unphysical>

  So the box centres are not merely *mis-placed*: for a hundred leaves they
  are points that *no event in the bin could occupy*, labelled with an energy
  transfer exceeding the beam energy. Every one of those leaves nonetheless
  contains real $pi^0$ obeying $y < 0.85$.

  The same test on the *count-weighted* means --- available today only in the
  @BSA output (@tab:schemas), which already computes them --- passes
  cleanly: $nu in [3.2, 7.1]$ GeV, $y in [0.25, 0.72]$, not one leaf outside
  the cuts. Identical data, identical bins; the abscissa is the only
  difference.

  Two consequences matter for reading this note:

  + *The "high-$z$" point is not at $z approx 0.70$.* It is an integral
    over $z in [0.37, 1.0]$, dominated by $z approx 0.4 - 0.5$ where the
    statistics are. Since the attenuation deepens with $z$, labelling that
    integral $z approx 0.70$ makes the suppression look far milder than it
    is at the $z$ the label claims --- so the *depth* of the attenuation in
    @tab:R-vs-z is understated by an amount that is itself unknown until the
    abscissae are fixed. It is a binning artefact, *not* a physics result.

  + *The apparent Cronin rise at $p_T^2 approx 1.9$ GeV#super[2] is an
    artefact.* The $R > 1$ value there is a genuine integral over
    $p_T^2 in [0.18, 3.64]$ holding $approx 20%$ of the $pi^0$, but the true
    $chevron.l p_T^2 chevron.r$ in that box is $approx 0.3$ GeV#super[2], not
    $1.9$. It must not be read as a high-$p_T$ measurement.

  The $R$ values themselves are *correct integrals over their boxes*; it is
  only their reported kinematic *position* that is wrong.
]

=== Status: fixed in code, pending a re-run <sec:binning-caveat-fix>

The fix has been implemented, on both sides of the chain:

/ C++: the analysis pass now accumulates the per-bin kinematic sums
  `4D/sum_Q2`, `4D/sum_xB`, `4D/sum_z` and `4D/sum_pT2` alongside the
  existing `4D/counts`, filled over exactly the $pi^0$ that enter each bin's
  count. Dividing gives $chevron.l X chevron.r$ per bin. This mirrors the
  `pT3D/sum_pT2` accumulators already used for the $p_T$ broadening
  (@sec:ptb-accumulators) and the count-weighted means the @BSA code already
  writes (`q2_mean`, `xb_mean`, `z_mean`, `pt2_mean`).

/ Python: every abscissa now resolves through a single helper that prefers
  the count-weighted mean and falls back to the box centre only when the
  sums are absent.

#important-box(title: "The results in this note still use box centres")[
  The fix changes nothing about @sec:results *yet*. The production files
  predate the accumulators, so no mean is available for them and the
  extraction falls back to box centres --- as it must, since the means
  cannot be reconstructed after the fact from the histograms alone.

  Two things follow. First, everything in @tab:R-vs-z and @fig:R-vs-z is
  still reported at a geometric box centre, and @sec:binning-caveat applies
  to it in full. Second, the fallback is no longer silent: the extraction
  prints a warning naming the problem, and every affected plot draws each
  bin's *full width* as an $x$ error bar, so a wide outer bin can no longer
  be mistaken for a point measurement.

  Making the results correct needs the farm pass re-run. Note the fix was
  implemented against the *factorized* C++ (@sec:binning-future), not the
  kd-tree that produced these results --- an equal-statistics rectangular
  grid has exactly the same unbounded outermost bin, so the same fix is
  required there and the next production gets it for free.
]

== Auxiliary grids <sec:aux-grids>

Two fixed rectangular grids coexist with the kd-tree and should not be
confused with it.

/ Normalisation grid: The @DIS denominator $N^"DIS"(Q^2, x_B)$ is
  histogrammed on the *configuration* product edges
  $Q^2 = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 11.0}$ and
  $x_B = {0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7}$ ($8 times 7 = 56$
  cells), *not* on the kd-tree quantiles. The Python must therefore
  integrate this histogram over each leaf's $(Q^2, x_B)$ box, which it does
  with an area-fraction weighting (@sec:normalisation).

/ Mixing pool grid: The event-mixing pool is keyed on the same
  configuration $Q^2$ and $x_B$ edges plus a photon-multiplicity class
  (@sec:mixing).

== What the refactor changes <sec:binning-future>

On 14 July 2026 the kd-tree was *removed* from the C++ and replaced by two
rectangular equal-statistics grids: Grid A over $(Q^2, x_B)$ and Grid B
over $(z, p_T^2)$, with the 4D bin being the product $A times B$ and Grid B
*independent of the Grid A cell*. Defaults are $8 times 7$ and
$5 times 5$, giving $1400$ bins. The multiplicity counts, the same/mixed
$m_(gamma gamma)$ spectra and the per-bin $phi_h$ histogram all share
$A times B$; $p_T$ broadening uses $A times (z "axis of" B)$.

This is a deliberate simplification --- factorized grids are reproducible,
have quotable global edges in every dimension, and remove the
non-determinism of @sec:binning-caveat. It also *sacrifices* the property
that motivated the kd-tree: because Grid B no longer adapts to the Grid A
cell, the $z$ and $p_T^2$ edges can no longer follow the strong correlation
between $p_T^2$ and $z$ that the nested tree captured (visible in
@tab:top-boxes, where the $p_T^2$ edges scale with $z$).

#wide-figure(
  "../figures/binning_grid_factorized.pdf",
  [The *new* factorized equal-statistics grid, from the LD#sub[2]
  outbending scan --- shown for reference, as it produced *none* of the
  results in this note. Cells failing a minimum-occupancy criterion are
  hatched. Contrast with @fig:kdtree: the grid lines here are global and
  rectangular, which is exactly the reproducibility gain, and exactly the
  adaptivity loss.],
  <fig:binning-factorized>,
  width: 88%,
)

None of the results in this note were produced with the new scheme. When the
next production runs, every number in @sec:results must be regenerated. The
bin-centre problem of @sec:binning-caveat has already been fixed in that
scheme (@sec:binning-caveat-fix) --- an equal-statistics rectangular grid
still has an unbounded outermost bin and so still needs count-weighted
abscissae, and it now emits them.
