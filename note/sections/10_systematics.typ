#import "../template/lib.typ": *

= Corrections and systematic uncertainties <sec:systematics>

== Corrections applied <sec:corrections>

#warning-box(title: "No physics correction is applied anywhere in the chain")[
  This is the most important statement in the note, and it is easy to miss
  because it is a statement about absence. Both codebases were searched
  exhaustively for acceptance, bin migration, radiative, contamination,
  charge normalisation, dead time, dilution, unfolding and efficiency
  corrections. *Every hit is a false positive* --- a histogram name, an
  unrelated binary, or a passive bank schema.

  The complete correction chain is:

  + *Combinatorial background* --- sideband-scaled mixed-event subtraction
    (@sec:subtraction). Applied to the multiplicity ratio *only*; not to
    $Delta chevron.l p_T^2 chevron.r$ (@sec:ptb-caveat) and not to $A_"LU"$
    (@sec:bsa-signal).
  + *Vertex correction* --- a track-level kinematic correction feeding the
    vertex *cut* (@sec:vertex). It corrects a cut boundary, not a yield.
  + *Luminosity, beam charge, dead time* --- not corrected, but *cancelled*
    by construction in the $N^"DIS"$-normalised ratio, to the extent they
    are common to the inclusive and semi-inclusive samples. They are, since
    both are counted from the same events.

  That is all. In particular there is *no Monte Carlo input of any kind* in
  this analysis: no acceptance model, no bin-migration unfolding, no
  radiative correction, and no simulation-based efficiency. The
  `CuSn_contamination` binary that measures Cu/Sn cross-contamination
  exists in the framework but is a *separate program* whose output does not
  feed the $pi^0$ chain.
]

== How much does the ratio actually cancel? <sec:ratio-cancellation>

$R_A$ and $Delta chevron.l p_T^2 chevron.r$ are *designed* so that the
dominant corrections cancel: numerator and denominator are reconstructed by
the same code, with the same cuts, on the same detector, so acceptance,
photon efficiency and radiative effects enter both alike and divide out to
first order. That cancellation is the load-bearing assumption of the whole
measurement, and it is *an assumption, not a measurement*.

It is also weaker here than it would be in a simultaneous-exposure
experiment, because *the targets were not exposed simultaneously*
(@sec:runs). Every ratio against LD#sub[2] compares runs taken days or weeks
apart. Four independent ways the cancellation fails are known to exist:

/ Run-to-run drift: the largest and most obvious. Beam current and
  position, target density, detector calibration, dead channels and trigger
  configuration all evolve over 515 runs, and none of it cancels between an
  LD#sub[2] block and a CxC or CuSn block recorded at a different time. The
  interleaving of @tab:run-blocks averages the long-term component down but
  removes nothing exactly.

/ Local occupancy: photon reconstruction efficiency depends on
  calorimeter occupancy, which differs between LD#sub[2] and Sn. The AI
  classifier explicitly uses neighbour multiplicity as input features
  (@sec:photon-id), so its response is occupancy-dependent *by
  construction*. This one would survive even a simultaneous exposure.

/ Within-leaf kinematics: a kd-tree leaf is a finite box, and the nuclear
  targets populate it differently. Any correction varying across the box
  therefore does not cancel exactly.

/ Vertex-dependent acceptance: Cu and Sn sit at different $z$
  (@tab:vz), so they see slightly different detector acceptance. LD#sub[2]
  is an extended cell, unlike the foils.

None of these has been quantified --- but the first, which is also the
largest, is *directly measurable from the existing data at no cost*, by
comparing the six LD#sub[2] blocks of @tab:run-blocks against one another.
Until that is done, the size of the run-to-run systematic on every number in
@sec:results is simply unknown.

#note-box(title: "One comparison that is clean")[
  Cu and Sn were exposed in the *same* runs (@tab:targets). The run-to-run
  systematic therefore cancels *exactly* in the Cu-to-Sn comparison, in a way
  it does not for any ratio against LD#sub[2]. Where the physics can be
  expressed as an $A$-dependence between Cu and Sn alone, it rests on much
  firmer ground than @tab:ptb's three-point fit against the deuterium
  reference.
]

== Systematic uncertainties evaluated: none <sec:no-systematics>

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, center, left),
    table.header([*Source*], [*Status*], [*Comment / how to get it*]),
    [*Run-to-run drift*], [not evaluated],
      [*Probably the largest.* The targets were taken in separate runs
      (@sec:runs), so nothing time-dependent cancels. Measurable today by
      ratioing the six LD#sub[2] blocks of @tab:run-blocks against each
      other. Cancels exactly for Cu-vs-Sn only.],
    [Background subtraction], [not evaluated],
      [The no-subtraction variant already exists for every leaf
      (@sec:no-subtraction); differencing it gives a bound *today*.],
    [Sideband definition], [not evaluated],
      [Vary $(0.17, 0.28)$; cheap.],
    [Mass window], [not evaluated],
      [Vary the $plus.minus 3 sigma$ integration window.],
    [Photon ID], [not evaluated],
      [Blocked: the skim freezes the photon selection (@sec:skim-cost).
      Requires a re-skim.],
    [@RGA model fallback], [not evaluated],
      [@sec:photon-id. Needs a model trained on @RGD, or a data-driven
      efficiency study.],
    [Acceptance], [not evaluated], [Requires @MC.],
    [Bin migration], [not evaluated], [Requires @MC.],
    [Radiative], [not evaluated], [Requires @MC.],
    [Cu/Sn contamination], [not evaluated],
      [A separate binary measures it; not propagated.],
    [Energy scale], [not evaluated],
      [The $4.4%$ peak shift of @sec:invariant-mass is a direct handle.],
    [Beam polarization], [not propagated],
      [@BSA only; $3.5%$ correlated scale, one line to add
      (@sec:bsa-polarization).],
    [Polarity comparison], [not done],
      [Inbending data exists and is unused (@sec:targets).],
    [False asymmetry], [not evaluated], [@BSA only.],
  ),
  caption: [Systematic uncertainties. *Every entry is unevaluated.* The
  entries in the upper block are accessible with the data and code already
  in hand; those requiring @MC are a larger undertaking.],
) <tab:systematics>

#important-box(title: "The statistical precision is misleading")[
  The median statistical uncertainty on $R_A$ is $1.80%$, reaching
  $plus.minus 0.001$ in $z$-binned averages (@sec:results). *These errors
  are not the uncertainty on the measurement.* They are the uncertainty on
  the measurement *given* that every correction in @tab:systematics is
  exactly zero and every cancellation in @sec:ratio-cancellation is
  perfect.

  Systematics will dominate the final uncertainty entirely. Quoting the
  statistical error alone --- as every plot and result file in the current
  production does --- overstates the precision by an unknown but certainly
  large factor. No result in this note should be interpreted as a
  measurement with a $0.1%$ uncertainty.
]

== Suggested order of attack <sec:systematics-plan>

Ordered by (impact $times$ ease), the cheapest first:

+ *Measure the run-to-run systematic* by ratioing the six LD#sub[2] blocks
  of @tab:run-blocks against one another. No new data, no simulation, no new
  code beyond splitting the LD#sub[2] sample by run. It is plausibly the
  dominant systematic and it is currently a blank.
+ *Propagate the polarization uncertainty* and obtain the measured Møller
  value. One line and one number; unblocks the @BSA entirely.
+ *Difference the subtracted and no-subtraction variants* into a background
  systematic. Both already exist for all $1450$ leaves.
+ *Re-run with the count-weighted abscissae* now emitted by the C++
  (@sec:binning-caveat). The code is in place; this needs only the farm
  pass. It plausibly changes the physics conclusion about the strength of
  the attenuation.
+ *Vary the sideband and the integration window.* Pure re-running.
+ *Add a dilution factor to the @BSA* from the $m_(gamma gamma)$ fits the
  multiplicity analysis already performs.
+ *Wire the mixing diagnostics* (@tab:mixing-validation) into a per-leaf
  quality mask.
+ *Run the inbending data* as a detector-systematic cross-check.
+ *Resolve the $p_T$-broadening mass window* (@sec:ptb-caveat) --- ideally
  by accumulating the moments per $m_(gamma gamma)$ bin.
+ *Commission an @RGD photon-ID model* and re-skim.
+ *Build the @MC chain* for acceptance, migration and radiative
  corrections.
