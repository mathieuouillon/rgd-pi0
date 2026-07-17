#import "../template/lib.typ": *

= Summary and outlook <sec:summary>

== What this analysis does <sec:summary-what>

Neutral pions are reconstructed via $pi^0 -> gamma gamma$ in @RGD
outbending data at $E_"beam" = 10.53$ GeV, on an LD#sub[2] / CxC / Cu / Sn
target set taken in *separate, block-alternating runs* (@sec:runs). A C++20 skim
applies electron, @DIS and AI-based photon selection and reduces the data
$approx 100 times$; a second C++ pass reconstructs $pi^0$, computes @SIDIS
kinematics and fills histograms binned on an adaptive four-dimensional
kd-tree of $approx 1450$ equal-occupancy leaves; Python extracts the physics.

Three observables are produced: the multiplicity ratio $R_A$ (yields from a
sideband-scaled mixed-event subtraction and a $plus.minus 3 sigma$ window
sum, normalised to the @DIS count per $(Q^2, x_B)$ box), the
$p_T$-broadening $Delta chevron.l p_T^2 chevron.r$ (accumulator moments on a
3D tree), and the beam-spin asymmetry $A_"LU"^(sin phi_h)$ (single-parameter
$sin phi_h$ fit in 12 $phi_h$ bins).

The statistical reach is the analysis's headline achievement: $approx 15.7$M
$pi^0$ on LD#sub[2] alone, $approx 1-2%$ statistical precision per
four-dimensional leaf, and all the expected hadronization signatures
present and correctly $A$-ordered (@sec:results).

== What blocks publication <sec:summary-blockers>

#figure(
  table(
    columns: (auto, 1fr, auto),
    align: (left, left, center),
    table.header([*Blocker*], [*What it is*], [*Section*]),
    [Run-to-run drift],
      [The targets were taken in *separate runs*, so nothing time-dependent
       cancels in $R_A$. Plausibly the dominant systematic, and entirely
       unmeasured --- though measurable today from the six LD#sub[2] blocks.],
      [@sec:ratio-cancellation],
    [Bin abscissae],
      [$52.6%$ of leaves are reported at a geometric box midpoint that is
       not where the data is --- for $6.9%$ of them, a point no event could
       occupy. The *depth* of the attenuation is not a measurement until
       this is fixed. *Fixed in code; needs the re-run.*],
      [@sec:binning-caveat],
    [No systematics],
      [Not one systematic is evaluated. Statistical errors of $0.1%$ are
       being quoted on an uncorrected measurement.],
      [@sec:systematics],
    [Beam polarization],
      [$P = 0.85$ is a self-declared placeholder; all @BSA numbers scale by
       $0.85\/P_"true"$.],
      [@sec:bsa-polarization],
    [Photon-ID model],
      [@RGD runs fall outside every model range; all photons are classified
       by an @RGA *inbending* pass-1 model via a silent fallback.],
      [@sec:photon-id],
    [$p_T$-broadening sample],
      [Moments accumulated over a $plus.minus 200$ MeV window with no
       background subtraction.],
      [@sec:ptb-caveat],
    [Provenance],
      [Production configuration and binning tables are not in any
       repository; results are untracked; run lists exist but are dead
       code.],
      [@sec:provenance-gaps],
    [Railed @BSA fits],
      [Two leaves at the fit bound with $sigma approx 10^(-6)$ poison any
       naive average.],
      [@sec:results-bsa],
  ),
  caption: [Open items, roughly in order of impact on the physics. The
  first two are the ones that change conclusions; the rest change numbers.],
) <tab:blockers>

== Recommended next steps <sec:outlook>

*Immediate --- no new data, no new code.*

+ *Measure the run-to-run systematic* from the six LD#sub[2] blocks of
  @tab:run-blocks. Same target, six times across the period: any difference
  is pure time-dependent systematic. This is the largest blank in
  @tab:systematics and the cheapest thing on this list.
+ Filter the railed @BSA leaves and add a fit-quality mask to the producing
  script.
+ Difference the existing subtracted and no-subtraction $R_A$ variants into
  a first background systematic (@sec:no-subtraction).
+ Archive the kd-tree `.txt` files and the production configuration
  alongside the results; put `results/`, `analyses/` and `common/` under
  version control.
+ Delete the $approx 5.2$ GB of orphaned caches (@sec:ptb-cache) and the
  superseded product-grid result files, so they cannot be picked up by
  mistake.

*Short term --- small code changes.*

+ *Re-run the production* to pick up the count-weighted abscissae, which the
  C++ now emits (@sec:binning-caveat-fix). The code is done; only the farm
  pass is outstanding. This is the step most likely to move the physics
  conclusion, since it is what fixes the depth of the attenuation.
+ *Project $R_A$ against $nu = Q^2 \/ (2 M_p x_B)$* at fixed $z$. The
  $(Q^2, x_B)$ grid already contains it (@sec:why-kinematics), so this needs
  no new data, no new binning and no simulation --- and at
  $chevron.l nu chevron.r approx 5.1$ GeV, where $nu \/ kappa$ is comparable
  to the nuclear size (@tab:lengths), the $nu$ dependence is what separates
  energy loss from pre-hadron absorption.
+ Obtain the measured Møller polarization and propagate its uncertainty.
+ Accumulate $sum p_T^2$ and $sum p_T^4$ per $m_(gamma gamma)$ bin so the
  moments can be sideband-subtracted like the yields.
+ Apply a per-leaf dilution factor to $A_"LU"$ from the $m_(gamma gamma)$
  fits already performed.

*Medium term.*

+ Run the inbending dataset --- both for statistics and as the
  detector-systematic cross-check the analysis currently lacks entirely.
+ Vary the sideband, the mass window and the integration window into a
  systematic band.
+ Commission an @RGD photon-ID model, or measure the photon
  efficiency from data; re-skim.

*Longer term.*

+ Build the @MC chain for acceptance, bin-migration and radiative
  corrections. This is the largest remaining piece and the one that turns
  the current ratios into a measurement.
+ Propagate the Cu/Sn cross-contamination already measured by the separate
  binary.
+ Extend to the charged-pion channels for a flavour-separated comparison
  using the same framework, and cross-check $R_A$ between the $pi^0$ and
  $pi^plus.minus$ channels --- their systematics are almost disjoint, which
  makes the comparison unusually powerful.

== Closing assessment <sec:closing>

The infrastructure is in good shape: the skim/analysis/extraction split is
clean, the event-mixing pool is carefully implemented and instrumented, the
$phi_h$ convention is verifiably Trento, and the adaptive binning is a
genuinely good idea well executed. The statistical precision available is
extraordinary for this channel.

The gap between that infrastructure and a publishable measurement is almost
entirely in *correction and validation*, not in reconstruction. That is a
comfortable position to be in --- but the gap is wide, and the current
statistical errors give a misleading impression of how narrow it is.

Two next actions dominate the rest. The bin abscissae are now fixed in code
and need only a re-run (@sec:binning-caveat-fix); until that lands, the
*depth* of the attenuation is not a measurement at all, only a mislabelled
integral, and nothing quantitative should be built on it. And because the
targets were taken in separate runs, the run-to-run systematic is both
plausibly the largest uncertainty and a complete blank --- yet the six
LD#sub[2] blocks of @tab:run-blocks measure it directly, from data already
on disk. Neither requires new beam time, new simulation, or a new method.

Beyond those, the dataset has one clear physics opportunity that costs
nothing to take: $R_A$ against $nu$, which the existing $(Q^2, x_B)$ grid
already contains and which is what actually distinguishes the two
hadronization mechanisms at these energies.
