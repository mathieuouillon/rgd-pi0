#import "../template/lib.typ": *

= Summary and outlook <sec:summary>

== What this analysis does <sec:summary-what>

Neutral pions are reconstructed via $pi^0 -> gamma gamma$ in @RGD
outbending data at $E_"beam" = 10.53$ GeV, on an LD#sub[2] / CxC / Cu / Sn
target set taken in separate, block-alternating runs (@sec:runs). A C++17
skim (`stageA_skim`) applies electron, @DIS and AI-based photon selection and
reduces the data $approx 100 times$; `make_grid` freezes a factorized
equal-statistics $(Q^2, x_B) times (z, p_T^2)$ grid (@sec:binning); a second
C++ pass (`stageB_bin`) reconstructs $pi^0$, computes @SIDIS kinematics, builds
the mixed-event background and fills four flat trees; Python (`pi0.*`) extracts
the physics. Every stage stamps its configuration and inputs into its output
and refuses to build physics on an input it cannot vouch for
(@sec:provenance-gaps).

Three observables are produced: the multiplicity ratio $R_A$ (yields from a
sideband-subtracted mixed-event subtraction, normalised to the @DIS count per
$(Q^2, x_B)$ cell), the $p_T$-broadening $Delta chevron.l p_T^2 chevron.r$
(compensated accumulator moments on the 3D grid), and the beam-spin asymmetry
$A_"LU"^(sin phi_h)$ (single-parameter $sin phi_h$ fit).

The chain has been run end to end on real data for the first time
(@sec:results). That run is diagnostic: it read a $2 times 10^6$-event
prefix of each SIDIS-train file, and its photons carry the @RGA fallback of
@sec:photon-id, so no number from it may be quoted. What it establishes is that
the machinery is correct --- the $p_T$ broadening comes out positive,
$A$-ordered and consistent with an independent prior analysis --- and that the
remaining work is statistics, correction and validation rather than
reconstruction.

== What blocks a measurement <sec:summary-blockers>

#figure(
  table(
    columns: (auto, 1fr, auto),
    align: (left, left, center),
    table.header([*Blocker*], [*What it is*], [*Section*]),
    [Statistics],
      [The results are a diagnostic $2 times 10^6$-event/file prefix, not the
       full luminosity, so $N_"DIS"$ is short by an unknown amount and every
       yield is flagged unpublishable. Removing the cap reuses the same
       machinery.],
      [@sec:statistics],
    [Photon-ID model],
      [@RGD runs fall outside every GBT model range; all photons are classified
       by an @RGA inbending pass-1 fallback. The single largest photon
       systematic, and unquantified.],
      [@sec:photon-id],
    [No systematics],
      [Not one systematic is evaluated, and no acceptance, bin-migration or
       radiative correction is applied --- there is no simulation chain yet.],
      [@sec:systematics],
    [Beam polarization],
      [`pi0.bsa` refuses to run without a measured $P$; there is no placeholder.
       $A_"LU"$ is not yet extracted.],
      [@sec:bsa-polarization],
    [Run-to-run drift],
      [The targets were taken in separate runs, so nothing time-dependent
       cancels in $R_A$. Plausibly the dominant systematic, and measurable
       today from the six LD#sub[2] blocks.],
      [@sec:ratio-cancellation],
    [$p_T$-broadening sample],
      [The moments are sideband-subtracted but still accumulated over a
       background-dominated window, which pulls
       $Delta chevron.l p_T^2 chevron.r$ toward zero by a target-dependent
       factor.],
      [@sec:ptb-caveat],
  ),
  caption: [Open items, roughly in order of impact on the physics. All are
  correction and validation; none is a reconstruction defect.],
) <tab:blockers>

== Recommended next steps <sec:outlook>

Immediate --- no new data.

+ Run the full-luminosity production. Remove the $2 times 10^6$-event cap
  and skim the complete SIDIS train; this is the same command with a larger
  budget and is what turns the diagnostic ratios into a statistics-limited
  measurement.
+ Measure the run-to-run systematic from the six LD#sub[2] blocks of
  @tab:run-blocks. Same target, six times across the period: any difference is
  pure time-dependent systematic. This is the largest blank in
  @tab:systematics and the cheapest thing on this list.
+ Difference the subtracted and no-subtraction $R_A$ variants into a first
  background systematic (@sec:no-subtraction).

Short term --- small code or data changes.

+ Project $R_A$ against $nu = Q^2 \/ (2 M_p x_B)$ at fixed $z$. The
  $(Q^2, x_B)$ grid already contains it (@sec:why-kinematics), so this needs
  no new data, no new binning and no simulation --- and at
  $chevron.l nu chevron.r approx 5.1$ GeV, where $nu \/ kappa$ is comparable
  to the nuclear size (@tab:lengths), the $nu$ dependence is what separates
  energy loss from pre-hadron absorption.
+ Obtain the measured Møller polarization, record it in `config/cuts.json`,
  and run `pi0.bsa` with its uncertainty propagated.
+ Vary the sideband, the mass window and the integration window into a
  systematic band.

Medium term.

+ Commission an @RGD photon-ID model, or measure the photon efficiency from
  data, and re-skim --- this is the highest-value fix, since it removes the
  single largest unquantified photon systematic.
+ Run the inbending dataset --- both for statistics and as the
  detector-systematic cross-check the analysis currently lacks entirely.

Longer term.

+ Build the @MC chain for acceptance, bin-migration and radiative
  corrections. This is the largest remaining piece and the one that turns
  the current ratios into a measurement.
+ Extend to the charged-pion channels for a flavour-separated comparison
  using the same framework, and cross-check $R_A$ between the $pi^0$ and
  $pi^plus.minus$ channels --- their systematics are almost disjoint, which
  makes the comparison unusually powerful.

== Closing assessment <sec:closing>

The infrastructure is in good shape: the skim / grid / bin / extract split is
clean, every cut and every bin edge is version-controlled and hashed into the
outputs, the event-mixing pool is carefully implemented and instrumented, the
$phi_h$ convention is verifiably Trento, and the whole chain now runs end to
end on data with a provenance guard that refuses to quote what it should not.

The gap between that infrastructure and a publishable measurement is almost
entirely in statistics, correction and validation, not in reconstruction.
That is a comfortable position to be in --- but the gap is wide. The diagnostic
run must be replaced by a full-luminosity production before any yield is a
number; the photon-ID fallback must be closed before the photon efficiency is
trusted; and because the targets were taken in separate runs, the run-to-run
systematic is both plausibly the largest uncertainty and a complete blank ---
yet the six LD#sub[2] blocks of @tab:run-blocks measure it directly, from data
already on disk. Beyond those, the dataset has one clear physics opportunity
that costs nothing to take: $R_A$ against $nu$, which the existing
$(Q^2, x_B)$ grid already contains and which is what actually distinguishes the
two hadronization mechanisms at these energies.
