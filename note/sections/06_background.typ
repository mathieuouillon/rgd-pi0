#import "../template/lib.typ": *

= Combinatorial background <sec:background>

Because the $pi^0$ candidate collection is effectively unselected in mass
(@sec:pairing), the $m_(gamma gamma)$ spectrum in every kinematic leaf
contains a $pi^0$ peak on a large combinatorial background of uncorrelated
photon pairs. Estimating and removing that background is the central
technical problem of the multiplicity measurement.

== Event mixing <sec:mixing>

The background shape is estimated by *event mixing*: pairing photons from
*different* events, which by construction cannot contain a true
$pi^0 -> gamma gamma$ and therefore trace the pure combinatorial shape
while retaining the single-photon kinematics and the detector acceptance.

=== The pool

Photons are pooled in bins of

$ ("bin") = (i_(Q^2), space i_(x_B), space i_"mult") , $ <eq:pool-key>

using the configuration product edges ($8 times 7$) and a photon
multiplicity class $i_"mult" in {1, 2, 3, >= 4}$, giving
$8 times 7 times 4 = 224$ pool bins. Matching on $Q^2$ and $x_B$ ensures
the mixed events have similar @DIS kinematics; matching on multiplicity
ensures the combinatorics are comparable.

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    table.header([*Property*], [*Implementation*]),
    [Depth], [$N = 50$ events per bin, @FIFO eviction],
    [Ordering], [Mix *then* insert --- an event never mixes with itself],
    [Pairing], [Current event's photons $times$ *every* photon in the pool],
    [Concurrency], [32 shards keyed on `bin % 32`, shared-mutex: concurrent
      reads while mixing, exclusive lock to update],
    [Storage], [5-double cache per photon, not full particles],
    [Output], [Histograms only --- mixed pairs are never written to the
      ntuple],
  ),
  caption: [The event-mixing pool. The mix-before-insert ordering is the
  detail that makes self-mixing impossible. Mixing is exhaustive rather
  than a fixed number of trials: once a bin is warm, each event's photons
  are paired against up to 50 events' worth of pooled photons, so the
  mixed-event statistics greatly exceed the same-event statistics and the
  background shape is not statistics-limited.],
) <tab:pool>

The finalisation step logs pool occupancy and *warns for every bin that
never reached depth $N$* --- a genuine quality-assurance feature, and the
right place to look for leaves where the background estimate is thin.

=== Two asymmetries in the mixed sample <sec:mixing-asymmetry>

#important-box(title: "The e–γ cut is one-sided for mixed pairs")[
  For a mixed pair only the *current-event* photon is tested against the
  $theta_(e gamma) > 8degree$ cut; the pooled photon is never re-tested
  against the current event's electron. This is *correct by construction*
  --- the pooled photon was already tested against its own event's
  electron, which is the physically meaningful reference.

  There is, however, a genuine inconsistency alongside it: *all* of the
  current event's photons are inserted into the pool, *including those that
  failed their own event's $e gamma$ cut*. So a photon that was rejected
  from the same-event spectrum can still contribute to the mixed-event
  spectrum. The effect is a small shape distortion of the background near
  the electron direction. It is not corrected and has not been quantified.
]

== Sideband normalisation and subtraction <sec:subtraction>

The mixed spectrum reproduces the background *shape* but carries an
arbitrary normalisation (it has far more entries). It is scaled to the
same-event spectrum in a high-mass sideband:

$ alpha = (sum_(m in "SB") N_"same" (m)) / (sum_(m in "SB") N_"mixed" (m)) ,
  quad "SB" = (0.17, 0.28) "GeV" , $ <eq:alpha>

and subtracted bin by bin:

$ S(m) = N_"same" (m) - alpha dot N_"mixed" (m) , quad
  sigma_S^2 (m) = N_"same" (m) + alpha^2 N_"mixed" (m) . $ <eq:subtraction>

The $m_(gamma gamma)$ histograms are 200 bins over $[0, 0.3]$ GeV, i.e.\
$1.5$ MeV per bin.

#note-box(title: "A single high-mass sideband")[
  The normalisation uses only the *upper* sideband $(0.17, 0.28)$ GeV. A
  low-mass sideband is available in principle but is not used --- at low
  $m_(gamma gamma)$ the spectrum is distorted by the opening-angle cut of
  @eq:theta-min, which removes precisely the small-angle pairs that
  populate it, so the mixed and same spectra are not expected to share a
  shape there.

  The uncertainty on $alpha$ is *not propagated*: @eq:subtraction treats it
  as exact. Given that the sideband contains a large number of counts this
  is a small effect, but it is a real omission and it is trivially fixable.
]

#warning-box(title: "Constants that look authoritative but are not")[
  The C++ defines `SIDEBAND_LOW/HIGH` $= 0.19 \/ 0.24$ GeV and
  `FIT_RANGE_LOW/HIGH` $= 0.10 \/ 0.16$ GeV. These appear in a log line and
  are *never applied to data*. The values that actually run are the Python
  ones: sideband $(0.17, 0.28)$, fit range $[0.08, 0.20]$. Do not quote the
  C++ constants as the analysis values --- they disagree with what is used.
]

== Cross-checks, implemented and not <sec:background-crosschecks>

A standalone diagnostic program applies five quantitative tests to the
subtraction on the LD#sub[2] sample:

#figure(
  table(
    columns: (auto, 1fr, auto),
    align: (left, left, left),
    table.header([*Test*], [*Definition*], [*Pass criterion*]),
    [Sideband $chi^2$],
      [$chi^2 \/ "ndf"$ of $N_"same" - alpha N_"mixed"$ in the sideband,
       with $"ndf" = n - 1$ for the normalisation constraint],
      [$< 5.0$],
    [Pull distribution], [Per-bin pulls compared against $cal(N)(0,1)$],
      [visual],
    [Scale stability], [$alpha$ versus kinematics], [visual],
    [Peak stability], [Gaussian fit to the subtracted peak],
      [$|mu - 135| <= 10$ MeV and $5 <= sigma <= 25$ MeV],
    [Sideband slope], [Linear slope of $N_"same" \/ (alpha N_"mixed")$],
      [$|"slope"| \/ sigma < 3$],
    [KS test], [Two-sample Kolmogorov--Smirnov in the sideband],
      [$alpha = 0.01$],
  ),
  caption: [Background-subtraction diagnostics. The sideband-slope test is
  the most informative of the set: with these statistics a $chi^2$ test is
  dominated by statistical fluctuations, whereas a *systematic slope* in
  the same-to-mixed ratio is the signature of a genuine shape mismatch.
  Failing bins are written to per-criterion PDF galleries for inspection.],
) <tab:mixing-validation>

#important-box(title: "These diagnostics are not wired into the measurement")[
  The validation program runs on LD#sub[2] only and produces plots. Its
  verdicts *do not feed back* into the multiplicity extraction: no leaf is
  rejected, and no systematic band is derived from them. The machinery to
  turn "the subtraction fails here" into an uncertainty exists and is
  unused.
]

#warning-box(title: "The rotated / swapped-background check does not exist")[
  Earlier drafts of this analysis describe a *rotated-background* or
  *swapped-photon* cross-check, in which one photon's azimuth is rotated by
  $180degree$ to destroy the decay correlation while preserving
  single-photon kinematics.

  *It is not implemented.* The two code branches that would consume it gate
  on a `mgg_swapped_all` histogram which is booked *nowhere* in the
  framework; both branches unconditionally print "Skipping". Any statement
  that the background was cross-checked with a rotated or swapped method
  would be false. Event mixing is the *only* background estimate in this
  analysis.
]
