#import "../template/lib.typ": *

= Introduction <sec:intro>

== Physics motivation <sec:motivation>

Hadronization --- the transition of a struck quark into the colour-neutral
hadrons that reach a detector --- is a genuinely non-perturbative part of
@QCD, and one of the least constrained. In elementary $e^+ e^-$ or $e p$
collisions the process can only be characterised through @FF:pl extracted
from data; the space--time development itself is not observable. Nuclei
change that. When the hard scattering happens on a nucleon bound inside a
nucleus, the struck quark and whatever colour-singlet pre-hadron it forms
must traverse a few femtometres of nuclear matter before escaping. The
nucleus acts as a femtometre-scale detector of the evolution, and the
comparison of hadron production on a heavy nucleus $A$ against a light
reference isolates the medium effect.

Two length scales control what is seen. The _production length_ $l_p$ is
the distance over which the struck quark remains a coloured object losing
energy by induced gluon radiation; the _formation length_ $l_f$ is the
distance after which the pre-hadron has evolved into a physical hadron with
its full cross section. If $l_p$ and $l_f$ are large compared with the
nuclear radius $R_A approx 1.2 A^(1\/3)$ fm, the relevant mechanism is
partonic energy loss and the observable consequence is a shift of the
hadron spectrum toward lower energy. If instead the pre-hadron forms inside
the nucleus, it is absorbed by hadronic final-state interactions and the
consequence is a genuine reduction of the yield. Distinguishing the two ---
and thereby measuring $l_p$ and $l_f$ --- is the central goal of nuclear
@SIDIS.

=== Why these kinematics <sec:why-kinematics>

Both lengths scale with the energy available to the struck quark, $nu$,
divided by the string tension $kappa approx 1$ GeV/fm: dimensionally
$l_f tilde nu \/ kappa$, times a model-dependent factor of order $(1-z)$
that is precisely what a measurement differential in $z$ is meant to pin
down. The measurement is interesting exactly when that length is comparable
to the nucleus --- large enough that the hadron sometimes escapes before
forming, small enough that it sometimes does not.

@RGD sits in that window. The @DIS cuts of @tab:dis-cuts bound the energy
transfer analytically: $W > 2$ GeV gives
$nu > (4 - M_p^2 + Q^2) \/ 2 M_p$, and $y < 0.85$ gives
$nu < 8.95$ GeV, so

$ 2.2 lt.tilde nu lt.tilde 9.0 "GeV" , $ <eq:nu-range>

with a count-weighted mean of $chevron.l nu chevron.r approx 5.1$ GeV over
the populated bins.#footnote[Computed from the count-weighted
$chevron.l Q^2 chevron.r$ and $chevron.l x_B chevron.r$ of the @BSA output,
the only current file carrying true per-bin means (@sec:binning-caveat).]

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, right, right, right),
    table.header([*Target*], [*$A$*], [*$R_A = 1.2 A^(1\/3)$*],
                 [*$chevron.l L chevron.r approx 3 R_A \/ 4$*]),
    [C],  [12],    [$2.75$ fm], [$2.1$ fm],
    [Cu], [63.5],  [$4.79$ fm], [$3.6$ fm],
    [Sn], [118.7], [$5.90$ fm], [$4.4$ fm],
  ),
  caption: [Nuclear sizes, against the hadronization length scale
  $nu \/ kappa approx 5$ fm at $chevron.l nu chevron.r = 5.1$ GeV. The two
  are *comparable* --- neither limit applies cleanly --- which is what makes
  the measurement worth doing and also what makes it hard: energy loss and
  absorption both contribute, and only their different dependence on $z$,
  $nu$ and $A$ separates them. $chevron.l L chevron.r$ is the mean path from
  a uniformly-sampled production point to the surface.],
) <tab:lengths>

One consequence for how this note is binned: $nu$ is not an independent
axis. It is fixed by $Q^2$ and $x_B$ through $nu = Q^2 \/ (2 M_p x_B)$, so
the $(Q^2, x_B)$ grid of @sec:binning *is* a $nu$ binning, read in a rotated
frame. Since the balance between energy loss and pre-hadron absorption is
governed by $nu$ --- through the $nu \/ kappa$ of @tab:lengths --- the
$nu$ dependence is the most direct handle this dataset has on which
mechanism dominates, and it is available without any new binning.

The @CLAS12 spectrometer is well matched to this problem: the @RGD run
period delivered a $10.53$ GeV electron beam at high luminosity onto a
liquid-deuterium cell and a set of solid nuclear foils, over a large
acceptance and with the statistics documented in @tab:statistics.

#important-box(title: "The targets were not exposed simultaneously")[
  Each @RGD run used *one* target. The reference and the nuclear foils were
  taken in *separate runs*, alternating in blocks over the run period
  (@sec:runs) --- there is no run in which LD#sub[2] and a nuclear target
  were exposed together.

  This matters for the whole measurement, because the ratio of
  @eq:multiplicity is then formed from data taken *at different times*, not
  under identical conditions. Luminosity drift, detector efficiency changes
  and calibration shifts between an LD#sub[2] block and a neighbouring
  nuclear block do *not* cancel in $R_A$ --- they enter it directly. The
  time-dependent systematics that a simultaneous multi-target exposure would
  have removed are, here, a real and so far unevaluated uncertainty
  (@sec:ratio-cancellation).

  The running was at least *interleaved*: the beam returned to LD#sub[2] six
  separate times across the period, so the reference samples much the same
  span of running conditions as the nuclear targets rather than sitting at
  one end of it. That averages long-term drift down; it does not cancel it.
  It also leaves a direct handle on the size of the effect, since the six
  LD#sub[2] blocks can be compared against each other
  (@sec:ratio-cancellation).
]

This note documents the $pi^0$ channel. The $pi^0 -> gamma gamma$ decay
(branching ratio $98.82%$) is reconstructed entirely in the electromagnetic
calorimeters, so the measurement is independent of the hadron
time-of-flight and tracking systems used for charged pions. It therefore
provides a systematically independent cross-check on the charged-pion
results, with a completely different set of dominant uncertainties: the
$pi^0$ has no tracking or $Delta t$ particle-identification uncertainty, but
it does carry a combinatorial background under the invariant-mass peak that
charged pions do not.

== Observables <sec:observables>

Three observables are extracted, in order of decreasing maturity.

=== Multiplicity ratio

The primary observable is the nuclear multiplicity ratio, the yield of
$pi^0$ per @DIS event on a nuclear target $A$ normalised to the same
quantity on deuterium:

$ R_A^(pi^0) (Q^2, x_B, z, p_T^2)
  = (N_A^(pi^0) (Q^2, x_B, z, p_T^2) \/ N_A^"DIS" (Q^2, x_B))
    / (N_D^(pi^0) (Q^2, x_B, z, p_T^2) \/ N_D^"DIS" (Q^2, x_B)) $ <eq:multiplicity>

Normalising to the @DIS count in the *same* $(Q^2, x_B)$ cell does two
things, and the second is the important one.

The bookkeeping part is that it removes the integrated luminosity, the beam
charge and the dead time, provided these are common to the inclusive and
semi-inclusive samples --- which they are, since both are counted from the
same events.

The physics part is that it removes the *initial state*. A nucleus is not a
free collection of nucleons: its quark distributions are modified relative
to the deuteron by shadowing at low $x_B$, antishadowing, and the EMC effect
at large $x_B$. Those are initial-state effects and have nothing to do with
hadronization --- but they change $N_A^(pi^0)$ and $N_A^"DIS"$ *by the same
factor*, because both are measured on the same nucleus at the same
$(Q^2, x_B)$ and so sample the same nuclear parton densities. In the ratio
of ratios they cancel. Writing schematically
$N_A^h prop f_(q\/A)(x_B, Q^2) times D_(q -> h)^A (z)$ and
$N_A^"DIS" prop f_(q\/A)(x_B, Q^2)$, the nuclear @PDF drops out and

$ R_A^(pi^0) approx D_(q -> pi^0)^A (z, p_T^2) \/ D_(q -> pi^0)^D (z, p_T^2) , $ <eq:ra-ff>

a ratio of *effective fragmentation functions*. This is what makes $R_A$ a
probe of hadronization rather than of nuclear structure, and it is the
reason the normalisation must be taken per $(Q^2, x_B)$ cell rather than
inclusively (@sec:normalisation). The cancellation is exact only to the
extent that the factorization above holds, which is itself part of what the
measurement tests.

$R_A < 1$ signals attenuation; $R_A > 1$ signals enhancement, expected at
low $z$ where the energy removed from the leading hadron reappears as
additional soft production.

=== Transverse-momentum broadening

The medium imparts transverse momentum to the propagating parton through
multiple soft scattering. Each scattering adds in quadrature, so the
transverse momentum executes a random walk and the *squared* momentum
accumulates linearly along the path. The resulting broadening,

$ Delta chevron.l p_T^2 chevron.r = chevron.l p_T^2 chevron.r_A - chevron.l p_T^2 chevron.r_D approx hat(q) dot chevron.l L chevron.r , $ <eq:broadening-intro>

measures the transport coefficient $hat(q)$ --- the mean squared transverse
momentum acquired per unit path length --- times the distance travelled in
the medium. It is theoretically cleaner than $R_A$: a difference rather than
a ratio, insensitive to the absolute normalisation, and free of the
fragmentation function to the extent that the intrinsic $p_T$ is
target-independent.

The $A$-dependence is the reason to measure it on three nuclei. If the
parton simply traverses a nucleus of uniform density, then
$chevron.l L chevron.r prop R_A prop A^(1\/3)$ (@tab:lengths), and with
$hat(q)$ a property of the medium rather than of the nucleus,

$ Delta chevron.l p_T^2 chevron.r prop A^(1\/3) . $ <eq:A13>

Departures from the $1\/3$ exponent are therefore informative rather than
embarrassing: they point to a density profile that is not uniform, to a path
length shortened because the pre-hadron formed early, or to an $hat(q)$ that
is not universal. The measured exponent is discussed in @sec:results-ptb.

=== Beam-spin asymmetry

With a longitudinally polarised beam and an unpolarised target, the
$sin phi_h$ moment of the helicity-difference cross section,
$A_"LU"^(sin phi_h)$, probes the twist-3 sector and is sensitive to
quark--gluon correlations. Its nuclear dependence is essentially
unmeasured. Because $A_"LU"$ is a ratio of yields taken _within a single
$phi_h$ bin_, any helicity-independent acceptance cancels exactly, making it
the most acceptance-robust of the three.

#important-box(title: "Status of this note")[
  This note documents the analysis *as it is currently implemented*, and is
  explicit about the difference between what the code computes and what is
  ready to be published. Several observables carry unresolved issues that
  are stated in place rather than deferred to a footnote:

  - the beam polarization is a *placeholder value*, so every @BSA number is
    provisional (@sec:bsa-polarization);
  - the reported kinematic position of the widest bins is *not* where the
    data sits, which affects how @tab:R-vs-z must be read
    (@sec:binning-caveat). This is now fixed in code, but the results
    predate the fix and still need the re-run (@sec:binning-caveat-fix);
  - the targets were taken in *separate runs*, so the run-to-run systematic
    does not cancel in the ratio and has not been measured
    (@sec:ratio-cancellation);
  - *no* acceptance, bin-migration, or radiative correction is applied
    anywhere in the chain, and no systematic uncertainty has been evaluated
    (@sec:systematics).

  The statistical precision is already at the $1-2%$ level. Systematics,
  none of which are yet quantified, will dominate the final uncertainty.
]

== What this measurement adds <sec:previous>

The qualitative picture --- attenuation growing with $A$ and with $z$, and a
$p_T$ broadening positive and increasing with $A$ --- is established. What
is not established is anything quantitative enough to separate partonic
energy loss from pre-hadron absorption, because that separation lives in the
*differential* behaviour: how the attenuation evolves with $z$ at fixed
$nu$, and with $nu$ at fixed $z$. Holding one fixed while scanning the other
demands both statistics and kinematic granularity.

That is what this dataset offers. The statistics ($15.7 times 10^6$ $pi^0$
on LD#sub[2] alone, @tab:statistics) support a genuinely four-dimensional
binning in $(Q^2, x_B, z, p_T^2)$ at $1-2%$ statistical precision per bin
(@sec:results), on three nuclei spanning an order of magnitude in $A$, at a
$nu$ that sits where the two mechanisms compete rather than where one of
them trivially wins (@sec:why-kinematics). The $pi^0$ channel adds a further
handle: its systematics are almost disjoint from the charged-pion channels
measured in the same run group with the same framework, so the two form a
strong mutual cross-check (@sec:outlook).
