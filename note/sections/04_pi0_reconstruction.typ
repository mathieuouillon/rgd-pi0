#import "../template/lib.typ": *

= Neutral pion reconstruction <sec:pi0-reco>

== Invariant mass <sec:invariant-mass>

For a photon pair the invariant mass is formed directly from the summed
four-momenta,

$ m_(gamma gamma) = sqrt((p_(gamma_1) + p_(gamma_2))^2)
                  = sqrt(2 E_(gamma_1) E_(gamma_2) (1 - cos theta_12)) , $ <eq:mgg>

with $E_gamma = |arrow(p)_gamma|$ and $theta_12$ the opening angle. The
nominal mass is $m_(pi^0) = 0.1349768$ GeV.

#important-box(title: "No photon energy correction is applied")[
  There is no photon energy correction, no calorimeter scale factor, no
  $pi^0$-mass constraint and no vertex constraint anywhere in the chain.
  Photon four-momenta pass from `REC::Particle` to the invariant mass
  untouched.

  This is visible in the data. The fitted peak position, median over the
  populated bins of a representative run, is

  #align(center)[
    #text(size: 9.5pt)[
      $mu = 129.06$ (LD#sub[2]), $129.50$ (CxC),
      $129.25$ (Cu), $129.40$ (Sn) MeV,
    ]
  ]

  i.e.\ *$4.8 - 5.9$ MeV below* the nominal $134.98$ MeV --- a $4.4%$
  energy-scale deficit, consistent across targets to better than $1$ MeV.
  The fitted width is $sigma approx 15.5 - 15.7$ MeV throughout.

  Because the $plus.minus 3 sigma$ integration window is centred on the
  *fitted* $mu$ rather than on the nominal mass (@sec:yield-extraction),
  the shift is largely absorbed and does not by itself bias the yield. Its
  target-independence further protects the ratio. It is nonetheless a clear
  indication that the calorimeter energy scale is uncalibrated for this
  sample, and the $approx 1$ MeV target-to-target spread is a real, if small,
  systematic on $R_A$ that has not been evaluated.
]

== Pair selection <sec:pairing>

Candidates are built by *greedy exclusive pairing*:

#figure(
  block(inset: (x: 6pt, y: 4pt))[
    #set text(size: 9.5pt)
    #set par(justify: false)
    + Consider all unused photon pairs $(i, j)$.
    + Among those with $|m_(i j) - m_(pi^0)| < Delta m$, take the pair
      minimising $|m_(i j) - m_(pi^0)|$.
    + Emit it; mark both photons used; repeat from step 1.
    + Stop when no admissible pair remains.
  ],
  caption: [The `find_all_unique_pi0_candidates` algorithm. Each photon is
  consumed by at most one $pi^0$, so an $N$-photon event yields at most
  $floor(N \/ 2)$ candidates. The loop is greedy and best-first, not a
  global optimum: an assignment maximising the total number of candidates,
  or the joint likelihood, could differ.],
) <fig:pairing>

#warning-box(title: "The ±200 MeV window is a pairing cut, not a mass selection")[
  The pairing cut is $|m_(gamma gamma) - 0.1349768| < Delta m$ with
  $Delta m$ = `pairing.mass_window_gev` $= 0.2$ GeV --- *deliberately wide*.
  Since $m_(gamma gamma) >= 0$ always, it reduces to

  $ m_(gamma gamma) < 0.335 "GeV" $

  with *no lower bound at all*: essentially every $gamma gamma$ pair below
  $335$ MeV enters the "$pi^0$ candidate" collection, so it is *not a
  mass-selected sample* in any meaningful sense.

  This is *by design, not by accident* --- the signal is extracted offline by
  fitting and subtracting the $m_(gamma gamma)$ spectrum (@sec:background),
  which requires the sidebands to be present; narrowing the window would
  destroy the background estimate. Two consequences must be stated plainly:

  - the $p_T$-broadening accumulators, which are filled from this collection
    *without* any offline subtraction, are computed over a largely
    combinatorial sample --- see @sec:ptb-caveat;
  - the greedy "closest to $m_(pi^0)$" ordering (@fig:pairing) is doing more
    work than the window is, since the window barely excludes anything.
]

== Pair-level cuts <sec:pair-cuts>

A further set of cuts is applied in the analysis stage, gating both the
$m_(gamma gamma)$ spectra and the physics histograms.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    table.header([*Cut*], [*Value*], [*Purpose*]),
    [Electron--photon angle], [$theta_(e gamma) > 8.0degree$],
      [Rejects photons collinear with the scattered electron
      (bremsstrahlung, radiative tails).],
    [Mass floor], [$m_(gamma gamma) > 0.001$ GeV],
      [Rejects nearly collinear pairs with unphysically small mass.],
    [Opening angle], [$theta_(gamma gamma) > theta_"min" (p)$],
      [Rejects merged clusters; @eq:theta-min.],
    [Energy fraction], [$0 < z < 1$], [Physical hadrons only; applies to the
      multiplicity fills.],
  ),
  caption: [Pair-level cuts. All are applied identically to same-event and
  mixed-event pairs, with one exception noted in @sec:mixing-asymmetry.],
) <tab:pair-cuts>

At high $pi^0$ momentum the two decay photons become increasingly
collinear and their calorimeter clusters merge. The opening-angle threshold
is therefore momentum-dependent:

$ theta_"min" (p) = 17.561 dot exp(-0.756 p) + 1.0 quad ["degrees"] , $ <eq:theta-min>

where $p = |arrow(p)_(gamma_1) + arrow(p)_(gamma_2)|$ is the *pair*
momentum. The threshold falls from $18.56degree$ at $p = 0$ through
$4.87degree$ at $p = 2$ GeV to an asymptote of $1.0degree$.

== SIDIS kinematics <sec:sidis-kinematics>

The virtual-photon frame is built per event from
$q = k_"beam" - k'_e$:

$ hat(q) = arrow(q) \/ |arrow(q)| , quad
  hat(y) = (arrow(q) times arrow(k)_"beam") \/ |arrow(q) times arrow(k)_"beam"| , quad
  hat(x) = hat(y) times hat(q) , $ <eq:frame>

which is right-handed ($hat(x) times hat(y) = hat(q)$) with
$hat(x) prop arrow(k)_perp$. The $pi^0$ variables follow:

$ z = E_(pi^0) / nu , $ <eq:z>
$ arrow(p)_perp = arrow(p)_(pi^0) - (arrow(p)_(pi^0) dot hat(q)) hat(q) ,
  quad p_T^2 = |arrow(p)_perp|^2 , $ <eq:pT2>
$ phi_h = arctan_2 (arrow(p)_perp dot hat(y), space arrow(p)_perp dot hat(x)) . $ <eq:phih>

#result-box(title: "φ_h follows the Trento convention — verified")[
  The frame construction of @eq:frame is not merely _asserted_ to be
  Trento; it was checked numerically. Evaluating the implemented
  $arctan_2$ form against the canonical definition

  $ phi_h = "sign"[(arrow(q) times arrow(k)) dot arrow(P)_h] dot
    arccos[((arrow(q) times arrow(k)) dot (arrow(q) times arrow(P)_h)) /
    (|arrow(q) times arrow(k)| |arrow(q) times arrow(P)_h|)] $

  over random $(e', pi^0)$ configurations at $E_"beam" = 10.53$ GeV
  reproduces it to $5 times 10^(-12)$ degrees for
  $|phi_h| in [1degree, 179degree]$. The convention may be quoted as Trento.

  The frame is also correct analytically, independent of any numerical test:
  $(arrow(q) times arrow(k)) times arrow(q) = |arrow(q)|^2 (arrow(k) - (arrow(k) dot hat(q)) hat(q))$,
  so $hat(x)$ lies along the incoming lepton's transverse component and
  $(hat(x), hat(y), hat(q))$ is right-handed, exactly as Trento requires.

  Note that $phi_h$ is stored in *degrees* throughout the C++ and the
  ntuple; the Python converts to radians before fitting.
]

#important-box(title: "Do not quote an all-φ agreement figure")[
  The restriction to $|phi_h| in [1degree, 179degree]$ above is not
  cosmetic, and an earlier draft of this note was wrong to omit it. It
  quoted a maximum deviation of $8.3 times 10^(-11)$ degrees over
  $20 ,000$ configurations, with no restriction on $phi_h$. That figure
  does not survive scrutiny:

  - it is *seed-dependent* --- over ten seeds at the same sample size it
    ranges from $1.7 times 10^(-11)$ to $9.3 times 10^(-10)$ degrees;
  - and it is *unbounded in sample size*, degrading from
    $3 times 10^(-11)$ at $N = 2 times 10^4$ to
    $5.3 times 10^(-8)$ degrees at $N = 2 times 10^7$.

  The cause is the *reference*, not the implementation. Every worst case
  sits at $phi_h approx 0$ or $plus.minus 180degree$, where $arccos$ is
  ill-conditioned: an error $epsilon$ in $cos phi_h$ becomes
  $tilde sqrt(2 epsilon)$ in $phi_h$. So an all-$phi$ "agreement to $X$"
  number measures the arccos reference and the random seed --- not the code
  under test --- and shrinking $X$ by taking more samples is impossible by
  construction.

  The bound that *does* hold uniformly in $phi_h$ comes from removing
  $arccos$ from the test entirely: construct a hadron with a known $phi_h$
  in the $(hat(x), hat(y), hat(q))$ frame and require the implementation to
  recover it. That round trip agrees to $5 times 10^(-13)$ degrees, uniform
  in $phi_h$, and it bounds the implementation rather than the reference.
]
