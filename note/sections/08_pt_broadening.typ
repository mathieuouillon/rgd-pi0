#import "../template/lib.typ": *

= Transverse-momentum broadening <sec:pt-broadening>

== Definition <sec:ptb-definition>

$ Delta chevron.l p_T^2 chevron.r_A = chevron.l p_T^2 chevron.r_A - chevron.l p_T^2 chevron.r_D $ <eq:broadening>

evaluated at matched leaf index of the *three-dimensional*
$(Q^2, x_B, z)$ tree (@tab:trees). The tree is 3D precisely so that
$p_T^2$ is *integrated* within each leaf rather than binned --- one cannot
take the mean of a quantity one has binned on.

== Accumulator means <sec:ptb-accumulators>

$chevron.l p_T^2 chevron.r$ is a plain accumulator mean --- *no fit, no
exponential model, no truncation, no $p_T^2$ range cut*. The C++ fills
three histograms with one bin per leaf,

$ H_"counts" [k] = N_k , quad
  H_(sum p_T^2) [k] = sum_(i in k) p_(T,i)^2 , quad
  H_(sum p_T^4) [k] = sum_(i in k) p_(T,i)^4 , $ <eq:ptb-accum>

from which the Python forms

$ chevron.l p_T^2 chevron.r = (sum p_T^2) / N , quad
  sigma_(chevron.l p_T^2 chevron.r)^2
    = (chevron.l p_T^4 chevron.r - chevron.l p_T^2 chevron.r^2) / N . $ <eq:ptb-mean>

The second relation is just the standard error of the mean of the $p_T^2$
distribution --- which is the whole reason $sum p_T^4$ is accumulated at
all. Errors combine in quadrature,
$sigma_Delta = sqrt(sigma_A^2 + sigma_D^2)$, treating target and reference
as independent (valid --- they are separate files). Leaves with fewer than
$20$ entries are dropped; in practice all $320$ survive.

#note-box(title: "Why this observable is technically the cleanest")[
  $Delta chevron.l p_T^2 chevron.r$ needs no fit, no window, and no
  normalisation. It is a difference of two means, so the @DIS
  normalisation, the luminosity and the beam charge all cancel identically
  rather than approximately. Its statistical error is analytic and exact.
  Everything that can go wrong with it is therefore in the *sample
  definition* --- which is where it does go wrong (@sec:ptb-caveat).
]

== The mass-window problem <sec:ptb-caveat>

#warning-box(title: "⟨p_T²⟩ is accumulated over a largely combinatorial sample")[
  The accumulators of @eq:ptb-accum are filled from the $pi^0$ *candidate*
  collection --- the greedy-paired, mass-windowed sample of
  @sec:pairing. With the shipped configuration that window is
  $plus.minus 200$ MeV, i.e.\ $m_(gamma gamma) < 0.335$ GeV with no lower
  bound.

  Unlike the multiplicity ratio, *there is no offline background
  subtraction here.* The multiplicity analysis recovers the signal by
  fitting and subtracting the $m_(gamma gamma)$ spectrum; the
  $p_T$-broadening analysis simply averages $p_T^2$ over everything in the
  collection. So $chevron.l p_T^2 chevron.r$ as computed is

  $ chevron.l p_T^2 chevron.r^"meas"
    = f_S chevron.l p_T^2 chevron.r^(pi^0) + (1 - f_S) chevron.l p_T^2 chevron.r^"bkg" $

  with $f_S$ the signal fraction --- and $f_S$ is *not* close to 1.

  The script's own docstring asserts the accumulators are $pi^0$-selected
  at $plus.minus 30$ MeV. *That assertion does not match the shipped
  configuration.*

  Mitigating considerations, which is why the numbers are still presented:
  the *difference* @eq:broadening cancels the background contribution to
  the extent that $chevron.l p_T^2 chevron.r^"bkg"$ and $f_S$ are
  target-independent, and the observed $Delta chevron.l p_T^2 chevron.r$ is
  positive, monotonic in $A$ and of a physically sensible magnitude
  (@sec:results-ptb) --- which would be a surprising accident if the
  measurement were pure background. But the dilution certainly *suppresses*
  the magnitude toward zero, and it is not corrected.

  *Resolution required.* Either confirm the production mass window from the
  farm-side configuration, or re-run the accumulators inside a true
  $plus.minus 3 sigma$ window, or --- best --- accumulate $sum p_T^2$ and
  $sum p_T^4$ *per $m_(gamma gamma)$ bin* so that the same sideband
  subtraction used for the yields can be applied to the moments. The third
  option is a small C++ change and would make this observable as clean as
  its statistics deserve.
]

== The orphaned cache <sec:ptb-cache>

A $3.6$ GB tree of `.npz` / `.npy` caches sits under `results/pT2_cache/`,
named
`pi0_multiplicity_<target>_<polarity>__mgg_<lo>_<hi>[__v<N>]`, encoding the
ROOT stem, the mass window and a schema version.

*Nothing reads it.* No script in either repository references
`pT2_cache`; the current $p_T$-broadening code reads the C++ histograms
directly and has no cache path at all. The caches are a fossil of a
superseded event-loop implementation which binned on the product grid
(their schema stores `pT`, not `pT2`, and `iz_mult`, not leaf indices), and
their stem `..._Outbending` no longer even matches today's `..._OB.root`,
so they could never hit.

Two facts worth recording before deletion. First, the generation bump is
explicable exactly: the LD#sub[2] v3 cache holds $3 ,644 ,002$ entries
against the older `.npz`'s $4 ,093 ,856$, a difference of $449 ,854$,
which is precisely the addition of the in-grid requirement
$i_(Q^2), i_(x_B), i_(z"mult") >= 0$. Second, the cache key contains *no
file modification time or size*, so an updated ROOT file would *not*
invalidate it --- the stale LD#sub[2] cache holds $3.6$M $pi^0$ where the
live file yields $15.7$M, a $4.3 times$ gap. Harmless today because the
caches are dead; silently wrong the moment anyone re-enables them. A
further $1.6$ GB of exact duplicates sits in the @BSA cache from an
`Outbending` $arrow.r$ `OB` rename, for $approx 5.2$ GB reclaimable in total.
