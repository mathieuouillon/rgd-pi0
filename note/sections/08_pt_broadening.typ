#import "../template/lib.typ": *

= Transverse-momentum broadening <sec:pt-broadening>

== Definition <sec:ptb-definition>

$ Delta chevron.l p_T^2 chevron.r_A = chevron.l p_T^2 chevron.r_A - chevron.l p_T^2 chevron.r_D $ <eq:broadening>

evaluated at matched bin index of the *three-dimensional* $(Q^2, x_B, z)$
grid (@tab:trees). The grid is 3D precisely so that $p_T^2$ is the
*observable* --- averaged within each bin --- rather than a binning
dimension: one cannot take the mean of a quantity one has binned on.

== Accumulator means <sec:ptb-accumulators>

$chevron.l p_T^2 chevron.r$ is a plain accumulator mean --- *no fit, no
exponential model, no truncation, no $p_T^2$ range cut*. Stage B fills the
`ptb3d` tree with, for each 3D bin $k$ *and each $m_(gamma gamma)$ bin*,

$ H_"counts" [k] = N_k , quad
  H_(sum p_T^2) [k] = sum_(i in k) p_(T,i)^2 , quad
  H_(sum p_T^4) [k] = sum_(i in k) p_(T,i)^4 , $ <eq:ptb-accum>

with a compensated (Neumaier) accumulator, so the sums are bit-identical for
any thread count. The Python forms

$ chevron.l p_T^2 chevron.r = (sum p_T^2) / N , quad
  sigma_(chevron.l p_T^2 chevron.r)^2
    = (chevron.l p_T^4 chevron.r - chevron.l p_T^2 chevron.r^2) / N . $ <eq:ptb-mean>

The second relation is just the standard error of the mean of the $p_T^2$
distribution --- which is the whole reason $sum p_T^4$ is accumulated at
all. Errors combine in quadrature,
$sigma_Delta = sqrt(sigma_A^2 + sigma_D^2)$, treating target and reference
as independent (valid --- they are separate files). The 3D grid has $280$
bins; in the diagnostic run of @sec:results, $197$ carry data for all three
targets.

#note-box(title: "Why this observable is technically the cleanest")[
  $Delta chevron.l p_T^2 chevron.r$ needs no yield fit, no window sum, and no
  normalisation. It is a difference of two means, so the @DIS normalisation,
  the luminosity and the beam charge all cancel identically rather than
  approximately, and its statistical error is analytic. In practice it is the
  observable that comes out cleanest (@sec:results-ptb).
]

== Background handling <sec:ptb-caveat>

The $p_T^2$ moments are accumulated from the $pi^0$ *candidate* collection ---
the greedy-paired, mass-windowed sample of @sec:pairing, which spans the full
`pairing.mass_window_gev` of `config/cuts.json` and so is mostly combinatorial. The
key improvement over a naive average is that the moments are binned in
$m_(gamma gamma)$ (@eq:ptb-accum), so the *same* sideband subtraction used for
the yields (@sec:subtraction) is applied to them:
$chevron.l p_T^2 chevron.r$ is reported over the sideband-subtracted signal
region, not the raw window.

#note-box(title: "Residual dilution")[
  Sideband subtraction removes the flat combinatorial background under the
  peak, but where the signal fraction is low the subtracted moment still
  carries a statistical penalty and some residual shape sensitivity. The
  *difference* @eq:broadening cancels whatever remains to the extent that it is
  target-independent, and the measured $Delta chevron.l p_T^2 chevron.r$ is
  positive, monotonic in $A$ and of a sensible magnitude (@sec:results-ptb).
  What is *not* yet done is a systematic on the residual: varying the sideband
  and the window into a band (@tab:systematics).
]
