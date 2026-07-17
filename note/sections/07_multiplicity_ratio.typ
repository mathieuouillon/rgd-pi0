#import "../template/lib.typ": *

= Multiplicity ratio <sec:multiplicity>

== Definition <sec:ra-definition>

Per kd-tree leaf,

$ R_A = (Y_A \/ N_A^"DIS") / (Y_D \/ N_D^"DIS") , $ <eq:ra>

with $Y$ the background-subtracted $pi^0$ yield in the leaf and $N^"DIS"$
the inclusive @DIS count in the leaf's $(Q^2, x_B)$ box. LD#sub[2] is
always the denominator; CxC, Cu and Sn are three
independent numerator files sharing one kd-tree, which is what makes the
leaf-by-leaf division meaningful.

The @DIS normalisation is integrated over $(Q^2, x_B)$ *only* --- not over
$z$ or $p_T^2$. This is the physically correct choice, and it is verifiable
in the output: all $25$ leaves ($5_z times 5_(p_T^2)$) of a given
$(Q^2, x_B)$ cell carry an identical $N^"DIS"$.

#note-box(title: "What counts as a DIS event")[
  The @DIS normalisation is filled *once per event* that has a good @DIS
  electron --- no photon requirement, no $pi^0$ requirement. It is a
  per-event count, not a per-$pi^0$ count. (A separate `Q2_vs_xB`
  histogram *is* filled per $pi^0$ candidate and must not be used as a
  normalisation; the two differ by the $pi^0$ multiplicity.)
]

== Yield extraction <sec:yield-extraction>

The yield is *not* the integral of a fitted Gaussian. The procedure is:

#figure(
  block(inset: (x: 6pt, y: 4pt))[
    #set text(size: 9.5pt)
    #set par(justify: false)
    + Subtract the sideband-scaled mixed spectrum: $S(m)$ from
      @eq:subtraction.
    + Fit $S(m)$ with a *pure Gaussian* --- no background term, since the
      background is already gone --- over $[0.08, 0.20]$ GeV.
    + *Iteratively restrict to the core*: refit over
      $[mu - 1.5 sigma, mu + 1.5 sigma]$, up to 2 iterations, stopping when
      $mu$ and $sigma$ both change by less than $1%$.
    + Sum the *subtracted histogram* over
      $mu plus.minus 3 sigma$: $Y = sum_(m in plus.minus 3 sigma) S(m)$,
      with $sigma_Y^2 = sum sigma_S^2 (m)$.
  ],
  caption: [Yield extraction. Steps 2--3 exist *only* to locate the
  integration window; the fit amplitude never enters the yield. This is a
  deliberate and defensible choice: it makes the yield insensitive to the
  Gaussian's inability to describe the peak's non-Gaussian tails, at the
  cost of a mild sensitivity to the fitted $mu$ and $sigma$.],
) <fig:yield-procedure>

The Gaussian is $A exp(-(m - mu)^2 \/ 2 sigma^2)$ with bounds
$A in [0, infinity)$, $mu in [0.10, 0.18]$, $sigma in [0.003, 0.05]$ GeV,
seeded at $(max S, 0.135, 0.012)$, fitted with `absolute_sigma=True` so the
errors are not rescaled by $chi^2 \/ "ndf"$.

#important-box(title: "A stale docstring to ignore")[
  The script's own header claims the yield is
  $N = A sigma sqrt(2 pi) \/ w_"bin"$ with the error from the full
  covariance matrix. *That is not what produces $R_A$.* The analytic
  integral and its covariance-propagated error are implemented, but used
  *only* in a pedagogical demonstration plot. The published yield is the
  $plus.minus 3 sigma$ window sum of @fig:yield-procedure.
]

Two quality gates apply: at least 5 raw same-event entries, and a
subtracted peak amplitude of at least 1 count. In the actual production
*all $1450$ leaves survive for all three targets* --- zero rejections. Fit
quality is good: median $chi^2 \/ "ndf"$ between $1.00$ and $1.05$, with
fewer than $0.7%$ of leaves above 2.

== Normalisation integration <sec:normalisation>

Because $N^"DIS"$ is histogrammed on the *configuration* product grid while
the yields live on kd-tree leaves (@sec:aux-grids), the denominator must be
integrated over each leaf's $(Q^2, x_B)$ box. The preferred path reads the
two-dimensional $N^"DIS"(Q^2, x_B)$ histogram and weights each cell by its
*separable* overlap fraction $f_(Q^2) times f_(x_B)$, assuming uniform
density within a cell.

#warning-box(title: "The factorized fallback is not safe")[
  If the 2D histogram is unavailable the code falls back to

  $ N^"DIS" (Q^2, x_B) approx (N(Q^2) dot N(x_B)) / N_"inc" , $

  i.e.\ it assumes $rho(Q^2, x_B) approx rho(Q^2) rho(x_B)$. *That
  factorization is badly wrong for @DIS*, where $Q^2$ and $x_B$ are
  strongly correlated by $Q^2 = 2 M_p nu x_B$ and the acceptance carves a
  diagonal band out of the plane.

  The script prints which branch it took. *That log line must be checked
  and recorded*: the 2D histogram is genuinely booked and filled by the
  C++, so the preferred branch should be taken --- but the note cannot
  assert which branch produced these numbers without the log.

  A related subtlety: the 2D normalisation histogram has *flow bins*, so
  events with $x_B < 0.1$ or $x_B > 0.7$ land outside the grid. The in-grid
  sum is therefore not the total @DIS count, and the leaf integration is
  only meaningful for leaves inside the product grid's span.
]

== Uncertainties <sec:ra-uncertainties>

The statistical uncertainty is propagated as

$ sigma_R / R = sqrt((sigma_(Y_D) / Y_D)^2 + (sigma_(Y_A) / Y_A)^2
  + 1 / N_D^"DIS" + 1 / N_A^"DIS") , $ <eq:ra-error>

treating the @DIS counts as Poisson. The measured budget is dominated
entirely by the yields:

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, right, right),
    table.header([*Term*], [*Median*], [*Share of rel. variance*]),
    [$(sigma_(Y_D) \/ Y_D)^2$], [$1.19 times 10^(-4)$], [$36.7%$],
    [$(sigma_(Y_A) \/ Y_A)^2$], [$1.99 times 10^(-4)$], [$61.1%$],
    [$1 \/ N_D^"DIS"$], [$2.50 times 10^(-6)$], [$0.8%$],
    [$1 \/ N_A^"DIS"$], [$4.22 times 10^(-6)$], [$1.3%$],
  ),
  caption: [Statistical error budget for CxC, medians over
  the $1450$ leaves. Median $sigma_R \/ R = 1.80%$. The @DIS normalisation
  terms are negligible --- the measurement is limited by the $pi^0$ yield
  extraction, which is in turn limited by the background subtraction.],
) <tab:error-budget>

#important-box(title: "Three known omissions in the error model")[
  + The uncertainty on the sideband scale $alpha$ is not propagated
    (@sec:subtraction).
  + The $plus.minus 3 sigma$ window sum ignores the uncertainty on the
    fitted $mu$ and $sigma$ that *define* the window.
  + *$R_"CxC"$, $R_"Cu"$ and $R_"Sn"$ are correlated* --- they share the
    LD#sub[2] numerator yield and @DIS count --- but the code treats them as
    independent. This does not affect any single $R_A$, but it *does* matter
    for any fit of the $A$-dependence, which is precisely what
    @sec:results-ptb does with $Delta chevron.l p_T^2 chevron.r$. An
    $A$-dependence fit using these errors will understate its own
    uncertainty.
]

== The no-subtraction cross-check <sec:no-subtraction>

A second variant is produced for every leaf: a Gaussian-plus-linear fit to
the *raw* same-event spectrum, from which $mu$ and $sigma$ are taken and
the *raw counts* summed over $plus.minus 3 sigma$.

The fitted linear background is *discarded*, not subtracted. This variant
is therefore not an alternative background subtraction --- it is a
*background-dilution bound*: a measurement of what $R_A$ would be with no
subtraction at all, which is pulled toward unity by the unsubtracted
combinatorial background.

It is plotted alongside the nominal but *never differenced into a
systematic band*, and never written to CSV. Turning the nominal-versus-raw
difference into a background systematic is the most obvious immediate
improvement available: the two variants already exist for every leaf.
