#import "../template/lib.typ": *

= Beam-spin asymmetry <sec:bsa>

== Definition <sec:bsa-definition>

With a longitudinally polarised beam on an unpolarised target, the
helicity asymmetry in a given $(phi_h, "bin")$ bin is

$ A_"LU" (phi_h) = 1/P dot (N^+ - N^-) / (N^+ + N^-) , $ <eq:alu>

$ sigma_(A_"LU") = 1/P sqrt((1 - (P dot A_"LU")^2) / (N^+ + N^-)) . $ <eq:alu-err>

Writing $r = (N^+ - N^-)\/(N^+ + N^-) = P A_"LU"$, the variance term is
exactly the binomial $(1 - r^2)\/N$ --- correct, and well defined for all
$|r| <= 1$. No charge or luminosity normalisation enters, and none is
needed: the beam charge cancels identically in the ratio provided the
helicity flips fast compared with any drift.

#result-box(title: "A_LU is acceptance-safe by construction")[
  @eq:alu is a ratio of yields taken *within a single $phi_h$ bin*, so any
  helicity-independent acceptance cancels *exactly* --- not approximately.
  This is what makes the @BSA the most robust observable in this note
  despite the analysis having no acceptance correction at all. The contrast
  with a raw $chevron.l cos phi_h chevron.r$-type moment is instructive:
  such a moment is taken *across* $phi_h$ bins, so the acceptance shape does
  not divide out and dominates the result.
]

== Helicity <sec:helicity>

The helicity is read from `REC::Event`'s *HWP-corrected* accessor and
propagated unchanged through the skim ($plus.minus 1$, stored as an integer
in `PI0::event`). Events with `helicity == 0` (undefined) are dropped.

#important-box(title: "HWP handling is inherited, not verified")[
  The chain uses the HWP-corrected `helicity` field throughout and never
  touches `helicityRaw`. The @HWP state is therefore handled *entirely
  upstream* by the CLAS12 cooking and @RCDB --- there is no per-run-period
  @HWP bookkeeping at any level of this analysis, and no HWP-in versus
  HWP-out closure test.

  Likewise the *sign* of $A_"LU"^(sin phi_h)$ rests on two inherited
  conventions: the cooking's helicity sign and the Trento $phi_h$ sign.
  The latter is verified (@sec:sidis-kinematics); the former is not
  validated anywhere in this analysis. Before publishing a *signed*
  asymmetry, the sign should be confirmed against a known CLAS12
  measurement on a common channel.

  No false-asymmetry or charge-asymmetry systematic has been evaluated.
]

== Fit model <sec:bsa-fit>

In each bin the asymmetry is binned in 12 equal $phi_h$ bins over
$[-180degree, 180degree]$ (matching the C++ histogram exactly) and fitted
with a *single-parameter* model:

$ f(phi_h) = A sin phi_h , quad A equiv A_"LU"^(sin phi_h) , $ <eq:bsa-model>

by $chi^2$ minimisation on the asymmetry points, with `absolute_sigma=True`
and $A$ bounded to $[-0.5, 0.5]$. With all 12 bins surviving,
$"ndf" = 11$.

#note-box(title: "Why the three-parameter ratio form was abandoned")[
  The canonical one-photon-exchange expression is

  $ A_"LU" (phi_h) = (A_"LU"^(sin phi_h) sin phi_h) /
    (1 + A_"LU"^(cos phi_h) cos phi_h + A_"LU"^(cos 2 phi_h) cos 2 phi_h) , $

  and an earlier iteration of this analysis fitted it. It was *deliberately
  abandoned*: at these per-bin statistics the denominator develops pole
  pathologies and, in the author's words recorded in the source, the fit
  landscape becomes meaningless. The unpolarised $cos phi_h$ and
  $cos 2 phi_h$ moments properly belong to $A_"UU"$ and would have to be
  measured separately.

  This is a sound decision, but it has a consequence that must be stated:
  the quoted $A_"LU"^(sin phi_h)$ carries *no correction for the
  $cos phi_h$ dilution of the denominator*. It is the $sin phi_h$ moment of
  the asymmetry, not of the structure-function ratio.
]

An unbinned maximum-likelihood alternative is documented in
@sec:appendix-mlm. It is not currently used.

== Signal extraction — there is none <sec:bsa-signal>

#warning-box(title: "The BSA applies no background subtraction")[
  The per-$(phi_h, "bin")$ yields $N^plus.minus$ are *raw counts* inside a
  fixed window $0.110 < m_(gamma gamma) < 0.160$ GeV
  ($135 plus.minus 25$ MeV). There is no mass fit and no sideband
  subtraction --- unlike the multiplicity analysis, which does both.

  If the combinatorial background is unpolarised, as expected, it dilutes
  the asymmetry *multiplicatively*:

  $ A_"LU"^"meas" = A_"LU"^"true" dot S / (S + B) . $ <eq:dilution>

  No purity or dilution factor is computed or applied, so *every quoted
  $A_"LU"$ is a lower bound on its true magnitude*, by a factor that is
  neither measured nor necessarily target-independent.

  The machinery to fix this already exists and is unused: the same- and
  mixed-event $m_(gamma gamma)$ histograms are booked and filled in the
  same C++ pass, binned on the same tree. The signal fraction $S\/(S+B)$
  could be read off per bin from the fits the multiplicity analysis
  already performs.
]

Quality gates: at least $200$ events per bin, at least $5$ per $phi_h$
bin, at least $6$ surviving $phi_h$ bins, and a converged fit. Stage B writes
the sparse `bsa` tree on the 4D grid; the asymmetry itself is not extracted in
the diagnostic production of @sec:results, because it has no measured
polarization (below).

#wide-figure(
  "../figures/phi_h_acceptance.pdf",
  [Raw $phi_h$ yield $N^+ + N^-$ for a representative bin --- the honest
  acceptance diagnostic. The shape is *not* physics: it is the CLAS12
  azimuthal acceptance, with the sector structure and the gaps between
  sectors plainly visible. This figure is the argument of
  @sec:bsa-definition made visually: $A_"LU"$ is formed *within* each of
  these $phi_h$ bins, so this entire shape divides out exactly, whereas a
  moment taken *across* these $phi_h$ bins would instead be dominated by
  this shape.],
  <fig:phih-acceptance>,
  width: 88%,
  page: 1,
)

== Beam polarization <sec:bsa-polarization>

#warning-box(title: "A_LU has no polarization, by design")[
  `pi0.bsa` *refuses to run* without a beam polarization: it is mandatory as
  `--polarization P` or `config/cuts.json` `/bsa/polarization/value`, and there
  is deliberately no default --- the old code's self-declared $P = 0.85$ was
  removed precisely because a placeholder that silently scales every $A_"LU"$
  is worse than a hard stop. Since $A_"LU" prop 1\/P$, the polarization is a
  fully-correlated scale on every asymmetry, and its uncertainty
  ($sigma_A^"pol" = A dot sigma_P \/ P$) is propagated once $P plus.minus
  sigma_P$ is supplied.

  Obtaining the measured @RGD Møller polarization, per run period, and
  recording it in `config/cuts.json` is the single blocking item for this
  observable. Until then the `bsa` tree is written and waiting but $A_"LU"$ is
  not extracted (@sec:results-bsa).
]
