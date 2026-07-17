#import "../template/lib.typ": *

= Current results <sec:results>

#important-box(title: "How to read this section")[
  Every number below is *statistical only* and carries the caveats of
  @sec:systematics --- no acceptance, no bin migration, no radiative
  correction, no evaluated systematic. The @BSA additionally depends on a
  placeholder polarization (@sec:bsa-polarization), and the quoted
  kinematic positions of the widest bins are not where the data sits
  (@sec:binning-caveat).

  All values were recomputed directly from the result files for this note
  rather than transcribed, and the fit-quality and error-budget statistics
  quoted throughout were reproduced independently.

  *These are internal, provisional numbers.* They are presented because
  they are informative about the state of the analysis, not because they
  are ready to circulate.
]

== Multiplicity ratio <sec:results-ra>

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    table.header(
      [*$z$ bin*], [*$N_"leaves"$*], [*$R_"CxC"$*], [*$R_"Cu"$*],
      [*$R_"Sn"$*],
    ),
    [$(0.00, 0.15]$], [330], [$1.0396 plus.minus 0.0014$],
      [$1.0996 plus.minus 0.0015$], [$1.1545 plus.minus 0.0014$],
    [$(0.15, 0.25]$], [395], [$1.0145 plus.minus 0.0010$],
      [$1.0248 plus.minus 0.0010$], [$1.0408 plus.minus 0.0010$],
    [$(0.25, 0.35]$], [285], [$1.0081 plus.minus 0.0010$],
      [$0.9955 plus.minus 0.0011$], [$0.9910 plus.minus 0.0010$],
    [$(0.35, 0.45]$], [115], [$0.9979 plus.minus 0.0015$],
      [$0.9733 plus.minus 0.0016$], [$0.9609 plus.minus 0.0015$],
    [$(0.45, 0.55]$], [35], [$0.9908 plus.minus 0.0026$],
      [$0.9633 plus.minus 0.0028$], [$0.9420 plus.minus 0.0026$],
    [$(0.55, 0.65]$], [20], [$1.0014 plus.minus 0.0035$],
      [$0.9816 plus.minus 0.0035$], [$0.9411 plus.minus 0.0032$],
    [$(0.65, 0.80]$], [270], [$0.9911 plus.minus 0.0009$],
      [$0.9385 plus.minus 0.0009$], [$0.8973 plus.minus 0.0008$],
  ),
  caption: [$R_A$ versus $z$, inverse-variance weighted over kd-tree
  leaves, statistical errors only. Binning is on the *reported* $z$ column,
  i.e.\ the geometric box midpoint. *The last row is not a measurement at
  $z approx 0.7$*: it is the top $z$ box, spanning $z in [0.37, 1.0]$ and
  dominated by $z approx 0.4-0.5$ (@sec:binning-caveat). The sparse
  $(0.45, 0.65]$ rows are the sampling of interior boxes near the top-box
  boundary and should not be over-interpreted.],
) <tab:R-vs-z>

The qualitative picture is textbook nuclear hadronization, and all three
signatures come out right:

/ Attenuation grows with $z$: $R$ falls monotonically for Cu and Sn, from
  above unity at low $z$ to a clear suppression at high $z$.
/ Attenuation orders with $A$: at every $z$ above $0.25$,
  $R_"CxC" > R_"Cu" > R_"Sn"$ --- the ordering expected from the growing
  path length in nuclear matter, and *not* an ordering the analysis was
  tuned to produce.
/ Low-$z$ enhancement: $R > 1$ below $z approx 0.25$, growing with $A$
  (Sn reaches $1.15$), consistent with energy removed from the leading
  hadron reappearing as additional soft production.

Full ranges: $R_"Cu" in [0.818, 1.352]$, $R_"Sn" in [0.776, 1.511]$,
$R_"CxC" in [0.897, 1.188]$ over all $1450$ leaves.

#wide-figure(
  "../figures/R_vs_z_per_Q2xB.pdf",
  [$R_A$ versus $z$ for one $(Q^2, x_B)$ cell, one of the 58 pages of the
  full result --- the page count is itself a check on @tab:trees, since
  $58 times 25 = 1450$. The three targets separate in the expected order at
  large $z$. *The rightmost point of each curve is the top $z$ box*, plotted
  at its geometric midpoint rather than at its true
  $chevron.l z chevron.r$; the visible gap between the fourth and fifth
  points is the signature of that artefact (@sec:binning-caveat), not a
  feature of the data.],
  <fig:R-vs-z>,
  width: 90%,
  page: 1,
)

== Transverse-momentum broadening <sec:results-ptb>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, right, right, right),
    table.header([*Target*], [*$A$*],
                 [*$Delta chevron.l p_T^2 chevron.r$ (GeV#super[2])*],
                 [*Unweighted mean $plus.minus$ s.d.*]),
    [CxC], [12.0], [$0.001679 plus.minus 0.000034$],
      [$0.0052 plus.minus 0.0046$],
    [Cu], [63.5], [$0.003021 plus.minus 0.000035$],
      [$0.0101 plus.minus 0.0088$],
    [Sn], [118.7], [$0.004662 plus.minus 0.000033$],
      [$0.0147 plus.minus 0.0122$],
  ),
  caption: [$Delta chevron.l p_T^2 chevron.r$ over the $320$ leaves of the 3D
  tree. The broadening is positive and monotonic in $A$, with essentially
  no negative leaves ($0.31%$ for CxC, none for Cu or Sn).
  $chevron.l p_T^2 chevron.r_(D)$ has mean $0.128$ GeV#super[2].],
) <tab:ptb>

#important-box(title: "Weighted and unweighted means differ by ~3×")[
  @tab:ptb deliberately quotes both. The inverse-variance weighted mean is
  dominated by high-statistics, low-$Delta$ leaves; the unweighted mean
  treats every leaf equally. They differ by a factor $approx 3$, which is not
  a small effect --- it says the broadening varies strongly across the
  kinematic range, as it should ($Delta chevron.l p_T^2 chevron.r$ grows with
  path length and with $z$). *Any quoted single number must state which
  average it is.* Neither is "the" broadening; the physics is in the
  differential distribution.
]

#subfig2(
  (
    ("../figures/pT_broadening_vs_z_per_Q2xB.pdf",
      [$Delta chevron.l p_T^2 chevron.r$ vs $z$, one $(Q^2, x_B)$ cell],
      <fig:ptb-vs-z>),
    ("../figures/pT_broadening_vs_A.pdf",
      [$chevron.l p_T^2 chevron.r$ vs $A$, one $z$ bin],
      <fig:ptb-vs-A>),
  ),
  [Transverse-momentum broadening. Left: one of the 64 $(Q^2, x_B)$ pages
  ($64 times 5 = 320$ leaves, confirming @tab:trees). Right: the
  $A$-dependence that @eq:A-scaling summarises. The $y$-axis of the left
  panel is clipped at $0.05$ GeV#super[2], which cuts three Sn and one
  CxC leaf --- worth widening before circulation.],
  <fig:ptb>,
)

Fitting the weighted values against $A$ gives

$ Delta chevron.l p_T^2 chevron.r prop A^(0.426) , $ <eq:A-scaling>

to be compared with the $A^(1\/3)$ of @eq:A13, expected if $hat(q)$ is a
property of the medium and the path length simply tracks the nuclear radius.
The measured exponent is notably *above* $1\/3$.

Taken at face value that would mean the broadening grows with $A$ *faster*
than the nucleus does --- so either $hat(q)$ itself rises with $A$ (denser
matter in heavier nuclei, beyond the constant-density picture behind
@tab:lengths), or the effective path length grows faster than $R_A$, or the
light target is losing broadening for a reason the heavy ones are not. It is
the kind of deviation the measurement exists to find. It is also, at present,
not established:

#warning-box(title: "Do not over-interpret the exponent")[
  @eq:A-scaling is a three-point fit using statistical errors that are
  known to be *underestimated in a correlated way*: the three targets share
  the LD#sub[2] reference, and that correlation is not modelled
  (@sec:ra-uncertainties). The exponent is also computed from the weighted
  mean, which per the box above is only one of two defensible averages ---
  the unweighted values give a different exponent.

  Layered on top, the $p_T^2$ moments are accumulated over a
  background-dominated sample (@sec:ptb-caveat), which suppresses
  $Delta chevron.l p_T^2 chevron.r$ toward zero by a factor that need not be
  the same for every target. The deviation from $A^(1\/3)$ is *interesting
  and worth pursuing*, but it is not yet a result.
]

== Beam-spin asymmetry <sec:results-bsa>

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, right, right),
    table.header([*Target*], [*$A_"LU"^(sin phi_h)$*], [*$sum N$*]),
    [LD#sub[2]], [$0.01378 plus.minus 0.00040$], [$15 ,726 ,897$],
    [CxC], [$0.01336 plus.minus 0.00050$], [$9 ,718 ,874$],
    [Cu], [$0.01242 plus.minus 0.00054$], [$8 ,517 ,384$],
    [Sn], [$0.01201 plus.minus 0.00048$], [$10 ,586 ,165$],
  ),
  caption: [$A_"LU"^(sin phi_h)$, inverse-variance weighted over the $235$
  clean leaves, at the placeholder $P = 0.85$. All four targets agree
  within $approx 1.5 sigma$: *there is no significant nuclear dependence* at
  this precision. The central values fall monotonically with $A$, but not
  significantly so. Values scale as $1\/P$ and will move when the measured
  polarization is applied.],
) <tab:bsa>

#warning-box(title: "Two railed fits destroy the naive average")[
  Two of the $237$ leaves --- ids *72* and *170*, both at low
  $Q^2 approx 1.19$, $x_B approx 0.10$ --- have fits pinned at the
  parameter bound $|A| = 0.5$ with
  $sigma_A = 3.86 times 10^(-6)$ and $chi^2 \/ "ndf" approx 10^(11)$.

  Inverse-variance weighting hands them weight $approx 6.7 times 10^(10)$,
  around $10^(10)$ times that of a real leaf. The unfiltered weighted
  averages are consequently

  #align(center)[
    #text(size: 9.5pt)[
      $A_"LU"("CxC") = +0.49997$, quad $A_"LU"("Cu") = -0.49997$
    ]
  ]

  --- i.e.\ *the railed leaves alone*, with the other 236 contributing
  nothing. @tab:bsa excludes them, giving the physical $approx 0.012 - 0.014$.

  *Never compute an average from this file without filtering.* The result
  file as written contains no quality flag, so nothing warns a downstream
  user. A fit-quality mask (reject $|A|$ at the bound, $sigma_A < 10^(-4)$,
  or $chi^2\/"ndf" > 10$) belongs in the producing script, not in each
  consumer.
]

Binned in $z$, the LD#sub[2] asymmetry rises mildly from $0.0106$ below
$z = 0.15$ to $0.0159$ near $z approx 0.3$ --- the expected trend. Fit
quality is good: median $chi^2\/"ndf"$ of $0.95 - 0.98$ with a maximum of
$2.7$ across the clean leaves.

#wide-figure(
  "../figures/ALU_grid_perQ2.pdf",
  [$A_"LU"^(sin phi_h)$ across the $(x_B, p_T^2)$ display grid for one
  $Q^2$ slice. The four targets overlap within their uncertainties
  throughout --- the visual form of the statement in @tab:bsa that no
  nuclear dependence is resolved. Values shown use the *placeholder*
  $P = 0.85$ and scale as $1\/P$ (@sec:bsa-polarization).],
  <fig:alu-grid>,
  width: 90%,
  page: 1,
)