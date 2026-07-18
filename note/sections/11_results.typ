#import "../template/lib.typ": *

= Current results <sec:results>

#important-box(title: "First end-to-end run of the rewritten pipeline — diagnostic only")[
  Every number in this section is the *first* output of the rewritten
  Stage A $arrow$ `make_grid` $arrow$ Stage B $arrow$ Stage C chain run end to
  end on data. It is *diagnostic only* --- three independent blockers forbid
  quoting any of it, and each is stamped into the result files:

  - *Photon fallback.* No GBT model covers @RGD, so photons were scored by the
    @RGA inbending model (`gbt.fallback_used = TRUE`). The photon sample is not
    the one the cuts describe.
  - *Truncated statistics.* Stage A read only the first $2 times 10^6$ events
    of each of the $95$ SIDIS-train files, so $N_"DIS"$ --- the normalisation
    denominator --- is a prefix of each run and every yield carries
    `--allow-truncated-inputs`.
  - *No systematics.* Errors are statistical only and correlated through the
    shared LD#sub[2] reference (@sec:ra-uncertainties); no acceptance, bin
    migration, radiative or evaluated systematic is applied (@sec:systematics).

  The sample is $approx 5.5 times 10^6$ @DIS events and $2.5 times 10^6$
  $pi^0$ candidates over the factorised $8 times 7$ $(Q^2, x_B) times 5 times
  5$ $(z, p_T^2)$ grid (`provenance_hash 2acf618294a6c3b0`,
  `cuts.sha256 801ba433…`). It is shown because it exercises the whole chain
  and is *informative about the state of the analysis*, not because it is
  ready to circulate.
]

The broadening is clean and $A$-ordered; the multiplicity ratio shows
attenuation in the inverse-variance weighted mean but not the median; the @BSA
is not yet extracted, as it needs a measured polarization
(@sec:bsa-polarization).

== Transverse-momentum broadening <sec:results-ptb>

The cleanest result, and the one that most directly validates the chain.

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, right, right, right),
    table.header([*Target*], [*$A$*],
                 [*$Delta chevron.l p_T^2 chevron.r$ (GeV#super[2])*],
                 [*$N_"bins"$*]),
    [C (CxC)], [12.0],  [$0.00160 plus.minus 0.00017$], [197],
    [Cu],      [63.5],  [$0.00277 plus.minus 0.00017$], [197],
    [Sn],      [118.7], [$0.00520 plus.minus 0.00016$], [197],
  ),
  caption: [$Delta chevron.l p_T^2 chevron.r$, inverse-variance weighted over
  the $197$ surviving 3D $(Q^2, x_B, z)$ bins per target, statistical errors
  only. Positive, monotonic in $A$, and each more than $4 sigma$ from zero.],
) <tab:ptb>

The broadening is positive, grows with $A$, and --- an external check --- agrees
to within $10-15%$ with an earlier, independent analysis of the same @RGD data
($0.00168$, $0.00302$, $0.00466$ GeV#super[2] for C, Cu, Sn). Two independent
implementations landing on the same values, from a different binning and
different statistics, is strong evidence the pairing $arrow$ donor-pool
$arrow$ sideband-subtraction chain is sound.

#subfig2(
  (
    ("../figures/results_dpt2_vs_z.pdf",
      [$Delta chevron.l p_T^2 chevron.r$ vs $chevron.l z chevron.r$],
      <fig:ptb-vs-z-new>),
    ("../figures/results_dpt2_vs_A.pdf",
      [$Delta chevron.l p_T^2 chevron.r$ vs $A$, with the power-law fit],
      <fig:ptb-vs-A-new>),
  ),
  [Transverse-momentum broadening (diagnostic). Points are inverse-variance
  weighted over the 3D bins; the $z$ axis uses the sideband-subtracted
  $chevron.l z chevron.r$. The right-panel line is @eq:A-scaling.],
  <fig:ptb-new>,
)

A power-law fit of the weighted values against $A$ gives

$ Delta chevron.l p_T^2 chevron.r prop A^(0.475) , $ <eq:A-scaling>

*above* the $A^(1\/3)$ of @eq:A13 expected if $hat(q)$ is a fixed property of
the medium, and consistent with the earlier analysis's $A^(0.43)$.

#warning-box(title: "Do not over-interpret the exponent")[
  @eq:A-scaling is a three-point fit on statistical errors that are
  *underestimated in a correlated way*: the three targets share the LD#sub[2]
  reference and that correlation is not modelled (@sec:ra-uncertainties). The
  $p_T^2$ moments are accumulated over a background-dominated sample
  (@sec:ptb-caveat), which pulls $Delta chevron.l p_T^2 chevron.r$ toward zero
  by a factor that need not be the same for every target. That two independent
  analyses land on the same exponent is encouraging; it is not a measurement.
]

== Multiplicity ratio <sec:results-ra>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, right, right, right),
    table.header([*$chevron.l z chevron.r$*],
                 [*$R_"C"$*], [*$R_"Cu"$*], [*$R_"Sn"$*]),
    [$0.11$],  [$0.940 plus.minus 0.007$], [$1.002 plus.minus 0.008$],
      [$1.077 plus.minus 0.008$],
    [$0.17$],  [$0.948 plus.minus 0.006$], [$0.828 plus.minus 0.006$],
      [$0.969 plus.minus 0.006$],
    [$0.23$],  [$0.929 plus.minus 0.005$], [$0.744 plus.minus 0.005$],
      [$0.925 plus.minus 0.005$],
    [$0.31$],  [$0.916 plus.minus 0.005$], [$0.852 plus.minus 0.005$],
      [$0.865 plus.minus 0.005$],
    [$0.53$],  [$0.912 plus.minus 0.004$], [$0.811 plus.minus 0.004$],
      [$0.674 plus.minus 0.004$],
    table.hline(),
    [all $z$], [$0.925 plus.minus 0.002$], [$0.826 plus.minus 0.002$],
      [$0.836 plus.minus 0.002$],
  ),
  caption: [$R_A$ vs the sideband-subtracted $chevron.l z chevron.r$,
  inverse-variance weighted over the $approx 810$ surviving 4D bins per
  target, statistical errors only. The bottom row is the average over all
  bins. The top-bin $chevron.l z chevron.r = 0.53$ is the *true mean* of the
  wide top $z$ box $[0.37, 1.0]$, not its midpoint (@sec:binning-caveat).],
) <tab:R-vs-z>

Read as an inverse-variance weighted average, the ratio shows genuine
attenuation, with the qualitative shape expected of nuclear hadronization:

/ Suppression is real: all three overall ratios lie well below unity --- C the
  least ($0.925$), Cu and Sn near $0.83$.
/ It grows with $z$ for the heavy target: $R_"Sn"$ falls monotonically from
  $1.08$ at $chevron.l z chevron.r = 0.11$ to $0.67$ in the top $z$ bin --- the
  expected high-$z$ attenuation, and the cleanest single trend in the ratio.
/ Low-$z$ enhancement: $R_"Sn" > 1$ at the lowest $z$, consistent with energy
  removed from the leading hadron reappearing as additional soft production.

#wide-figure(
  "../figures/results_RA_vs_z.pdf",
  [$R_A$ vs the sideband-subtracted $chevron.l z chevron.r$ (diagnostic). Sn
  separates downward at large $z$; the rightmost point of each curve is the top
  $z$ box at its true $chevron.l z chevron.r approx 0.53$, not the geometric
  midpoint. Statistical errors only, correlated through LD#sub[2].],
  <fig:R-vs-z>,
  width: 78%,
)

#warning-box(title: "Weighted mean and median disagree — read before quoting R_A")[
  The attenuation above is in the *inverse-variance weighted* mean. The
  *median* $R$ over the same bins is $0.99$ (C), $0.95$ (Cu), $0.97$ (Sn)
  --- consistent with unity. The two disagree because the suppression lives in
  the high-statistics bins that dominate the weighted mean, while the many
  thin, noisy bins scatter symmetrically about one and dominate the median.
  With this note's own caution about weighted averages in mind (the @BSA
  section makes the same point), the $R_A$ signal is *suggestive, not
  established*: the ordering above C is not resolved ($R_"Cu" approx R_"Sn"$),
  and the whole picture rests on error weights from a truncated,
  background-subtracted, fallback-photon sample. More statistics --- and a $z$
  binning that resolves the top box rather than folding $z in [0.37, 1.0]$ into
  one bin --- are the way to settle it.
]

== Beam-spin asymmetry <sec:results-bsa>

Not yet extracted. `pi0.bsa` refuses to run without a measured beam
polarization (@sec:bsa-polarization) --- there is deliberately no placeholder,
the old code's self-declared $P = 0.85$ having been removed --- and the @RGD
value is not yet recorded in `cuts.json`. The Stage B `bsa` tree is written and
waiting; the asymmetry will be filled in here once $P plus.minus sigma_P$ is
supplied.
