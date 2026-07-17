#import "../template/lib.typ": *

= Event selection <sec:event-selection>

Selection proceeds in three stages --- electron identification, @DIS
kinematics, photon identification --- all of them in the C++ skim. Nothing
downstream can loosen them.

#note-box(title: "Provenance of the distributions in this section")[
  The figures throughout @sec:event-selection are *illustrative*: they apply
  the selection documented here to a single LD#sub[2] run (18431), truncated to
  its first $150 000$ events. That is enough to show the *shape* of each cut and
  what it removes, and nothing more --- it is one run of one target, not the
  production sample, so no rate, ratio or efficiency may be read from it. The
  photons in particular are classified by the @RGA fallback model
  (@sec:photon-id), so none of these figures is a photon-efficiency measurement.
  The "accepted" / "rejected" populations are shown at the *candidate* level,
  before and after the cut under discussion.
]

== Trigger electron <sec:trigger-electron>

The scattered electron is the first `REC::Particle` row satisfying
$"pid" = 11$, $"status" < 0$, $2000 <= |"status"| < 4000$ (i.e.\ in the
@FD), with a non-zero momentum vector.

== Electron identification <sec:electron-id>

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    table.header([*Cut*], [*Value*], [*Comment*]),
    [$chi^2_"PID"$], [$-5 < chi^2_"PID" < 5$], [Standard @EB PID quality.],
    [Momentum], [$p > 2.0$ GeV], [See the warning below --- *not* $0.8$ GeV.],
    [Vertex $z$], [target-dependent], [@sec:vertex.],
    [@SF], [$mu(p) plus.minus 3.5 sigma(p)$], [Per sector, per polarity;
      @eq:sf.],
    [@PCAL fiducial], [$l_v > 9.0$ cm, $l_w > 9.0$ cm], [Loose level.
      Photons use $14.0$ cm (@sec:photon-id).],
    [@DC edge], [$R_1 > 1.680$, $R_2 > 2.0$, $R_3 > 8.750$ cm],
      [`REC::Traj` edge distance, layers 6 / 18 / 36.],
  ),
  caption: [Electron selection, applied in this order with short-circuit
  evaluation. All comparisons are strict. There is *no* @ECAL (ECIN/ECOUT)
  fiducial cut, *no* polar-angle cut, and *no* @HTCC photoelectron cut ---
  the @PCAL $l_v \/ l_w$ pair is the only calorimeter fiducial, and the
  @SF cut is the only $e\/pi$ separation.],
) <tab:electron-cuts>

#warning-box(title: "The momentum cutflow label is wrong")[
  The cutflow row is printed as `"Momentum > 0.8 GeV"`, but the applied cut
  is `p > MIN_ELECTRON_MOMENTUM` with
  `MIN_ELECTRON_MOMENTUM = 2.0` GeV. The label is a stale string, repeated
  in the electron-analysis algorithm. *The effective cut is 2.0 GeV.* Do
  not quote the label --- and note that any cutflow table copied from the
  log inherits the error.
]

#wide-figure(
  "../figures/sel_electron_cutflow.pdf",
  [Where trigger electrons are lost, in the order the cuts run
  (@tab:electron-cuts), with short-circuit evaluation --- so each bar counts
  only the candidates a cut is the *first* to reject. The momentum cut
  dominates, and it is drawn at its true $2.0$ GeV value, not the stale
  $0.8$ GeV log label. The @PCAL fiducial removes almost nothing at the loose
  $9$ cm level; the sampling-fraction and @DC edge cuts trim the tails.],
  <fig:sel-cutflow>,
  width: 82%,
)

=== Sampling fraction <sec:sampling-fraction>

The sampling fraction is built from all three @ECAL layers,

$ "SF" = (E_"PCAL" + E_"ECIN" + E_"ECOUT") / p , $ <eq:sf-def>

and compared against a per-sector, per-polarity parameterisation of its
mean and width, each a third-order polynomial in $p$:

$ mu(p) = a_mu + b_mu p + c_mu p^2 + d_mu p^3 , quad
  sigma(p) = a_sigma + b_sigma p + c_sigma p^2 + d_sigma p^3 . $ <eq:sf>

The cut is $mu(p) - 3.5 sigma(p) < "SF" < mu(p) + 3.5 sigma(p)$, requiring
a valid @PCAL sector in $[1,6]$. The coefficients for both polarities are
tabulated in @sec:appendix-sf. They are target-averaged: one set per
sector and polarity, applied to LD#sub[2], C, Cu and Sn alike.

#wide-figure(
  "../figures/sel_sampling_fraction.pdf",
  [Sampling fraction @eq:sf-def against electron momentum. The accepted band
  is the $mu(p) plus.minus 3.5 sigma(p)$ envelope of @eq:sf; rejected
  candidates (below the band, and everything with $p < 2$ GeV) are the
  minimum-ionising and mis-identified tracks the cut removes. The band widens
  slightly toward low $p$, as the polynomial $sigma(p)$ intends.],
  <fig:sel-sf>,
  width: 68%,
)

== Vertex and target separation <sec:vertex>

The vertex cut both selects the target and, for the solid foils, *defines*
which nucleus the event belongs to. Cu and Sn are two foils in one
assembly; they are distinguished only by $v_z$.

For the solid targets the raw vertex is first corrected by a per-cell
(sector, $p$-bin, $phi$-bin) cubic-in-$theta$ parameterisation, blended
between two Gaussian hypotheses with a sigmoid weight. The corrected
vertex is then required to lie within $3 sigma$ of the appropriate peak:

$ |v_z^"corr" - mu| < 3.0 sigma . $ <eq:vz>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, left, center, center),
    table.header([*Target*], [*Rule*], [*Peak $mu$ (cm)*], [*Window (cm)*]),
    [LD#sub[2]], [raw $v_z$, correction disabled], [---], [$(-15.0, 5.0)$],
    [Cu], [corrected, upstream peak only], [$-7.861$], [$(-9.106, -6.616)$],
    [Sn], [corrected, downstream peak only], [$-2.916$], [$(-4.026, -1.806)$],
    [CxC], [corrected, either peak], [$-7.887$ / $-2.906$],
      [union of both $plus.minus 3 sigma$],
  ),
  caption: [Vertex windows for *outbending* running, the only polarity used
  here. Widths are $sigma = 0.415$ (Cu), $0.370$ (Sn), $0.395 \/ 0.373$
  (CxC). The CxC target is two carbon
  foils occupying the same two $z$ slots as the Cu and Sn foils, which is
  why it accepts either peak. Inbending values are given in
  @sec:appendix-vz.],
) <tab:vz>

#note-box(title: "Two naming traps in the vertex code")[
  + The internal names `lo`/`hi` refer to the *peak parameter set*, not to
    the $z$ ordering: `mu_hi` $= -7.861$ is the *more negative*, i.e.\
    further upstream, position. The doc-comment calling the `hi` peak
    "downstream" is backwards. The code's assignment (Cu $arrow.l.r -7.86$,
    Sn $arrow.l.r -2.92$) is self-consistent and agrees with the legacy
    lookup table; trust the numbers, not the prose.
  + The service *overrides* the peak parameters read from the parameter
    file with hard-coded values --- e.g.\ the file gives
    $sigma_"hi" = 0.385$ for CuSn outbending where the service uses
    $0.415$, an $8%$ difference. The correction polynomials come from the
    file; the cut window does not. @tab:vz quotes the values that are
    actually applied.
]

A separate, legacy raw-$v_z$ table exists and is used by the standalone
`electron_analysis` binary. It is *not* the $pi^0$ path, and for the solid
targets it gives materially different windows. Do not use
`electron_analysis` as a proxy for this analysis's electron selection.

#wide-figure(
  "../figures/sel_vz.pdf",
  [Corrected electron vertex for the LD#sub[2] run, log scale. The accepted
  band (blue) drops sharply at the raw LD#sub[2] window edges $(-15, 5)$ cm
  (@tab:vz, dotted); the Cu ($plus.minus 3 sigma$ about $-7.86$) and Sn (about
  $-2.92$) windows are overlaid to show they are *disjoint* --- for the solid
  targets $v_z$ alone assigns the nucleus, so the two foils never share an
  event. The tails beyond the axis ($4.8%$, annotated) are misreconstructed
  tracks. For LD#sub[2] the correction is disabled, so this is the raw $v_z$.],
  <fig:sel-vz>,
  width: 70%,
)

== DIS kinematics <sec:dis>

@DIS variables are computed from the scattered electron against a
*stationary proton*, for every target:

$ Q^2 = 4 E_"beam" E'_e sin^2 (theta_e \/ 2) , $ <eq:Q2>
$ nu = E_"beam" - E'_e , quad y = nu \/ E_"beam" , $ <eq:nu>
$ x_B = Q^2 / (2 M_p nu) , quad W = sqrt((q + P_"target")^2) . $ <eq:xB>

with $M_p = 0.938272$ GeV and $E_"beam" = 10.53$ GeV. Two implementation
details are worth stating because they are inherited rather than chosen:
$Q^2$ uses the massless/angle form rather than $-(k - k')^2$ (equivalent to
better than the resolution at these energies), and $x_B$ uses the *proton*
mass for the nuclear targets too --- conventional for multiplicity ratios,
where the per-nucleon variable is what cancels, but worth stating
explicitly.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    table.header([*Variable*], [*Applied cut*], [*Motivation*]),
    [$Q^2$], [$> 1.0$ GeV#super[2]], [Hard-scattering regime.],
    [$W$], [$> 2.0$ GeV], [Above the resonance region.],
    [$y$], [$< 0.85$], [Suppresses large radiative corrections and poor
      resolution.],
  ),
  caption: [Effective @DIS selection. All cuts are strict. The
  configuration also sets $Q^2 < 100$ GeV#super[2] and $W < 100$ GeV, but at
  $E_"beam" = 10.53$ GeV the kinematic maxima are $Q^2 approx 9.6$
  GeV#super[2] and $W approx 4.6$ GeV, so neither upper bound is ever
  reached. There is no $x_B$ cut.],
) <tab:dis-cuts>

#note-box(title: "`dis_region` is dead configuration")[
  Both $pi^0$ binaries construct the @DIS algorithm with
  `dis_region = standard_dis_region()`, which would imply
  $Q^2 in [1, 20]$ and $W in [2, 25]$. *That field is never read.* The cut
  code consults only the flat `q2_min/q2_max`, `w_min/w_max`, `y_min/y_max`
  members, which the configuration file sets to $1\/100$, $2\/100$,
  $0\/0.85$.

  The distinction has no physics consequence here --- none of $20$, $25$,
  $100$ is reachable at this beam energy, so all four candidate bounds
  select the same events. It matters only for correctness of the record: a
  note quoting "$Q^2 < 20$, $W < 25$" would be describing a cut that does
  not run. The same is true of a handful of other declared-but-unread
  configuration keys, collected in @sec:appendix-dead-config.
]

#subfig2(
  (
    ("../figures/sel_dis_phase_space.pdf", [$Q^2$ versus $x_B$, with the three cuts drawn as the boundaries they are.]),
    ("../figures/sel_electron_kinematics.pdf", [Scattered-electron $p$ versus $theta_e$, with the $p > 2$ GeV cut.]),
  ),
  [The @DIS acceptance of the accepted sample. *Left:* the $W = 2$ GeV and
  $y = 0.85$ cuts are curves in the $(x_B, Q^2)$ plane, not axis limits, and
  the events fill exactly the wedge between them and $Q^2 > 1$; the strong
  $Q^2$--$x_B$ correlation is why a rectangular product grid is mostly empty
  and an adaptive binning was chosen (@sec:binning). *Right:* the accepted
  electrons sit above $p > 2$ GeV along the usual @DIS $p$--$theta$ band.],
  <fig:sel-dis>,
)

#wide-figure(
  "../figures/sel_dis_1d.pdf",
  [The @DIS variables of the accepted sample against their cuts
  (@tab:dis-cuts). $Q^2$ and $W$ press against their lower bounds; $y$ reaches
  only $~0.82$, comfortably inside the $0.85$ ceiling, so on this sample the
  $y$ cut removes nothing --- consistent with its role as a soft
  radiative-correction guard rather than an active edge. $nu$ has no cut of
  its own.],
  <fig:sel-dis1d>,
  width: 84%,
)

== Photon identification <sec:photon-id>

Photons are selected once, at skim time, by a Catboost @GBT classifier
rather than by cut-based identification.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    table.header([*Cut*], [*Value*], [*Comment*]),
    [PID], [$= 22$], [@EB neutral, calorimeter-matched.],
    [AI classifier], [$sigma("model") > 0.78$],
      [$sigma(x) = 1 \/ (1 + e^(-x))$; see below.],
    [Energy], [$E_gamma >= 0.2$ GeV], [Inside the AI pre-filter.],
    [@PCAL energy], [$E_"PCAL" > 0$], [Shower must start in the preshower.],
    [Polar angle], [$5degree <= theta_gamma <= 35degree$], [Forward
      detector.],
    [@PCAL fiducial], [$l_v > 14.0$ cm, $l_w > 14.0$ cm], [Harder than the
      electron's $9.0$ cm.],
    [$beta$], [$0.9 < beta < 1.1$], [Consistency with light speed.],
  ),
  caption: [Photon selection. The energy, @PCAL energy and polar-angle
  requirements sit inside the classifier's purity pre-filter, so they apply
  only on the AI path --- which is the production path. Photon energy is
  taken as $E_gamma = |arrow(p)_gamma|$ (massless). There is *no* timing
  cut anywhere in the $pi^0$ chain; $beta$ is the only timing-derived
  handle.],
) <tab:photon-cuts>

The classifier takes a *45*-component feature vector: seven variables for
the candidate itself ($E_gamma$, $E_"PCAL"$, $theta$, the shower second
moments $m_(2u)$ and $m_(2v)$, and the angular distance and energy
difference to the scattered electron), then five variables --- angular
distance, energy difference, @PCAL energy and the two second moments ---
for each of the three nearest photons, two nearest charged hadrons and two
nearest neutral hadrons ($5 times 7 = 35$), and finally three counts of
neighbouring photons within $Delta R < 0.1$, $0.2$ and $0.35$ rad. Angular
distances are measured between calorimeter hit positions, taking the
@PCAL position where available and falling back to @ECIN then @ECOUT. Its purpose is to
reject the neutral background, principally neutrons and split-off clusters,
that a simple energy threshold cannot.

#warning-box(title: "RG-D photons are classified by an RG-A model")[
  The classifier is selected from a lookup table keyed on run range and
  cooking pass. The table covers @RGA and @RGC only:

  #align(center)[
    #text(size: 9.5pt)[
      runs 5032--5666, 6156--6783, 11093--11571 (@RGA), 16042--16772 (@RGC)
    ]
  ]

  @RGD runs span *18305--19131* (@tab:runs) and therefore match *no*
  entry. The lookup falls through to a silent default:

  #align(center)[
    #text(size: 9.5pt, font: ("JetBrains Mono", "Menlo"))[
      `return ApplyCatboostModel_RGA_inbending_pass1(data);`
    ]
  ]

  So *every RG-D photon in this analysis --- outbending data, nuclear
  targets --- is identified by a model trained on RG-A inbending pass-1
  data.* No warning is logged. The pass index is likewise hard-coded to 1
  and never set from the configuration.

  This is the single largest known systematic on the photon efficiency, and
  it is entirely unquantified. It does not obviously bias the *ratio*
  $R_A$, since the same model is applied to every target and the efficiency
  should cancel to first order --- but that cancellation is an assumption,
  not a measurement, and it fails to the extent that photon efficiency
  depends on the local occupancy, which does differ between LD#sub[2] and
  Sn. Quantifying it requires either a model trained on @RGD or a
  data-driven efficiency study.
]

#wide-figure(
  "../figures/sel_gbt_score.pdf",
  [The classifier score $sigma("model")$ for photon candidates, log scale,
  coloured by the $0.78$ threshold. The distribution is bimodal --- a background
  population piling up near zero, a signal peak near one --- and the threshold
  sits in the valley between them. Clusters rejected by the purity pre-filter
  (energy, @PCAL energy, polar angle) are never scored and are counted in the
  title, not drawn at a false zero. Of the $15\,841$ clusters above threshold, a
  further $~1700$ are removed by the later $beta$ and $14$ cm @PCAL fiducial
  cuts (@tab:photon-cuts), so $14\,127$ photons enter the slim. *This is the
  @RGA fallback model applied to @RGD data* (see the box above); the clean
  separation shown here is the model's behaviour on data it was not trained on,
  not evidence that it is correct for @RGD.],
  <fig:sel-gbt>,
  width: 68%,
)

#wide-figure(
  "../figures/sel_photon.pdf",
  [Selected photons. *Left:* the energy spectrum, falling steeply from the
  $0.2$ GeV pre-filter floor. *Centre:* the polar angle, spanning the
  $5degree$--$35degree$ forward-detector window. *Right:* the angle to the
  scattered electron; the $theta_(e gamma) > 8degree$ pair-level cut
  (@sec:pair-cuts) removes the small population overlapping the electron, and
  is drawn here even though it is applied at pairing rather than in the skim.],
  <fig:sel-photon>,
  width: 92%,
)

== What the skim costs <sec:skim-cost>

The skim writes only two banks (@sec:appendix-banks): per-event @DIS
kinematics plus the trigger electron, and the four-momenta of the selected
photons. This is a $~100 times$ reduction and it is what makes the
downstream iteration fast.

The price is that *the photon selection is frozen*. The skim carries no
PID, status, $beta$ or vertex column for photons, so the downstream
analysis cannot re-cut or vary them. Re-tuning the AI threshold, the
$beta$ window or the @PCAL fiducial --- each of which is a systematic that
ought to be varied --- requires re-running the skim over the full dataset.
The slim reader also rebuilds photons with $beta$ and vertex zeroed, so the
slim path cannot even reproduce the $beta$ cut it inherited.
