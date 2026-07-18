#import "../template/lib.typ": *

= Dataset and data flow <sec:dataset>

== Beam and targets <sec:targets>

The analysis uses @RGD production data: a longitudinally polarised electron
beam of $E_"beam" = 10.53$ GeV incident on the RG-D multi-target assembly.
Four target configurations enter this note:

#figure(
  table(
    columns: (auto, auto, auto, 1fr),
    align: (left, center, center, left),
    table.header(
      [*Target*], [*$A$*], [*Role*], [*Comment*],
    ),
    [LD#sub[2]], [2], [reference], [Liquid deuterium cell; the denominator of
      every ratio in this note.],
    [CxC], [12], [nuclear], [Double carbon foil. Written
      `CxC` throughout the code and in the result files.],
    [Cu], [63.5], [nuclear], [Solid foil.],
    [Sn], [118.7], [nuclear], [Solid foil.],
  ),
  caption: [Target set. *Each run used a single target*: LD#sub[2], CxC and
  CuSn were taken in separate, non-overlapping runs (@sec:runs), so the ratio
  of @eq:multiplicity compares data recorded at different times. Cu and Sn
  are the exception --- they are two foils in one CuSn assembly, exposed in
  the *same* runs and separated offline by the electron vertex $v_z$
  (@sec:vertex). Stage A writes one slim file per target.],
) <tab:targets>

The Cu/Sn exception is worth keeping in mind when reading @sec:results: of
the three nuclear-to-deuterium ratios, $R_"Cu"$ and $R_"Sn"$ share their
running conditions with each other exactly, while $R_"CxC"$ does not share
them with either. The Cu-to-Sn *comparison* is therefore free of the
run-to-run systematic that afflicts every ratio against LD#sub[2].

*Torus polarity.* Every result in this note comes from *outbending* data
only. The run-list filter (@sec:runs) selects the outbending runs; the
inbending set is present in `config/runs.json` but no production reads it.
The inbending dataset therefore remains an untapped statistical reserve, and
--- more importantly --- a polarity comparison is one of the few
detector-systematic handles available; it has not been exercised.

== Run ranges <sec:runs>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, right),
    table.header([*Target*], [*Polarity*], [*Run range*], [*$N_"runs"$*]),
    [LD#sub[2]],            [Outbending], [18419 -- 19058], [135],
    [CxC], [Outbending], [18440 -- 18850], [112],
    [CuSn],                 [Outbending], [18561 -- 19131], [268],
    [LD#sub[2]],            [Inbending],  [18305 -- 18336], [24],
    [CxC], [Inbending],  [18339 -- 18346], [8],
    [CuSn],                 [Inbending],  [18347 -- 18394], [27],
  ),
  caption: [@RGD run lists as recorded in `config/runs.json`. Cu and Sn share
  the `CuSn` run list --- they are two foils in one assembly, separated
  offline by vertex (@sec:vertex). *These lists are applied*: `pi0.farm`
  refuses any input file whose run is not in the target's list for the
  requested polarity, because the SIDIS-train directory tree holds inbending
  and outbending runs together and the slim schema records no polarity
  (@sec:provenance-gaps). The file also records a `special_DC_runs_CuSn` set
  (18358--18368) and six small LD#sub[2] trigger-configuration sets, none of
  which enters the outbending production.],
) <tab:runs>

The three outbending lists are *pairwise disjoint* --- no run appears under
two targets --- confirming that each run carried a single target. Together
they are $135 + 112 + 268 = 515$ distinct runs spanning 18419--19131.

Ordering those 515 runs and labelling each by its target resolves them into
*13 contiguous single-target blocks*:

#figure(
  block(inset: (y: 3pt))[
    #set text(size: 9pt, font: ("JetBrains Mono", "Menlo"))
    #set par(justify: false)
    #table(
      columns: (auto, auto, auto, auto, auto, auto),
      align: (left + horizon,) * 6,
      stroke: none,
      inset: (x: 5pt, y: 2.5pt),
      [LD#sub[2]], [18419--18439], [20],  [CxC],  [18440--18524], [67],
      [LD#sub[2]], [18528--18559], [29],  [CuSn], [18561--18642], [61],
      [LD#sub[2]], [18644--18656], [13],  [CuSn], [18660--18755], [72],
      [CxC],       [18756--18762], [5],   [LD#sub[2]], [18764--18790], [27],
      [CxC],       [18796--18850], [40],  [LD#sub[2]], [18851--18873], [18],
      [CuSn],      [18875--18966], [76],  [LD#sub[2]], [19021--19058], [28],
      [CuSn],      [19062--19131], [59],  [], [], [],
    )
  ],
  caption: [The outbending running, in run order: target, run span, number of
  runs. The beam returns to LD#sub[2] *six* times, so the reference is spread
  across the whole period rather than concentrated at one end of it. The
  nuclear blocks are long (40--76 runs), so conditions can drift appreciably
  *within* a block as well as between blocks.],
) <tab:run-blocks>

#result-box(title: "A free systematic probe: six LD₂ blocks")[
  The six separate LD#sub[2] blocks are the most useful thing in
  @tab:run-blocks. They are the *same target* measured at six points spread
  across the run period, so any difference between them is pure
  time-dependent systematic --- there is no physics that should make
  deuterium behave differently in April than in June.

  Forming $R$ of one LD#sub[2] block against another therefore measures,
  directly and with no simulation, the run-to-run systematic that
  @sec:ratio-cancellation says is unevaluated. It needs no new data and no
  new method: it is the existing machinery with the LD#sub[2] sample split
  by run. This is the single cheapest systematic available to this analysis
  and it should be done before anything else in @tab:systematics.
]

== Processing chain <sec:dataflow>

The analysis is one repository --- `github.com/mathieuouillon/rgd-pi0`, C++20
and Python, LGPL-3.0 --- run in four stages. Each stage stamps its
configuration and its inputs into its output, and refuses to build physics on
an input it cannot vouch for (@sec:provenance-gaps).

#figure(
  block(inset: (y: 4pt))[
    #set text(size: 9.5pt)
    #table(
      columns: (auto, 1fr, auto),
      align: (left, left, left),
      table.header([*Stage*], [*What it does*], [*Output*]),
      [A. Skim], [`stageA_skim` (hipo4 + RDataFrame) reads @HIPO DST files ---
        the SIDIS train under
        `/cache/clas12/rg-d/production/pass1/recon/<T>/dst/train/SIDIS/` ---
        applies electron, @DIS and @GBT photon selection, and writes a slim
        `events` TTree. Driven by `pi0.batch` (@SWIF2) on the farm or
        `pi0.local_batch` interactively.],
        [slim `events` ROOT],
      [B. Grid], [`make_grid` scans the slims and writes the frozen
        equal-statistics binning: Grid A $(Q^2, x_B)$ and Grid B $(z, p_T^2)$
        (@sec:binning). The edges and their hash are committed to
        `config/binning/`.],
        [`grid_A_q2_xb.json`, `grid_B_z_pt2.json`],
      [C. Bin], [`stageB_bin` reads the slims, reconstructs $pi^0 -> gamma
        gamma$, computes @SIDIS kinematics, builds the frozen mixed-event
        donor pool (@sec:background) and fills four flat trees on the grid.],
        [binned ROOT (`spectra`, `ptb3d`, `n_dis`, `bsa`)],
      [D. Extract], [`pi0.ratio`, `pi0.broadening`, `pi0.bsa`, `pi0.qa`
        (uproot) read the binned trees and extract the physics: yields,
        $R_A$, $Delta chevron.l p_T^2 chevron.r$, $A_"LU"$.],
        [CSV + PDF],
    )
  ],
  caption: [The four processing stages. A--C are C++20 (static library
  `pi0_core` plus the three executables); D is Python (`pi0` package). The
  division matters for reading this note: *every kinematic cut lives in
  `config/cuts.json`* and is applied in stages A--C; the Python layer applies
  no kinematic selection beyond a $m_(gamma gamma)$ window and fit-quality
  gates.],
) <tab:dataflow>

The skim exists because it is a large reduction: only the trigger electron's
four-momentum, the @DIS kinematics, the helicity and the selected photons
survive. That makes the downstream analysis fast to iterate, at a price
documented in @sec:skim-cost --- the photon selection is baked in and cannot
be varied downstream.

*Program freezing.* `pi0.batch` snapshots the executable, its shared
libraries and the configuration into a per-workflow directory and runs @SWIF2
jobs against that copy, so a queued workflow is decoupled from ongoing
development of the source tree.

== Statistics <sec:statistics>

The numbers in @sec:results are a *diagnostic first run* of the whole chain,
not the full luminosity: Stage A read only the first $2 times 10^6$ events of
each of the 95 SIDIS-train files (a few hours on one interactive node), so the
counts below are a prefix of the available data and carry the truncation
caveat throughout.

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, right, right),
    table.header([*Target*], [*$N^"DIS"$*], [*$pi^0$ (same-event, binned)*]),
    [LD#sub[2]], [$1 ,987 ,072$], [$882 ,869$],
    [CxC], [$1 ,157 ,245$], [$521 ,954$],
    [Cu], [$1 ,045 ,316$], [$461 ,529$],
    [Sn], [$1 ,299 ,764$], [$577 ,531$],
  ),
  caption: [Statistics of the diagnostic run. $N^"DIS"$ is the inclusive
  normalisation of @eq:multiplicity (events read per target); the $pi^0$
  column is the same-event $gamma gamma$ pairs binned by Stage B, *before*
  background subtraction --- an upper bound on the signal, not a $pi^0$
  count. Under the $2 times 10^6$-event/file cap these total
  $approx 5.5 times 10^6$ @DIS events and $2.5 times 10^6$ $pi^0$ candidates.
  A full-luminosity run reuses the same machinery with the cap removed.],
) <tab:statistics>

#wide-figure(
  "../figures/kinematics_pi0.pdf",
  [Kinematic coverage of the reconstructed $pi^0$ sample. The $Q^2$ versus
  $x_B$ correlation band is why the factorized $(Q^2, x_B)$ grid of
  @sec:binning leaves its corner cells empty by construction: a product grid
  cannot follow a correlated distribution the way the equal-statistics *axes*
  can, and that trade is made deliberately in exchange for global, quotable,
  version-controlled edges.],
  <fig:kinematics>,
  width: 92%,
)

== Provenance <sec:provenance-gaps>

Reproducibility is a design goal of this codebase, not an afterthought, and
the pieces a published note needs are in the repository rather than on a
scratch disk.

#result-box(title: "What every result can be traced to")[
  + *The run list is applied, and recorded.* `config/runs.json` holds the
    @RGD run lists (@tab:runs), and `pi0.farm` filters every input against
    them: a file whose run is not in the requested target/polarity list is
    rejected with a reason. Inbending and outbending runs share one directory
    tree and the slim schema records no polarity, so this filter is what keeps
    the two torus polarities from being silently mixed.

  + *The configuration is committed.* Every threshold lives in
    `config/cuts.json`; a cut value appearing in code is a bug. Its SHA-256
    (`801ba433…`) is stamped into every slim, every binned file and every
    result CSV, and each stage refuses to run if the config it loaded is not
    the one its inputs were made with.

  + *The binning is committed and hashed.* `make_grid` writes
    `config/binning/grid_{A,B}.json` with a provenance block --- the dataset,
    the targets, the event count and a hash of the edges
    (`provenance_hash 2acf618294a6c3b0` for the grids used here) --- and that
    hash is stamped into every output binned on them. The edges are global,
    quotable and under version control.

  + *Unpublishable inputs are refused, loudly.* Stage A stamps
    `gbt.fallback_used` and any `--max-events` truncation into its provenance;
    Stage B propagates one `provenance_stageA_NNN` block per input; the Python
    stage reads them and refuses to produce a physics number from a fallback,
    truncated or placeholder-grid input unless explicitly told the result is
    diagnostic (`--allow-unpublishable`). Every diagnostic number in
    @sec:results carries that stamp.
]
