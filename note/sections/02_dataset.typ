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
  (@sec:vertex). The C++ stage writes one file per target.],
) <tab:targets>

The Cu/Sn exception is worth keeping in mind when reading @sec:results: of
the three nuclear-to-deuterium ratios, $R_"Cu"$ and $R_"Sn"$ share their
running conditions with each other exactly, while $R_"CxC"$ does not share
them with either. The Cu-to-Sn *comparison* is therefore free of the
run-to-run systematic that afflicts every ratio against LD#sub[2].

*Torus polarity.* Every result in this note comes from *outbending* data
only. All four production files carry the `_OB` suffix, and no inbending
file is read by any script. A `Polarity.INBENDING` enumerator exists in the
configuration module but is unused in the $pi^0$ chain. The inbending
dataset therefore remains an untapped statistical reserve, and --- more
importantly --- a polarity comparison is one of the few
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
  caption: [@RGD run lists as recorded in `clas12/Runs.hpp`, with run counts
  obtained by parsing the header directly. Cu and Sn share the `CuSn` run
  list --- they are two foils in one assembly, separated offline by vertex
  (@sec:vertex). *These lists are not applied by any binary* (see
  @sec:provenance-gaps); they document the run ranges but do not constitute
  the production selection. The header additionally records a
  `special_DC_runs_CuSn` set (18358--18368, 10 runs) and six small LD#sub[2]
  trigger-configuration sets spanning 18419--18535, none of which is used to
  include or exclude anything.],
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

The analysis runs in three stages, in two separate repositories.

#figure(
  block(inset: (y: 4pt))[
    #set text(size: 9.5pt)
    #table(
      columns: (auto, 1fr, auto),
      align: (left, left, left),
      table.header([*Stage*], [*What it does*], [*Output*]),
      [1. Skim], [`pi0_multiplicity_skim` reads @HIPO DST files from tape
        (the SIDIS train under
        `/mss/clas12/rg-d/production/pass1/recon/<target>/dst/train/SIDIS/`),
        applies electron, @DIS and photon selection, and writes a slim
        two-bank @HIPO file. Submitted to the JLab farm as a @SWIF2
        workflow, one job per few input files, staged from `/mss` and
        returned to `/volatile`.],
        [`PI0::event`, `PI0::photons`],
      [2. Analysis], [`pi0_multiplicity` reads the slim skim, reconstructs
        $pi^0 -> gamma gamma$, computes @SIDIS kinematics, and fills the
        binned histograms and the `pi0_pairs` ntuple. Run per target via
        `run_pi0_batch.py`.],
        [`pi0_multiplicity_<T>_OB.root`],
      [3. Extraction], [Python scripts read the ROOT histograms/ntuple and
        extract the physics: yields, $R_A$, $Delta chevron.l p_T^2 chevron.r$,
        $A_"LU"$.],
        [CSV + PDF],
    )
  ],
  caption: [The three processing stages. Stages 1--2 are C++20
  (`clas-framework`); stage 3 is Python (`analyses/pi0/`). The division
  matters for reading this note: *every kinematic cut lives in stages 1--2*,
  and the Python layer applies no kinematic selection at all beyond a
  $m_(gamma gamma)$ window and fit-quality gates.],
) <tab:dataflow>

The skim exists because it is a large reduction: only the trigger
electron's four-momentum, the @DIS kinematics, the helicity and the
selected photons survive. That makes the downstream analysis fast to
iterate, at a price documented in @sec:skim-cost --- the photon selection
is baked in and cannot be varied downstream.

*Program freezing.* The submission script snapshots the executable, its
shared libraries and the configuration into a per-workflow directory and
runs jobs against that copy, so a queued workflow is decoupled from
ongoing development of the source tree. This is good practice and is worth
recording: it means the code state of a given production is recoverable
from the frozen directory, not from `git`.

== Code provenance <sec:code-provenance>

Both C++ working directories on the analysis machine are checkouts of the
same repository (`gitlab.com/MathieuOuillon/clas-framework`) at *different
commits*, and they differ in exactly the place that matters most for this
note.

#figure(
  table(
    columns: (auto, auto, auto, 1fr),
    align: (left, left, left, left),
    table.header([*Checkout*], [*HEAD*], [*Date*], [*$pi^0$ binning*]),
    [`clas-framework`], [`e8334b1`], [16 Jun 2026],
      [kd-tree (`AdaptiveBinner4D`)],
    [`clas-analysis-1`], [`29150c0`], [14 Jul 2026],
      [factorized $(Q^2, x_B) times (z, p_T^2)$ grids],
  ),
  caption: [The two checkouts. `e8334b1` is a direct ancestor of `29150c0`:
  the newer tree is three commits ahead, one of which
  (`eac29a1`, _"replace kd-tree pi0 binning with factorized
  $(Q^2,x_B) times (z,p_T^2)$ grids"_) *deletes* `AdaptiveBinner4D.hpp`
  outright.],
) <tab:checkouts>

#important-box(title: "Which code this note describes")[
  *This note documents the kd-tree binning*, because that is what produced
  every number and figure presented here. The evidence is unambiguous: all
  result files carry a `leaf` index with data-driven quantile box edges
  (@sec:binning), and they are dated 8--16 June 2026 --- contemporaneous
  with `e8334b1` (16 June), and a month before the refactor (14 July).

  The kd-tree machinery *no longer exists* in `clas-analysis-1`. Anyone
  reading that checkout to understand the current results will find a
  factorized grid that has produced none of them. The refactor is
  summarised in @sec:binning-future; it changes the binning for the *next*
  production and will require every number in @sec:results to be
  regenerated.
]

== Statistics <sec:statistics>

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, right, right),
    table.header([*Target*], [*$N^"DIS"$*], [*$pi^0$ pairs in window*]),
    [LD#sub[2]], [$23 ,098 ,120$], [$15.73 times 10^6$],
    [CxC], [$13 ,642 ,531$], [$9.72 times 10^6$],
    [Cu], [$12 ,277 ,770$], [$8.52 times 10^6$],
    [Sn], [$15 ,215 ,380$], [$10.59 times 10^6$],
  ),
  caption: [Integrated statistics. @DIS counts are the inclusive
  normalisation of @eq:multiplicity, summed over the 58 populated
  $(Q^2, x_B)$ cells. The $pi^0$ column is the number of $gamma gamma$
  pairs inside $0.110 < m_(gamma gamma) < 0.160$ GeV summed over all
  kd-tree leaves entering the @BSA, and is *not* background-subtracted ---
  it is an upper bound on the signal, not a $pi^0$ count. The resulting
  statistical precision on $R_A$ is $1-2%$ per leaf.],
) <tab:statistics>

#wide-figure(
  "../figures/kinematics_pi0.pdf",
  [Kinematic coverage of the reconstructed $pi^0$ sample. The $Q^2$ versus
  $x_B$ correlation band is the reason the factorized normalisation
  fallback of @sec:normalisation is unsafe, and the reason an adaptive
  binning was chosen over a product grid in the first place: a rectangular
  grid spanning this acceptance is mostly empty.],
  <fig:kinematics>,
  width: 92%,
)

== Provenance gaps <sec:provenance-gaps>

Three pieces of provenance that a published note requires are *not
recoverable* from either repository and must be supplied before
circulation.

#warning-box(title: "Unrecoverable provenance")[
  + *The run list is recorded but never applied.* `clas12/Runs.hpp` contains
    explicit @RGD run lists (@tab:runs), but they are *dead code*: the only
    consumer would be `FileUtils::filter_by_run_numbers`, which has no
    caller anywhere in the repository. No binary performs run-by-run
    selection, bad-run vetoing, or trigger-configuration filtering. The runs
    that actually entered a production are therefore whatever file list was
    passed on the command line --- which is not recorded. The lists in
    @tab:runs should be treated as *documentation of the RG-D run ranges*,
    not as a description of the selection that ran.
    Likewise `rcdb_backup.sql` (27 MB) and `RGD_charge_final-MH-Nov24.xlsx`
    sit at the root of the ifarm analysis folder but are read by *no* $pi^0$
    script --- the charge spreadsheet is used only by the unrelated
    nuclear-transparency analysis.

  + *The production configuration is not in the repository.* No shipped
    config produces the file names that the Python actually reads
    (`pi0_multiplicity_<T>_OB.root`; the configs emit `..._Outbending.root`),
    and the produced @BSA binning tree has 250 leaves where the shipped
    config specifies a stratified $5 times 5 times 3 times 3 = 225$ split.
    The configuration that made the production files is therefore *not* the
    one in the repository.

  + *The binning tables live only on ifarm.* The kd-tree leaf boxes
    (`pi0_adaptive_binning*.txt.boxes.csv`) are written by the C++ pass into
    `/work` and exist nowhere else. They are an input to every Python script
    and are not archived with the results. The binning reported in
    @sec:binning was reconstructed from the output CSVs
    (@sec:binning-reconstruction).

  Additionally, none of `results/`, `analyses/`, `common/` or `plots/` is
  under version control --- the ifarm analysis folder has two commits and
  the analysis directories are untracked. The only provenance for every
  number quoted in @sec:results is a file modification time.
]
