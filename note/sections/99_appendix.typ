#import "../template/lib.typ": *

#set heading(numbering: "A.1")
#counter(heading).update(0)

= Sampling-fraction parameters <sec:appendix-sf>

Coefficients of @eq:sf, per sector and polarity. The cut is
$mu(p) - 3.5 sigma(p) < "SF" < mu(p) + 3.5 sigma(p)$ with
$"SF" = (E_"PCAL" + E_"ECIN" + E_"ECOUT") \/ p$. These are
target-averaged: one set per sector and polarity, applied to LD#sub[2], C,
Cu and Sn alike.

#sf-param-table(
  (
    ("Outbending", (
      ("1", "2.300896E-01", "4.862074E-03", "-6.441645E-04", "1.247456E-05",
             "2.226749E-02", "-3.186742E-03", "3.087387E-04", "-1.247691E-05"),
      ("2", "2.286452E-01", "4.985634E-03", "-3.668113E-04", "-1.110160E-05",
             "2.150934E-02", "-2.549354E-03", "1.855833E-04", "-4.394923E-06"),
      ("3", "2.258053E-01", "9.073840E-03", "-1.240745E-03", "3.579321E-05",
             "2.399538E-02", "-4.208668E-03", "4.503779E-04", "-1.783519E-05"),
      ("4", "2.281016E-01", "6.073295E-03", "-5.916479E-04", "2.569914E-06",
             "2.217527E-02", "-2.788589E-03", "2.266746E-04", "-7.717066E-06"),
      ("5", "2.302059E-01", "3.020237E-03", "-4.307581E-05", "-2.216133E-05",
             "2.277357E-02", "-3.941980E-03", "4.455085E-04", "-1.912118E-05"),
      ("6", "2.262583E-01", "6.849767E-03", "-8.027577E-04", "1.288497E-05",
             "2.253050E-02", "-3.334553E-03", "3.374127E-04", "-1.413483E-05"),
    )),
    ("Inbending", (
      ("1", "2.311088E-01", "4.309933E-03", "-5.451627E-04", "4.965823E-06",
             "2.235239E-02", "-3.117249E-03", "2.887205E-04", "-1.130317E-05"),
      ("2", "2.304337E-01", "3.585716E-03", "-1.672223E-04", "-7.317175E-06",
             "2.170118E-02", "-2.428037E-03", "1.302147E-04", "-2.503248E-09"),
      ("3", "2.269123E-01", "8.797582E-03", "-1.171416E-03", "3.051342E-05",
             "2.423467E-02", "-4.301928E-03", "4.742880E-04", "-1.983022E-05"),
      ("4", "2.288870E-01", "6.786520E-03", "-8.022602E-04", "1.557251E-05",
             "2.313626E-02", "-3.260271E-03", "3.018634E-04", "-1.156770E-05"),
      ("5", "2.288645E-01", "4.637362E-03", "-3.662892E-04", "-4.233509E-06",
             "2.310895E-02", "-3.892876E-03", "4.115173E-04", "-1.656769E-05"),
      ("6", "2.272093E-01", "6.581981E-03", "-8.197177E-04", "1.632178E-05",
             "2.321748E-02", "-3.665685E-03", "3.804588E-04", "-1.587616E-05"),
    )),
  ),
  [Sampling-fraction parameters for both torus polarities. Only the
  *Outbending* block is used by the results in this note
  (@sec:targets).],
  <tab:sf-params>,
)

= Vertex windows, inbending <sec:appendix-vz>

For completeness; unused here, since only outbending data is analysed.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    table.header([*Target*], [*$mu_"lo"$*], [*$sigma_"lo"$*],
                 [*$mu_"hi"$*], [*$sigma_"hi"$*]),
    [CuSn / Cu / Sn], [$-2.812$], [$0.315$], [$-7.765$], [$0.321$],
    [CxC], [$-2.713$], [$0.307$], [$-7.670$], [$0.313$],
  ),
  caption: [Inbending vertex peaks, in cm. The cut is
  $|v_z^"corr" - mu| < 3 sigma$. As in @tab:vz, `lo`/`hi` name the
  parameter set, not the $z$ ordering: the `hi` peak is the *upstream*
  (more negative) foil. Sn is the $-2.8$ peak, Cu the $-7.8$ peak.],
) <tab:vz-inbending>

= Skim bank definitions <sec:appendix-banks>

Both banks are generated from `data/bankdefs/hipo4/pi0skim.json`, group
`30000`.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, center, left),
    table.header([*Field*], [*Type*], [*Meaning*]),
    table.cell(colspan: 3)[*`PI0::event`* --- item 1, one row per event],
    [`run`], [`int32`], [Run number, from `RUN::config`],
    [`event`], [`int64`], [Event number],
    [`helicity`], [`int32`], [Beam helicity, HWP-corrected, $plus.minus 1$],
    [`Q2`], [`float`], [@DIS $Q^2$ (GeV#super[2])],
    [`xB`], [`float`], [Bjorken $x$],
    [`nu`], [`float`], [Energy transfer $nu$ (GeV)],
    [`W`], [`float`], [Invariant mass $W$ (GeV)],
    [`y`], [`float`], [Inelasticity],
    [`ex`, `ey`, `ez`], [`float`], [Trigger electron $p_x, p_y, p_z$ (GeV)],
    [`eE`], [`float`], [Trigger electron $E$ (GeV)],
    table.cell(colspan: 3)[*`PI0::photons`* --- item 2, one row per selected
      photon],
    [`px`, `py`, `pz`], [`float`], [Photon momentum (GeV)],
    [`E`], [`float`], [Photon energy (GeV), $= |arrow(p)|$],
  ),
  caption: [Slim skim contents, $approx 87$ bytes per event. The photon bank
  carries *no* PID, status, $beta$ or vertex column --- the AI, fiducial
  and $beta$ selections are already applied and cannot be varied
  downstream (@sec:skim-cost). The skim keeps events with $>= 1$ photon
  even though $pi^0$ reconstruction needs $>= 2$, deliberately: single-photon
  events still feed the event-mixing pool.],
) <tab:banks>

Kinematics are stored as `float` and re-widened to `double` downstream, so
the slim and full chains differ at the $approx 10^(-7)$ relative level ---
irrelevant physically, relevant only if the two are diffed bit-for-bit.

= Declared-but-unread configuration <sec:appendix-dead-config>

Keys and fields that appear authoritative but have no effect. None changes
any physics today --- the values coincide with the hard-coded constants ---
but each is a trap for anyone attempting a systematic variation by editing
the configuration.

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    table.header([*Item*], [*Status*]),
    [`DISSelectionConfig::dis_region`],
      [Set to `standard_dis_region()` by both binaries; *never read*. The
       flat `q2_min/max` etc.\ are what apply (@sec:dis).],
    [`[Pi0Reconstruction] min_photon_energy`],
      [`select_photons()` is never called from `process()`; the effective
       floor is the $0.2$ GeV inside the AI purity filter.],
    [`[PhotonSelection] min_photon_energy`],
      [Unused when `use_ai = true` (the production setting); the AI branch
       uses a hard-coded $0.2$ GeV. Same value, different source.],
    [`pcal_lv_min`, `pcal_lw_min`, `dc_r{1,2,3}_edge_min`],
      [Declared on the electron algorithm and set in six TOMLs; *never
       read*. Cuts come from `ElectronCutsService` constants
       ($9.0\/9.0$; $1.680\/2.0\/8.750$).],
    [`[Pi0Analysis] W_edges`, `nu_edges`],
      [Declared and registered but never read --- the mixing pool key
       dropped $W$ and $nu$ (@sec:mixing). The configuration comment still
       describes a five-dimensional pool key. Stale.],
    [`[IO]` section],
      [All five keys read by nothing; the framework prints a
       "does not match any registered algorithm" warning on every run. The
       ifarm/XRootD tuning advice in the comments is not functional.],
    [`[analysis] beam_energy`],
      [Does not reach the @DIS kinematics, which use the hard-coded
       `BeamConfig` default. `[Pi0Analysis] beam_energy` *does* set the
       $gamma^*$ frame. Three independent beam energies that agree only
       because all default to $10.53$.],
    [`enable_*_cut = false`],
      [Does *not* skip the cut --- it makes the pass flag `false` and
       filters *every* event. Counter-intuitive; do not use it to disable a
       cut.],
    [`phi_to_sector`],
      [*Fixed.* It was declared to take radians while the call site passed
       `Particle::phi()`, which returns *degrees*; the $[0,5]$ clamp then
       saturated and collapsed six sectors to $tilde 2$. It now takes
       degrees, matching `Particle::phi()`, and wraps with `fmod` so any
       finite input is handled. The fix also corrected a second, independent
       error: truncation toward zero had put the wedges at $[0,60)$ rather
       than centring them on $0, 60, ... 300$ as
       `VzCorrector::kSectorPhiCenter` requires, making sector 0 twice as
       wide as its neighbours. Both are covered by unit tests. This only
       ever affected the *enhanced* mixing-pool cross-check --- not the
       baseline pool, and no physics observable in this note.],
    [Cutflow labels],
      [`"Momentum > 0.8 GeV"` where the cut is $2.0$ GeV;
       `"DC Chi2/NDF"` where the cut is a DC *edge* distance. Do not
       transcribe a cutflow table from the logs.],
  ),
  caption: [Configuration and code items that do not do what they appear to
  do.],
) <tab:dead-config>

= Unbinned maximum-likelihood method <sec:appendix-mlm>

Documented for reference; *not currently used* --- the production @BSA uses
the binned $chi^2$ fit of @eq:bsa-model.

As an alternative to binning in $phi_h$, the amplitudes can be extracted
from an unbinned likelihood using every event individually. This removes
any dependence on the $phi_h$ binning and makes optimal use of the
statistics --- an advantage in sparsely populated four-dimensional leaves.

For a longitudinally polarised beam on an unpolarised target, the
probability for event $i$ with azimuth $phi_(h,i)$ and helicity
$h_i = plus.minus 1$ is

$ P(h_i | phi_(h,i); A, B, C) = 1/2 (1 + h_i P A_"LU" (phi_(h,i))) , $ <eq:mlm-prob>

with

$ A_"LU" (phi_h) = (A sin phi_h) / (1 + B cos phi_h + C cos 2 phi_h) , $ <eq:mlm-alu>

identifying $A equiv A_"LU"^(sin phi_h)$, $B equiv A_"LU"^(cos phi_h)$,
$C equiv A_"LU"^(cos 2 phi_h)$. The likelihood over $N$ events in a leaf is
the product of @eq:mlm-prob, so the amplitudes follow from minimising

$ -ln cal(L)(A, B, C) = - sum_(i=1)^N ln (1 + h_i P A_"LU" (phi_(h,i))) , $ <eq:mlm-nll>

dropping the constant $-N ln 2$. Uncertainties come from the Hessian at the
minimum.

Because the unpolarised $phi_h$-dependent terms appear identically for both
helicity states, the acceptance largely cancels in the helicity-difference
probability, so the extracted amplitudes are to first order
acceptance-insensitive --- the same argument that protects the binned
$A_"LU"$ (@sec:bsa-definition).

#note-box(title: "Caveat inherited from the binned fit")[
  @eq:mlm-alu is the *three-parameter ratio form*, which
  @sec:bsa-fit records as having been abandoned for exactly the reason that
  matters here: at these per-leaf statistics its denominator develops pole
  pathologies. An unbinned fit does not by itself cure that --- it uses the
  same functional form. The likelihood approach is worth revisiting
  *together with* a decision about whether $B$ and $C$ should be fitted at
  all, or fixed externally.
]

= Result-file schemas <sec:appendix-schemas>

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    table.header([*File*], [*Columns*]),
    [`pi0_multiplicity_ratio_<T>.csv` \ (1450 rows)],
      [`leaf`; `q2, xb, z, pt2` (*geometric box midpoints* ---
       @sec:binning-caveat); `y_LD2, ye_LD2`; `y_<T>, ye_<T>`;
       `N_DIS_LD2, N_DIS_<T>`; `R, R_err`;
       `mu_LD2, sigma_LD2, chi2ndf_LD2`; `mu_<T>, sigma_<T>, chi2ndf_<T>`],
    [`pi0_pT_broadening.csv` \ (320 rows)],
      [`leaf`; `q2_center, xb_center, z_center`; `mean_pT2_<T>, err_<T>`;
       `dpT2_<T>, dpT2_err_<T>`],
    [`pi0_BSA.csv` \ (237 rows)],
      [`leaf`; `{Q2,xB,z,pT2}_{lo,hi}` (kd-tree box);
       `q2_mean, xb_mean, pt2_mean, z_mean` (*count-weighted* --- the
       pattern the other files should adopt); then per target
       `A_<T>, A_err_<T>, chi2_ndof_<T>, ndof_<T>, n_tot_<T>`],
  ),
  caption: [Current result-file schemas. Note the @BSA file already carries
  count-weighted kinematic means while the multiplicity file carries box
  midpoints --- the inconsistency at the heart of @sec:binning-caveat.],
) <tab:schemas>

#warning-box(title: "Superseded files on disk")[
  Older, *broken* result files sit at the top level of `results/` alongside
  the current ones, distinguishable by a `bin_id` (not `leaf`) column and a
  `pt` (not `pt2`) column, and by row counts near $800$ rather than $1450$.
  They are from the product-grid era and are numerically pathological ---
  $R$ reaching $7 times 10^18$, $R_"err"$ to $2 times 10^35$. *Do not read
  them.* Only the files inside the per-observable subdirectories are
  current. Several plots dated before the 16 June CSVs are likewise stale
  and should be regenerated before use.
]
