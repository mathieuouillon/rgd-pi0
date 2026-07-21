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
  Outbending block is used by the results in this note
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
  parameter set, not the $z$ ordering: the `hi` peak is the upstream
  (more negative) foil. Sn is the $-2.8$ peak, Cu the $-7.8$ peak.],
) <tab:vz-inbending>

= Slim schema <sec:appendix-banks>

Stage A writes a single ROOT `TTree` named `events`, one row per event, read
downstream by `uproot` without ROOT.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, center, left),
    table.header([*Branch*], [*Type*], [*Meaning*]),
    table.cell(colspan: 3)[Per event --- scalar branches],
    [`run`], [`int`], [Run number, from `RUN::config`],
    [`helicity`], [`int`], [Beam helicity, HWP-corrected, $plus.minus 1$
      ($0$ = undefined)],
    [`q2`, `xb`, `nu`, `w`, `y`], [`double`], [@DIS kinematics
      (@eq:Q2 -- @eq:xB)],
    [`ex`, `ey`, `ez`, `ee`], [`double`], [Trigger electron $p_x, p_y, p_z, E$
      (GeV)],
    table.cell(colspan: 3)[Per photon --- variable-length array branches],
    [`gpx`, `gpy`, `gpz`], [`double[]`], [Selected photon momenta (GeV)],
    [`g_e_gamma_deg`], [`double[]`], [Each photon's angle to the scattered
      electron (deg)],
  ),
  caption: [The slim `events` tree. It carries no per-photon PID, status,
  $beta$ or vertex --- the @GBT, fiducial and $beta$ selections are already
  applied and cannot be varied downstream (@sec:skim-cost). The skim keeps
  events with $>= 1$ photon even though $pi^0$ reconstruction needs $>= 2$,
  deliberately: single-photon events still feed the event-mixing pool. The
  file also carries a `provenance` directory and one `provenance_stageA` block
  stamping the config hash, target, run and @GBT model (@sec:provenance-gaps).],
) <tab:banks>

Stage A is single-threaded and streaming, and its output is byte-identical
regardless of the `--threads` value.

= Unbinned maximum-likelihood method <sec:appendix-mlm>

Documented for reference; not currently used --- the production @BSA uses
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
  @eq:mlm-alu is the three-parameter ratio form, which
  @sec:bsa-fit records as having been abandoned for exactly the reason that
  matters here: at these per-leaf statistics its denominator develops pole
  pathologies. An unbinned fit does not by itself cure that --- it uses the
  same functional form. The likelihood approach is worth revisiting
  together with a decision about whether $B$ and $C$ should be fitted at
  all, or fixed externally.
]

= Output schemas <sec:appendix-schemas>

Stage B writes four flat trees, all decodable by their index columns without
ROOT.

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    table.header([*Tree*], [*Columns (dense unless noted)*]),
    [`spectra`],
      [`bin4d`, `imgg`, `n_same`, `n_mixed`, `sum_q2`, `sum_xb`, `sum_z`,
       `sum_pt2` --- $n_"4D" times n_"mgg"$ rows],
    [`ptb3d`],
      [`bin3d`, `imgg`, `counts`, `sum_pt2`, `sum_pt4` ---
       $n_"3D" times n_"mgg"$ rows],
    [`n_dis`], [`cell_a`, `n_dis` --- one row per Grid A cell],
    [`bsa`],
      [`bin4d`, `imgg`, `iphi`, `helicity`, `counts` --- sparse
       (zero-count cells omitted); decode via the index columns],
  ),
  caption: [Stage B binned-file schema. The `sum_*` columns are the
  count-weighted kinematic sums that let every abscissa be reported at its
  sideband-subtracted mean rather than a box midpoint (@sec:binning-caveat).],
) <tab:schemas>

Stage C writes one CSV per observable and target, each with a `#`-commented
provenance header (config hash, `provenance_hash`, and any blockers) above the
columns:

/ `ratio_<T>.csv`: `bin4d, i_q2, i_xb, i_z, i_pt2, cell_a`; the
  sideband-subtracted means `q2_mean, xb_mean, z_mean, pt2_mean`; the yields
  `Y_A, sY_A, Y_D, sY_D`; `N_DIS_A, N_DIS_D`; and `R, sR`.
/ `broadening_<T>.csv`: `bin3d, i_q2, i_xb, i_z`; `pt2_A, spt2_A, pt2_D,
  spt2_D`; `delta, sdelta`.
/ `bsa_<T>.csv`: per 4D bin, `A, A_err, chi2, ndof, n` --- written once a
  measured polarization is supplied (@sec:bsa-polarization).
