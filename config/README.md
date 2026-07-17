# `config/` — the only place cut values live

`cuts.json` is the **single source of truth** for every cut value in the RG-D
π⁰ → γγ nuclear multiplicity analysis.

**If a cut value appears anywhere in code, that is a bug.** Not a style
problem, not a cleanup task for later — a bug, to be reported and fixed like
any other. This applies to C++ constants, Python literals, TOML keys, default
arguments, docstrings, log lines, and plot labels alike.

## Why this file exists

This project replaces an older analysis whose cut values were scattered across
C++ constants, TOML files, and config keys that were declared, set, and never
read. The result was not merely untidy. It was **wrong in ways nobody noticed
for months**, and every one of those failures had the same shape: two places
claimed to hold a value, and the one people read was not the one that ran.

The `_comment` fields in `cuts.json` record each of these individually. The
pattern is worth stating once, in the aggregate, because it is the entire
argument for this file:

| What the reader saw | What actually ran |
|---|---|
| Cutflow label `"Momentum > 0.8 GeV"` | `p > 2.0` GeV |
| Docstring "±30 MeV around 135 MeV" | `|m_γγ − 0.135| < 0.2` → *m*<sub>γγ</sub> < 0.335 GeV, no lower bound |
| C++ `SIDEBAND_LOW/HIGH = 0.19/0.24` | Python sideband `(0.17, 0.28)` |
| C++ `FIT_RANGE_LOW/HIGH = 0.10/0.16` | Python fit range `[0.08, 0.20]` |
| `dis_region = standard_dis_region()` → Q² ∈ [1,20], W ∈ [2,25] | flat `q2_min/w_min`; the region field was never read |
| TOML `pcal_lv_min`, `dc_r{1,2,3}_edge_min` in six files | `ElectronCutsService` constants; the TOML changed nothing |
| vz parameter file `sigma_hi = 0.385` | service override `0.415` (8% difference) |
| `[Pi0Analysis] W_edges`, `nu_edges` | never read; the pool key had dropped W and ν |
| `[analysis] beam_energy = 10.53` | never reached the DIS kinematics — three independent beam energies agreeing only by coincidence |

Note the last three rows especially. A stale *label* is embarrassing; a stale
*value* that silently differs from the applied one — the vz sigma, the beam
energy — is a physics error waiting for the day the two stop coinciding.

## Rules

1. **Every cut value is read from `cuts.json` at runtime.** No compiled-in
   defaults, no fallback constants, no "the config didn't have it so we used
   the struct default". A missing required key is a **hard error at startup**,
   not a silent default. The old code's struct defaults are exactly how
   `pi0_mass_window` came to have a `0.03` default and a `±30 MeV` docstring
   that every shipped config overrode to `0.2`.

2. **Fail loudly, never substitute.** Where the old code hit a missing case it
   silently fell through to a wrong-but-plausible answer — most damagingly the
   GBT model lookup, which matched no entry for RG-D runs (18305–19131; the
   table tops out at 16772) and quietly returned the RG-A **inbending pass-1**
   model for every RG-D outbending photon, with no warning logged. Hence
   `photon.allow_rga_fallback: false`: when a run matches no model, the program
   must **refuse to run** and name the run. A program that stops is recoverable.
   A program that lies is not.

3. **Comments are part of the file.** `_comment` fields carry provenance and
   the trap attached to each non-obvious number. When you change a value, change
   its comment. A number whose comment no longer describes it has become exactly
   the thing this file exists to prevent.

4. **Don't "clean up" the redundancy.** Some values look duplicated and are not:
   the electron PCAL fiducial (9.0 cm) and the photon one (14.0 cm); `core_k`
   (1.5, locates the peak) and `n_sigma_window` (3.0, integrates it); the CxC
   vertex peaks (−7.887/0.395) and the Cu/Sn ones (−7.861/0.415). These are
   independent measurements that happen to be adjacent. Unifying them is a
   physics change disguised as a refactor.

5. **Absent values are decisions.** `q2_max` and `w_max` are omitted
   deliberately — the old `100/100` were unreachable at 10.53 GeV (kinematic
   maxima ≈ 9.6 GeV² and ≈ 4.6 GeV) and never rejected an event. An unreachable
   bound in a config invites a reader to believe a cut exists where none does.
   Do not restore them. The same goes for the vacuous `y_min = 0.0`.

## Format

Plain JSON — parses with any standard parser, no JSONC, no comment syntax, no
trailing commas. Two conventions a consumer must handle:

- **`null` means unbounded**, never zero and never a large finite sentinel.
  It appears in `extraction.fit_bounds.amp_max` (+∞) and in the open-ended
  `>=4` photon-multiplicity class. Mapping `null` to `0` would silently pin the
  Gaussian amplitude to zero.
- **Keys beginning with `_` are documentation**, not data. A loader should
  ignore them — but a *reviewer* should not.

Bounds are strict (`<`, `>`) unless a comment says otherwise. Three are not,
and each says so: `photon.min_energy_gev` (`>=`), `photon.theta_{min,max}_deg`
(inclusive both ends), and `electron.trigger.status_abs_max` (exclusive upper).

## Verifying a change

`cuts.json` overdetermines several quantities, which makes it self-checking.
After any edit, confirm these still hold — they are stated in the analysis note
and reproduce exactly from the file:

| Derived from `cuts.json` | Must equal |
|---|---|
| Cu window: −7.861 ± 3 × 0.415 | (−9.106, −6.616) cm |
| Sn window: −2.916 ± 3 × 0.370 | (−4.026, −1.806) cm |
| (0.3 − 0.0) / 200 | 1.5 MeV per bin |
| 17.561·exp(−0.756 p) + 1.0 at p = 0, 2 GeV | 18.56°, 4.87° |
| 8 × 7 × 4 pool bins | 224 |
| 0.1349768 + 0.2 | *m*<sub>γγ</sub> < 0.335 GeV |

If an edit breaks one of these, either the edit is wrong or the note needs
updating — resolve it before committing, and never by adjusting only one side.

## Provenance

Values are transcribed from the analysis note at
`/Users/mathieuouillon/Documents/Typst/pi0-analysis`, principally
`sections/03_event_selection.typ` (electron, DIS, photon, vertex),
`04_pi0_reconstruction.typ` (pairing, pair-level cuts),
`06_background.typ` (mixing, sideband, *m*<sub>γγ</sub> histogram) and
`07_multiplicity_ratio.typ` (yield extraction), cross-checked against the
superseded tree at `/Users/mathieuouillon/Documents/tmp/clas-analysis-1`.

**The note is upstream of this file.** Where they disagree, the note wins and
this file is corrected — with one exception, which is the whole point of the
exercise: where the note documents what the *old code did wrong*, this file
encodes the *fix* and says so in its comment. `photon.allow_rga_fallback` is
the model case.

Two values currently in the file come from **neither** source and are flagged
inline as needing confirmation:

- **`mixing.donors_per_bin = 200`** — the note (`tab:pool`) and the old code
  (C++ default, binary, and TOML alike) all say **50**. Nothing in either
  source mentions 200.
- **`mixing.seed_mode = "file_hash"`** — no old-code equivalent. The old mixing
  was fully deterministic and used no RNG, so there was nothing to seed. The
  unseeded-RNG problem the note *does* record was in the kd-tree binning, which
  was removed from the C++ on 14 July 2026.

Both are plausible new design decisions for this project. Neither is an
inherited value, and their `_comment` fields say so. Resolve them and delete
the flags.
