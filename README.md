# RGD_pi0_analysis

Analysis of neutral pion ($\pi^0 \to \gamma\gamma$) production in CLAS12 Run Group D
data.

Photon pairs are reconstructed from HIPO data files, and the $\pi^0$ signal is
extracted by **fitting the $m_{\gamma\gamma}$ spectrum** — the pairing mass window
(`cuts.json`, `pairing.mass_window_gev`) is deliberately wide, so a `GGPair` is a
gamma-gamma pair and *not* necessarily a $\pi^0$.

## Layout

| Path                  | What                                                                       |
| --------------------- | -------------------------------------------------------------------------- |
| `config/`             | **The only place cut values live.** `cuts.json` + the frozen binning grids. A cut value appearing in code is a bug. |
| `src/core/`           | Pure analysis code → static lib `pi0_core`: constants, kinematics, γγ pairing, binning. No HIPO. Unit-testable without a data file. |
| `src/config/`         | `Cuts::load()` — parses `cuts.json`, fails loudly on a missing key.        |
| `src/photonid/`       | The CatBoost photon classifier: 45-feature builder, run→model map, models. **Models are vendored verbatim — never edit.** |
| `src/vertex/`         | The vertex correction and the per-target `v_z` windows.                    |
| `src/selection/`      | Electron and photon selection; the sampling-fraction tables.               |
| `src/stageA_skim/`    | **Stage A** — HIPO → slim TTree.                                           |
| `src/stageB_bin/`     | **Stage B** — slim → binned spectra + kinematic sums. Includes the frozen donor pool. |
| `src/tools/`          | `dump_columns` (RHipoDS column probe), `make_grid` (equal-statistics grids). |
| `python/pi0/`         | **Stage C** — extraction: yields, $R_A$, $\Delta\langle p_T^2\rangle$, $A_{LU}$, QA. |
| `note/`               | The physics analysis note (Typst). `cd note && typst compile main.typ`. Its figures come from `pi0.plots_selection` (cutflow, GBT score, sampling fraction), `pi0.plots_kinematics` (M(γγ), z, pT², φ_h, the binning grid) and `pi0.plots_results` (R_A, Δ⟨pT²⟩). |
| `data/Vz/`            | Vertex-correction parameters. Vendored verbatim.                           |
| `tests/`              | Catch2 v3 unit tests (C++). `python/tests/` for pytest.                    |
| `external/hipo-cpp/`  | **Submodule, pinned. Do not modify.** See below.                           |
| `meson/native-macos.ini` | Native file for macOS + Homebrew. See below.                            |

Includes within the project are project-relative — the build adds `src/` as an
include dir:

```cpp
#include "core/Kinematics.hpp"
```

## Requirements

- C++17 compiler
- meson >= 1.2, ninja
- ROOT >= 6.36 (tested against 6.40.02)
- vdt, fmt, lz4 (`brew install vdt fmt lz4`)

## Build

Two things are required and **both are easy to get wrong**, so read this section
rather than guessing:

```sh
export ROOT_INCLUDE_PATH=/opt/homebrew/include
meson setup build --native-file meson/native-macos.ini
ninja -C build
```

### 1. `--native-file meson/native-macos.ini`

`pkg-config` on macOS here resolves to MacPorts' (`/opt/local/bin/pkg-config`)
with an unset `PKG_CONFIG_PATH`, so Homebrew's `fmt.pc` and `liblz4.pc` are not
found. hipo-cpp declares both deps with a wrap fallback, so instead of failing,
meson **silently downloads** fmt and lz4 into `external/`. The native file points
pkg-config at Homebrew and stops that.

### 2. `export ROOT_INCLUDE_PATH=/opt/homebrew/include`

ROOT 6.40's `<ROOT/RVec.hxx>` includes `<vdt/vdtMath.h>`. The C++ compile is
fine — ROOT's CMake config exports a `VDT::VDT` target that carries the right
`-I`. But **rootcling is not the C++ compiler** and never sees meson's `-I`
flags: hipo-cpp's dictionary target hardcodes its rootcling command line. Without
this variable the build dies generating the `libHipoDataFrame` dictionary:

```
fatal error: 'vdt/vdtMath.h' file not found
```

`ROOT_INCLUDE_PATH` is the mechanism rootcling intends for this, and upstream
hipo-cpp uses it too (`external/hipo-cpp/meson/this_hipo.sh.in`). It is read at
**build** time, so it must be exported in the shell that runs `ninja` — not just
`meson setup`. It cannot be moved into the native file: meson machine files have
no `[env]` section.

`meson.build` checks this at configure time and fails with a pointed message
rather than letting you discover it mid-build.

### Tests

```sh
meson test -C build
```

Catch2 v3 is fetched via `external/catch2.wrap` on first configure (needs
network once). To skip tests entirely:

```sh
meson setup build --native-file meson/native-macos.ini -Dbuild_tests=false
```

### Interactive ROOT

The built programs need no environment (meson bakes the rpath into the build
tree). But an interactive ROOT session that has to load the HipoDataFrame
dictionary does need `ROOT_INCLUDE_PATH`; `meson devenv -C build` exports it for
you:

```sh
meson devenv -C build
```

## Running the analysis

The analysis is four stages. Each is a separate program, each takes exactly one
input, and **each stamps its provenance into its output** — inputs, config hash,
grid hashes, the photon model used, whether a fallback was taken. The next stage
reads that block back and refuses to build physics on something unpublishable.
That refusal is the design, not an obstacle: see [Refusals](#refusals-and-what-they-mean).

```
  HIPO DST                                       one file per run
     │
     │  stageA_skim      electron + DIS + GBT photon selection
     ▼
  slim.root             TTree "events": DIS kinematics + selected photon momenta
     │
     ├─ make_grid       (once) equal-statistics grid edges → config/binning/*.json
     │
     │  stageB_bin       γγ pairing, frozen donor pool, binning
     ▼
  binned.root           per-bin m_γγ spectra (same + mixed), kinematic sums,
     │                  N_DIS, BSA counts
     │  python -m pi0.*  yield extraction, R_A, Δ⟨pT²⟩, A_LU, QA
     ▼
  CSV + figures
```

Each program prints a cutflow and its provenance block on success, and takes no
environment to run.

#### Multi-threading and progress

`stageA_skim`, `stageB_bin` and `make_grid` all take **`--threads N` (default
1)**. The default is the single-threaded path, byte-for-byte the old output, and
**the output does not depend on `N`** — that is the invariant, not an aspiration:

- **Stage B** splits the events into a *fixed* number of partitions decided by
  the event count (never by the thread count), sums each partition's kinematic
  moments with a compensated accumulator, and merges them in partition-index
  order. So `sum_q2` and the other abscissa moments are **bit-identical for any
  `--threads`**, including 1. Verified across `--threads 1/4/8`.
- **Stage A** scores photons (the GBT, the costly part) in parallel but fills the
  slim and QA trees in entry order, so the slim is byte-identical across thread
  counts. This matters because Stage B's donor pool is order-dependent — a
  nondeterministic slim order would make the mixed-event background irreproducible.
- **make_grid** parallelises across input files; its edges are quantiles of the
  pooled sample, so order — and thread count — cannot reach the result.

Two honest caveats:

- **`--threads N>1` buffers events in memory** (Stage A and Stage B), so it is
  for bounded runs (`--max-events`, or a per-file slim), not for pointing at a
  raw 200 GB HIPO. **The farm uses `--threads 1`**, which streams in O(1) memory,
  so production is unaffected.
- **Stage A shows little speedup on a single HIPO file.** Its cost is dominated
  by a one-time serial scan RHipoDS does to decompress the file's record index
  (reading 50 events costs about the same as reading 4000). `--threads` helps
  only when GBT scoring dominates — many photons across many events — and the
  scan is paid once regardless. It never *hurts*, and the farm parallelises
  across files anyway.

Progress goes to **stderr**: a live bar on a terminal, an occasional
grep-friendly line in a batch log (so `log.txt` stays readable). stdout is
reserved for the cutflow and provenance.

---

### Step 0 — `dump_columns` (optional, diagnostic)

Prints the RDataFrame column names RHipoDS exposes for a HIPO file. Use it when a
bank name is in doubt: it is the ground truth for the naming contract the rest of
the project depends on (HIPO `::` and `.` each become `_`, so `REC::Particle.px`
→ `REC_Particle_px`).

```sh
./build/src/tools/dump_columns <file.hipo> [prefix]
```

| Argument | Meaning |
| --- | --- |
| `<file.hipo>` | any HIPO file |
| `[prefix]` | filter, default `REC_`. Pass `""` to list every column (order 1000, depending on the banks the file carries). |

---

### Step 1 — `stageA_skim`: HIPO → slim

Applies the electron selection (including the target's $v_z$ window), the DIS
cuts, and the GBT photon selection. Writes one entry per surviving DIS event.

```sh
./build/src/stageA_skim/stageA_skim \
    --input  /path/to/rec_clas_018500.hipo \
    --output slim/LD2_018500.root \
    --config config/cuts.json \
    --target LD2
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--input <file.hipo>` | yes | **One** HIPO file. Entry numbering restarts per file — do not chain. |
| `--output <slim.root>` | yes | Slim ROOT file to create (overwritten if present). |
| `--config <cuts.json>` | yes | Every threshold comes from here. |
| `--target <LD2\|CxC\|Cu\|Sn>` | yes | Selects the vertex window. **For Cu/Sn this *is* the target assignment** — the two foils are distinguished only by $v_z$. |
| `--run <N>` | no | Override the run used to pick the GBT model. Normally read from `RUN::config`; pass only if the file's run is absent or wrong. |
| `--max-events <N>` | no | Stop after N tag-0 events. For smoke tests. The output is a **prefix**, not a sample — stamped into the provenance so a partial run cannot pass for a full one. |
| `--threads <N>` | no | Parallelise the GBT scoring. Default 1 (the farm path, streaming). Output is byte-identical to `--threads 1` for any N. `N>1` buffers events — bounded runs only. |
| `--partition-target <N>` | no | Approx events per reproducible partition (default 200000). A granularity knob, not a correctness knob. |

**Exit codes:** `0` ok · `1` usage · `2` no run number found · `3` **no GBT model
covers this run** · `4` general · `5` file holds more than one run.

Only tag-0 (physics) events are read; RHipoDS applies that filter and it is
correct. A file may hold many non-tag-0 events that are correctly excluded.

> **Exit 3 is expected on RG-D.** No trained GBT photon model exists for runs
> 18305–19131 — the map stops at 16772. See [the RG-D photon gap](#the-rg-d-photon-gap).

Run one job per HIPO file; the stage is deliberately single-file so that the
donor pool downstream is built per file.

---

### Step 2 — `make_grid`: freeze the binning (once per production)

Computes equal-statistics (quantile) edges from real slim files and writes them to
JSON. **Do this once**, commit the result, and never let it drift: those two files
*are* the binning of the production, and the reason they are in version control is
that the superseded analysis's kd-tree edges lived only on `/work` and are now
unrecoverable.

```sh
./build/src/tools/make_grid/make_grid \
    --input slim/LD2_018500.root --input slim/LD2_018501.root \
    --config config/cuts.json \
    --out-a config/binning/grid_A_q2_xb.json \
    --out-b config/binning/grid_B_z_pt2.json
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--input <slim.root>` | yes | Repeat for more. Edges are quantiles of the pooled sample; order affects only the provenance listing. |
| `--config <cuts.json>` | yes | Grid shape, pairing cuts and the z window come from here. Its hash is stamped into both outputs and checked against what each input recorded at skim time. |
| `--out-a <json>` | yes | Grid A = $(Q^2, x_B)$, per **event**. |
| `--out-b <json>` | yes | Grid B = $(z, p_T^2)$, per **π⁰**. |
| `--na WxH` | no | Override `binning.grid_a`, e.g. `8x7`. See the warning below. |
| `--nb WxH` | no | Override `binning.grid_b`, e.g. `5x5`. Ditto. |
| `--max-events <N>` | no | Stop after N events across all inputs. A **prefix**, not a random subset — the edges are then not the edges of the full dataset. Stamped into the output. |
| `--threads <N>` | no | Read input files in parallel. Default 1. Edges are identical for any N (they are quantiles of the pooled sample). `--max-events` forces single-threaded so truncation stays a deterministic prefix. |
| `--cap-z <v>` | no | Clamp Grid B's top $z$ edge and **discard** π⁰ above it. |
| `--cap-pt2 <v>` | no | The same for $p_T^2$, in GeV². |

> **`--na`/`--nb` do not write back to `cuts.json`.** They change the grid the C++
> emits but leave the config's shape saying something else, and `pi0.config` then
> refuses the pair outright — *"cuts.json binning/grid_a/n_q2 = 8 but the grid file
> gives 2. The configuration is internally inconsistent; fix it, do not guess."*
> To change the shape for a production, **edit `cuts.json` and re-run `make_grid`
> without the overrides**. The flags are for quick experiments only.

> **The caps default to off and should usually stay off.** With the
> count-weighted abscissae Stage B accumulates, a wide outer bin is honestly
> *positioned* even though it is coarse — so capping throws away data to solve a
> problem that is already solved. See `config/cuts.json` `/binning/_cap_comment`.

Grid B is computed **once over all π⁰**, not per Grid A cell. That factorization
is a deliberate simplification against the old nested kd-tree, and it is recorded
in the output JSON.

---

### Step 3 — `stageB_bin`: slim → binned

Two passes over the slim file. Pass 0 builds the **frozen donor pool** with a
seeded pre-pass; pass 1 bins, with mixing reduced to a stateless lookup against a
`const` pool. That is what makes the background bit-for-bit reproducible.

```sh
./build/src/stageB_bin/stageB_bin \
    --input  slim/LD2_018500.root \
    --input  slim/LD2_018501.root \
    --output binned/LD2.root \
    --config config/cuts.json \
    --grid-a config/binning/grid_A_q2_xb.json \
    --grid-b config/binning/grid_B_z_pt2.json
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--input <slim.root>` | yes | A Stage A file (TTree `events`). **Repeatable** — see below. |
| `--output <binned.root>` | yes | Binned ROOT file to create. |
| `--config <cuts.json>` | yes | Every threshold and every axis. |
| `--grid-a <json>` | yes | The frozen $(Q^2,x_B)$ grid. Its hash is stamped into the output. |
| `--grid-b <json>` | yes | The frozen $(z,p_T^2)$ grid. Ditto. |
| `--max-events <N>` | no | Stop after N entries **in total across all inputs**, not per input. **Both passes are truncated identically** — a pool built from more events than were binned would estimate the background from a sample the spectra never saw. |
| `--threads <N>` | no | Parallelise pass 1. Default 1. The four output trees are byte-identical for any N — a fixed-partition compensated reduction, independent of thread count. |
| `--partition-target <N>` | no | Approx events per partition (default 200000). Granularity/tuning only; does not affect the output. |
| `--allow-truncated-inputs` | no | Accept an input that Stage A itself marked truncated. Off by default; see below. |

**Give a target all of its slims in one run.** The donor pool is a reservoir per
$(Q^2, x_B, N_\gamma)$ bin, and a bin no input filled has **no** mixed background
at all — not a thin one, a missing one. One farm slim fills **33 of the 224**
shipped pool bins; two fill 56. Running Stage B once per file and combining
afterwards is not an option anyway: `pi0.io.load` requires dense rows, so `hadd`
of two binned files is rejected.

Three properties make chaining safe, all of them checked:

- **Order does not matter.** The chain is sorted by **content SHA-256** before
  anything reads it, so a shell glob is safe: `--input a --input b` and
  `--input b --input a` give bit-identical outputs. Sorting by *path* would have
  made the pool depend on where the files live — copy a target's slims to another
  mount and the offer sequence, and so the background, would move. The resolved
  order and the seed are stamped into the output.
- **A mixed chain is refused, not warned about.** Inputs must agree on target,
  config hash, GBT model, beam energy and polarity. Nothing downstream of the
  read knows a photon's target — a donor is four floats — so chaining LD₂ with Sn
  would make one target's photons the mixed background subtracted from the
  other's events, corrupting both halves of $R_A$ while the spectra looked
  entirely ordinary.
- **A truncated input is refused.** Stage A marks its own output when
  `--max-events` was used. Chained with full files, its short N_DIS — the $R_A$
  normalisation denominator — would be invisible.

**Output** — four flat TTrees, readable by uproot with no ROOT and no dictionary:

| Tree | Rows | Contents |
| --- | --- | --- |
| `spectra` | dense, 1400 × 200 | `bin4d, imgg, n_same, n_mixed, sum_q2, sum_xb, sum_z, sum_pt2` — the kinematic sums are **same-event only**. |
| `ptb3d` | dense, 280 × 200 | `bin3d, imgg, counts, sum_pt2, sum_pt4` — binned in $m_{\gamma\gamma}$ so $\Delta\langle p_T^2\rangle$ can be sideband-subtracted. |
| `n_dis` | 56 | `cell_a, n_dis` — filled **once per event**, not per π⁰. The $R_A$ denominator. |
| `bsa` | **sparse** | `bin4d, imgg, iphi, helicity, counts`. Decode **via the index columns**; a positional read is silently wrong. |

Index decoding (also implemented in `python/pi0/config.py`, the only place that
knows it):

```
cell_a = iq2 * n_xb  + ixb           bin4d = cell_a * (n_z * n_pt2) + cell_b
cell_b = iz  * n_pt2 + ipt2          bin3d = cell_a * n_z + iz
m_γγ bin centre = mgg_histogram.min_gev + (imgg + 0.5) * bin_width
```

The axis sizes come from `cuts.json` — `binning.grid_a`, `binning.grid_b`,
`mgg_histogram`. With the shipped 8×7 and 5×5 that is 56 Grid A cells, 1400 4D
bins, 280 3D bins and 200 mass bins over [0, 0.3] GeV. **Read them from the
config, don't hard-code them**; changing the shape changes every index above.

The `sum_*` columns are what make the count-weighted abscissae work: each 4D bin
carries its own $\langle Q^2\rangle, \langle x_B\rangle, \langle z\rangle,
\langle p_T^2\rangle$ per mass bin, so the reported abscissa is
sideband-subtracted exactly like the yield rather than being the geometric centre
of the box. In the superseded analysis 52.6% of leaves were reported at a centre
that was not where the data was, and 6.9% at a point no event could occupy.

---

### Step 4 — `python -m pi0.*`: extraction

The Python stage reads the Stage B ROOT files with **uproot** — no ROOT, no C++,
no dictionaries. It runs anywhere, including a laptop with none of Step 1–3's
dependencies installed.

```sh
cd python
python -m venv .venv && . .venv/bin/activate
pip install -e .            # or: pip install -e '.[test]' && pytest
```

That also installs four console scripts — `pi0-qa`, `pi0-ratio`,
`pi0-broadening`, `pi0-bsa` — equivalent to the `python -m pi0.*` forms used
below. With `uv`, `uv run python -m pi0.<module> ...` needs no venv step.

> **`--config` takes the config *directory* here**, not the file — the Python
> stage reads `cuts.json` *and* both grid JSONs and checks their hashes against
> what the input recorded. The C++ programs take the file. And `pi0.bsa` takes
> `--file`, while `pi0.ratio` and `pi0.broadening` take `--target` / `--ld2`.

#### `pi0.qa` — diagnostics (start here)

```sh
python -m pi0.qa --file ../binned/LD2_018500.root --config ../config --outdir figs/
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--file` | yes | A Stage B file. |
| `--config` | yes | The `config/` **directory**. |
| `--outdir` | yes | Where the figures go. |
| `--allow-unpublishable` | no | Proceed despite fatal provenance blockers. Figures get a **watermark**. |
| `--shared-window` | no | Locate every bin's mass window with the fit to the *summed* spectrum instead of its own. **Not the production path** — a fallback for files whose per-bin statistics cannot support a fit. Labelled on the output. |

The figure worth looking at is the **abscissa-vs-bin-centre** plot: it shows what
the count-weighted mean bought over the geometric centre.

#### `pi0.ratio` — the multiplicity ratio $R_A$

```sh
python -m pi0.ratio --target ../binned/Sn.root --ld2 ../binned/LD2.root \
                    --config ../config --out ra_Sn.csv
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--target` | yes | Stage B file for the **nuclear** target. |
| `--ld2` | yes | Stage B file for **LD2** — always the denominator. |
| `--config` | yes | The `config/` directory. |
| `--out` | yes | Output CSV. |
| `--allow-unpublishable` | no | As above. Output is diagnostic only. |

Needs **two different targets**. Passing the same file twice is refused —
dividing a target by itself measures nothing.

#### `pi0.broadening` — $\Delta\langle p_T^2\rangle$

```sh
python -m pi0.broadening --target ../binned/Sn.root --ld2 ../binned/LD2.root \
                         --config ../config --out dpt2_Sn.csv
```

Same options as `pi0.ratio`. **Sideband-subtracted**, unlike the superseded
analysis, whose moments were accumulated over a ±200 MeV window with no
subtraction and are therefore diluted by an unknown factor.

#### `pi0.bsa` — $A_{LU}$

```sh
python -m pi0.bsa --file ../binned/LD2.root --config ../config --out bsa.csv \
                  --polarization 0.86 --polarization-err 0.03
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--file` | yes | A Stage B file. |
| `--config` | yes | The `config/` directory. |
| `--out` | yes | Output CSV. |
| `--polarization P` | **yes*** | The measured beam polarization. **No default.** |
| `--polarization-err σ_P` | **yes*** | Its uncertainty — a fully-correlated scale error on every $A$. |
| `--allow-uncorrected-dilution` | no | Emit bins whose $S/(S+B)$ could not be measured. Those $A_{LU}$ are then **lower bounds**. Off by default; such bins are dropped. |
| `--allow-unpublishable` | no | As above. |

\* unless `cuts.json` `/bsa/polarization/value` is set. It is `null` on purpose.
`pi0.bsa` **refuses to run** without a polarization rather than substituting one —
the old code's `0.85` was a self-declared placeholder, and every $A_{LU}$ it
published scales by $0.85/P_\text{true}$.

---

## Running it at scale — the farm and `local_batch`

Steps 1–4 above are the pipeline for one file. LD₂ alone is 30,866 HIPO files
spanning 160 runs — of which **135 are the outbending production** (see the run
lists below). So there are two drivers:

| | | |
| --- | --- | --- |
| `python -m pi0.batch` | JLab farm, **SWIF2** | fans Stage A out over ~28k jobs |
| `python -m pi0.local_batch` | one machine, N subprocesses | the same selection, no scheduler |

Both are **pure standard library** — no uproot, no numpy, no ROOT, no build.
That is deliberate: you submit from an ifarm login node, which has no business
carrying the extraction stage's dependencies.

### The run lists — read this before your first submission

`config/runs.json` holds the RG-D run lists, transcribed from clas-framework's
`clas12/Runs.hpp`. Both drivers **refuse any file whose run is not in the
requested (polarity, target) list**.

This is not bookkeeping. Inbending and outbending runs live in the **same `/mss`
directory tree** — LD₂ inbending is 18305–18336, outbending starts at 18419 — and
**the slim schema does not record polarity**, so a plain recursive scan mixes
torus polarities into one dataset and nothing downstream can notice. The note
(`tab:runs`) records that the original lists were *"not applied by any binary"*;
this is what applies them.

Measured against clas-framework's real LD₂ file list (`ld2_files.dat`, 30,866
files, named as the input in `swif2_pi0_skim.py`'s own docstring):

```
files accepted : 27971
files rejected :  2895   (9.4%)
                  2739   all 24 inbending LD2 runs — wrong torus polarity
                   156   run 18432 — in NEITHER list: it is
                         LD2_trigger_rgd_v2_2_Q2_2_5, a Q2 > 2.5 trigger config
```

Run 18432 is the one to worry about: a Q²-biased trigger sculpts the Q² spectrum,
and **Q² is a Grid A binning axis**.

Four targets, three run lists — Cu and Sn are two foils in one assembly exposed
in the *same* CuSn runs and separated only by the electron `vz` window, so
`Cu.farm.json` and `Sn.farm.json` point at an identical `/mss` tree and differ
only in `target`.

### Farm config

One per target, in `config/farm/`. It carries **no cut values** — only where the
data is and what resources to ask for. A threshold here is a bug (there is a test
asserting it).

```jsonc
{ "farm": {
    "target": "LD2",              // LD2 | CxC | Cu | Sn  (NOT CuSn — that is a run list)
    "polarity": "outbending",     // selects the run list; flipping it selects inbending
    "inputs": ["/mss/clas12/rg-d/production/pass1/recon/LD2/dst/recon/"],
    "files_per_job": 1,
    "exclude_runs": [],
    "output_dir": "/volatile/clas12/ouillon/rgd-pi0/LD2",
    "log_dir":    "/volatile/clas12/ouillon/rgd-pi0/LD2/logs",
    "swif2": {
      "workflow": "rgd_pi0_stageA_LD2_outbending",
      "account": "clas12", "partition": "production",
      "cores": 1, "ram_gb": 4, "disk_gb": 20, "time": "04:00:00",
      "modules": { "load": ["gcc/13.2.0", "root/6.36.04"] },
      "tags": { "stream": "pi0", "pass": "pass1" }
    },
    "env": {}
} }
```

Two nesting details are inherited from clas-framework's *code* (its struct layout
suggests otherwise, and copying that gets them backwards): **`env` is a child of
`farm`**, not of `farm.swif2`; and **`modules` is an object with a `load` array**,
not a bare array. A bare array parses fine and silently loads nothing, so the job
fails on the node for want of ROOT rather than at submit — `load_farm_config`
refuses it instead.

`time` is converted to SWIF2's seconds form (`04:00:00` → `14400s`) at synthesis,
and a malformed value is **refused at parse time**. clas-framework returns
malformed input unchanged and lets SWIF2 reject it after you have submitted.

### Farm: Stage A

```sh
python -m pi0.batch config/farm/LD2.farm.json                    # dry-run
python -m pi0.batch config/farm/LD2.farm.json --submit
```

| Option | Meaning |
| --- | --- |
| `--config <dir>` | config **directory** (default `config`). Needs `cuts.json` + `runs.json`. |
| `--submit` | actually run the generated script. Default is dry-run. |
| `--stage a\|b` | which stage to drive (default `a`). |
| `--file-list <f>` | take the file list from a `.dat` instead of scanning. **Prefer this on a login node** — scanning 30k tape paths is slow, and the list is then a record of exactly what the production consumed. |
| `--files-per-job N` | override `farm/files_per_job`. |
| `--workflow NAME` | override `farm/swif2/workflow`. |
| `--max-jobs N` | synthesize at most N jobs. Labelled a test submission in the output. |
| `--scripts-dir <d>` | where wrappers + the swif2 script go (default `batch_scripts`). |
| `--exe <path>` | `stageA_skim` **as seen on the node** — must be absolute and node-visible. See below. |
| `--allow-rga-fallback-production` | submit despite the pre-flight blocker below. |

The dry-run writes the wrapper scripts and the swif2 script and **nothing else** —
unlike clas-framework's, whose dry-run performs its full snapshot copy (hundreds
of MB, clobbering the previous workflow's frozen program dir) before it ever
checks `--submit`.

**`--exe` must be absolute and on a filesystem the nodes mount** (`/work`,
`/volatile`, `/cache`, `/home`, `/group`, `/scigroup`, `/u`). A SWIF2 job runs in
a scratch directory containing *only* the staged inputs, so a relative path like
`./build/src/stageA_skim/stageA_skim` resolves to nothing there and **every job
exits 127**. `pi0.batch` refuses both cases up front and offers you the resolved
path. So: build on `/work`, and

```sh
python -m pi0.batch config/farm/LD2.farm.json \
    --exe /work/clas12/users/<you>/rgd-pi0/build/src/stageA_skim/stageA_skim --submit
```

No libraries are staged, deliberately: meson bakes an rpath into the build tree,
so the binary runs with no `LD_LIBRARY_PATH` provided the tree itself is
node-visible. (clas-framework's `swif2_pi0_skim.py` snapshots the exe and every
`*.so` into a frozen directory — necessary there because it copies by basename
into a flat dir, which loses the rpath.)

**The pre-flight is the most valuable thing here.** On RG-D, `stageA_skim` exits 3
for *every* run, so a production submitted with `photon.allow_rga_fallback = false`
is ~28,000 identical failures:

```
============================================================================
REFUSING TO SUBMIT
============================================================================
  * photon.allow_rga_fallback is false, and NO GBT PHOTON MODEL COVERS RG-D.
      The model map stops at run 16772; RG-D is 18305-19131. stageA_skim will
      exit 3 on every one of these 27971 jobs. ...
```

With the fallback on, it submits — and says so in a way you cannot miss, because
every output is then stamped `gbt.fallback_used = TRUE` and the Python stage will
refuse to publish from it.

Each job runs the skim **once per file** (`stageA_skim` takes one file per
invocation) and copies back `slim_<i>.root` plus a log. The wrapper deliberately
does **not** use `set -e`: one bad file must not abandon the rest of the chunk
with a zero exit.

#### Only tape is staged

`-input` **copies**. Whether that is right depends entirely on where the file is:

| Input | What happens | Why |
| --- | --- | --- |
| `/mss/…` | staged as `input_<j>.hipo` | it is on **tape** and must be brought to disk |
| `/cache`, `/work`, `/volatile` | **not staged** — the skim opens the path directly | already a mounted cluster filesystem the worker reads |

This is not a micro-optimisation. An RG-D train skim is **~200 GB**
(`SIDIS_018431.hipo`, measured), so staging one puts a 202 GB copy into a job
that asked for 20 GB of disk and reads 2000 events. The first real submission did
exactly that and both jobs sat in `preparing` — SWIF2's staging state — where
they would have stayed:

```
site_job_disk_bytes = 20000000000     ← 20 GB requested
local               = input_0.hipo    ← a 202 GB file being copied in
job_attempt_status  = preparing
```

`pi0.farm.job_inputs()` is the single place that decides, and both the swif2
script and the wrapper call it, so they cannot disagree about which name a file
has. `disk_gb` therefore only needs to cover the **outputs** unless you are
skimming from `/mss`.

#### Fast turnaround for a smoke test

`priority` is the fast queue; production runs belong on `production`. Both are
CLI overrides rather than config edits:

```sh
python -m pi0.batch config/farm/LD2.farm.json --config config \
    --exe /work/clas12b/users/<you>/rgd-pi0/build/src/stageA_skim/stageA_skim \
    --max-jobs 2 --max-events 2000 \
    --partition priority --time 01:00:00 --submit
```

Monitor and retry the usual way:

```sh
swif2 status     -workflow rgd_pi0_stageA_LD2_outbending
swif2 diagnose   -workflow rgd_pi0_stageA_LD2_outbending   # why jobs are not progressing
swif2 show-job   -workflow rgd_pi0_stageA_LD2_outbending -name <job>
swif2 retry-jobs -workflow rgd_pi0_stageA_LD2_outbending   # resubmit problem jobs
swif2 cancel     -workflow rgd_pi0_stageA_LD2_outbending -delete
```

These are checked against `swif2 -help` on ifarm. There is no `list-jobs` —
clas-framework's docs name one and it errors with *"is not a valid swif
command"*; the real listing command is `swif2 list`.

### Farm: Stage B

```sh
python -m pi0.batch config/farm/LD2.farm.json --stage b          # print the command
python -m pi0.batch config/farm/LD2.farm.json --stage b --submit  # run it
```

Stage B is **not** a SWIF2 workflow, on purpose. It is one job per target, it
needs *every* Stage A output, and it runs **once over all of them** so the donor
pool is drawn from the whole target rather than one file at a time. There is
nothing for a scheduler to overlap, and expressing "wait for all of A" as SWIF2
antecedents would mean one job with ~28,000 of them. (clas-framework has no
antecedent support at all — its `JobSpec` has four fields and emits a flat
create/add-job/run.)

`stageB_bin` takes **repeated `--input`** for this, and validates that every input
agrees on target, config hash and photon model. Chaining an LD₂ slim into an Sn
run would put Sn photons into LD₂'s mixed background and quietly corrupt both
halves of R_A — so it refuses.

### local_batch

```sh
python -m pi0.local_batch --farm config/farm/LD2.farm.json --concurrent 8
python -m pi0.local_batch --input /path/to/hipo --target LD2 --concurrent 8
```

| Option | Meaning |
| --- | --- |
| `--farm <f>` | a farm JSON: gives inputs, target, and the run-list filter. |
| `--input <p>` | a directory or one file, instead of `--farm`. Needs `--target`. |
| `--target`, `--polarity` | with `--input`. |
| `--no-run-filter` | skip the run-list filter. **Only** for non-production data (a test file, one run you are debugging). |
| `--concurrent N` | how many skims run at once (default: half your cores). |
| `--outdir <d>` | slims land at `<d>/<run>/slim_<stem>.root`, logs at `<d>/logs/`. |
| `--max-events N` | per file. A prefix, not a sample — stamped into every output. |
| `--max-files N` | process at most N files. |
| `--dry-run` | plan and print; run nothing. |
| `--stage-b` | after Stage A, run the single Stage B over every slim produced. |

**There is no `--slices`.** clas-framework's `local_batch` slices one big HIPO
file into record ranges because its analysis binary honours `--record-range`.
`stageA_skim` has no such flag and takes one file per invocation *by contract*
(it exits 5 on a multi-run file, because its provenance header records a single
run and a single model, and recording one for a file holding several would be a
lie). With ~30k files the file *is* the natural unit, so `--concurrent` is the
only parallelism knob and it means what it says.

Two things this does that the tool it mirrors does not:

- **A failure names the file, the exit code, what that code means, and the log.**
  clas-framework's prints `Slices: 100 (99 OK, 1 failed)` and stops — you grep
  100 logs by hand.
  ```
  FAILURES, by exit code:
    exit 3 -- NO GBT PHOTON MODEL COVERS THIS RUN -- expected on RG-D unless
              photon.allow_rga_fallback is true   (1 file(s))
        /path/rec_clas_022083.evio.00000-00009.hipo
          log: slim/logs/rec_clas_022083.evio.00000-00009.log
  ```
- **Ctrl-C kills the children.** clas-framework's has no signal handling; its
  children die only because the terminal happens to deliver SIGINT to the whole
  foreground process group, which stops being true the moment anything puts them
  in their own group.

`--stage-b` **refuses if any Stage A job failed**: the slim set would be
incomplete, and N_DIS — the R_A normalisation denominator — would be short by
those files with nothing in the output recording it.

## Refusals, and what they mean

The pipeline refuses rather than producing a plausible-looking wrong number. If
something refuses, that is the tool working. The refusals, and the fix for each:

| Refusal | Why | What to do |
| --- | --- | --- |
| `stageA_skim` exit **3**, "no GBT photon model covers this file" | The run is outside every trained model's range. RG-D is. | See below. Set `photon.allow_rga_fallback: true` **only** for a study you will report as such. |
| `ProvenanceError: N fatal provenance blocker(s)` | The file was made with placeholder grids, a fallback model, a drifted config, or was truncated. | Fix the cause. `--allow-unpublishable` forces it, watermarks every figure, and the numbers are then diagnostic only. |
| "both files carry `target='LD2'`" | $R_A$ needs two different targets. | Pass a nuclear target and LD2. |
| "No beam polarization" | There is no measured RG-D Møller number. | `--polarization P --polarization-err σ_P`, or fill the config key. |
| `below_min_stats` on every bin | Too few π⁰ for per-bin fits. | More data. `--shared-window` (QA only) will show you the shape meanwhile. |
| `negative_variance` | Subtraction drove $\langle p_T^4\rangle < \langle p_T^2\rangle^2$. | Low statistics. The gate rejects it rather than emit an imaginary error bar. |

### The RG-D photon gap

**No trained GBT photon model exists for RG-D** (runs 18305–19131). The model map
stops at run 16772 — this is an upstream CLAS12 gap, not a local one, and it is
identical in Iguana. The superseded analysis fell back to an **RG-A inbending
pass-1** model *silently*; upstream at least logs a warning.

Here the fallback is **opt-in, defaults to refusing, and is stamped into the
provenance** when taken. Applying an RG-A inbending model to RG-D outbending
nuclear-target data is a real assumption: the feature vector includes neighbour
multiplicity, and calorimeter occupancy differs between LD₂ and Sn — precisely
the direction in which photon efficiency would fail to cancel in a target ratio.

Getting a real RG-D model is the single highest-value fix available to the photon
selection. See `src/photonid/models/PROVENANCE.md`.

## Do not modify `external/hipo-cpp`

`external/hipo-cpp` is a **pristine git checkout of hipo-cpp v4.4.1**
(code.jlab.org/hallb/clas12/hipo-cpp). Never edit it, and never commit changes
inside it — `git -C external/hipo-cpp status` must stay clean.

It is consumed read-only as a meson subproject. Everything the top-level build
needs is pulled out with `get_variable()`:

- `hipo_dep` — hipo4/meson.build's own `declare_dependency()` (libhipo4).
- `dataframe_lib` + `dataframe_headers` — the dataframes extension exposes only a
  raw `library()` object and a header list; it has **no** `declare_dependency()`,
  so the top-level `meson.build` assembles `HipoDataFrame_dep` from them.

If you think you need to patch hipo-cpp, add the workaround to the top-level
`meson.build` instead, or push the fix upstream.

## Conventions

- Namespace `pi0` for everything. 4-space indent, no tabs.
- Headers `.hpp`, sources `.cpp`, ROOT-dependent programs `.cxx`.
- Prefer `double`.
- **Angles are in DEGREES at any public boundary.** If a function takes radians
  its parameter *must* be named `*_rad` and say so in its doc comment. (The
  analysis this replaces shipped two separate degree/radian bugs. Be explicit.)

## Licence

**LGPL v3.0** — see [`COPYING.LESSER`](COPYING.LESSER) (the LGPL v3 supplement)
and [`COPYING`](COPYING) (the GPL v3 text it builds on). LGPL v3 consists of
both texts together.

This is not a free choice. `src/photonid/models/` contains five generated
CatBoost models that are byte-for-byte CLAS12/Iguana products distributed under
LGPL v3, and they are incorporated into this work rather than linked as a
separate library — so the combined work is distributed under the same terms.
See [`NOTICE`](NOTICE) for the full attribution of every vendored component, and
[`src/photonid/models/PROVENANCE.md`](src/photonid/models/PROVENANCE.md) for the
models specifically.

## Status

Pre-publication. **Nothing produced by this code is quotable yet**, and the
outputs say so themselves — every file carries a provenance block recording its
inputs, config hash, grid hashes, the photon model used and whether the RG-A
fallback was taken.

The whole chain has been run end-to-end on real RG-D data once, as a
**diagnostic first pass**: Stage A skimmed the SIDIS train under a 2M-event/file
cap (95 slims, ≈5.5M DIS events, ≈2.5M $\pi^0$ candidates), `make_grid` froze the
binning from those slims, and Stage B/C produced $R_A$ and
$\Delta\langle p_T^2\rangle$. The machinery is validated — the $p_T$ broadening
increases monotonically with $A$, as expected — but the numbers are
**diagnostic**, and the truncation and the RG-A photon fallback are stamped into
every one of them.

Known blockers, all recorded in the provenance or the config rather than in
someone's memory:

- Results so far are the **diagnostic run above, not full luminosity** — Stage A
  read a 2M-event prefix of each file. A full run reuses the same machinery with
  the cap removed.
- **No GBT photon model exists for RG-D** (the map stops at run 16772). The RG-A
  fallback refuses by default, must be opted into explicitly, and stamps
  `gbt.fallback_used` into every output — which the Python stage then refuses to
  publish.
- The **beam polarisation** is not set; the BSA cannot be quoted without it.
- **No systematic uncertainty is evaluated.**

The binning grids in `config/binning/` are **no longer placeholders**: they are
the frozen equal-statistics edges `make_grid` computed from the diagnostic slims,
committed and hashed. A full-luminosity re-fit would move the edge positions but
not the schema.
