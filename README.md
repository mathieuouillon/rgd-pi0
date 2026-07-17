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
| `src/core/`           | Pure analysis code → static lib `pi0_core`. No HIPO. No ROOT beyond `ROOT::VecOps::RVec`. Unit-testable without a data file. |
| `src/tools/`          | ROOT-dependent programs (`.cxx`).                                          |
| `tests/`              | Catch2 v3 unit tests for `src/core`.                                       |
| `external/hipo-cpp/`  | **Vendored dependency — do not modify.** See below.                        |
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

### Running the built programs

`meson devenv -C build` exports `ROOT_INCLUDE_PATH` for you, which matters if an
interactive ROOT session needs to load the HipoDataFrame dictionary:

```sh
meson devenv -C build
```

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

Known blockers, all recorded in the provenance or the config rather than in
someone's memory:

- `config/binning/*.json` are **placeholders**. Real grids come from `make_grid`
  run over real RG-D slims.
- **No GBT photon model exists for RG-D** (the map stops at run 16772). The
  RG-A fallback refuses by default and must be opted into explicitly.
- The **beam polarisation** is not set; the BSA cannot be quoted without it.
- **No systematic uncertainty is evaluated.**
