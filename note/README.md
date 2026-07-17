# RG-D π⁰ Analysis Note (Typst)

Physics analysis note for the CLAS12 RG-D neutral-pion analysis: multiplicity
ratios, transverse-momentum broadening and beam-spin asymmetries on
LD2 / CxC / Cu / Sn.

## Quick start

```bash
typst compile main.typ      # build main.pdf
typst watch main.typ        # rebuild on save
```

Requires Typst ≥ 0.15. Packages (`subpar`, `glossarium`, `showybox`,
`physica`) download automatically on first compile. Fonts: STIX Two
Text/Math, JetBrains Mono.

## Layout

```
main.typ                       # Title metadata + section manifest
template/                      # Reused from RGDCommonAnalysisNote
  lib.typ                      #   Public surface: re-exports + #note() show-rule
  theme.typ                    #   Fonts, colors, page geometry, headings
  figures.typ                  #   wide-figure / subfig2 / subfig3 / subfig2x2
  callouts.typ                 #   note-box / important-box / result-box / warning-box
  glossary.typ                 #   Acronyms (extended with pi0-specific entries)
  refs.typ, tables.typ
sections/
  01_introduction.typ          # Physics motivation, observables
  02_dataset.typ               # Beam, targets, runs, data flow, code provenance
  03_event_selection.typ       # Electron ID, DIS cuts, photon ID
  04_pi0_reconstruction.typ    # Invariant mass, pairing, SIDIS kinematics
  05_binning.typ               # kd-tree adaptive binning
  06_background.typ            # Event mixing, sideband subtraction
  07_multiplicity_ratio.typ    # R_A extraction
  08_pt_broadening.typ         # Delta<pT^2>
  09_bsa.typ                   # A_LU, azimuthal moments
  10_systematics.typ           # Corrections and systematics (i.e. their absence)
  11_results.typ               # Current results, with caveats
  12_summary.typ               # Summary, blockers, next steps
  99_appendix.typ              # SF params, banks, dead config, MLM, schemas
figures/                       # Figure PDFs, copied from the farm plots/ tree
```

## Sources

| What | Where |
|---|---|
| C++ framework (kd-tree era — **produced the results**) | `~/Documents/tmp/clas-framework` @ `e8334b1` |
| C++ framework (current — kd-tree **removed**) | `~/Documents/tmp/clas-analysis-1` @ `29150c0` |
| Python extraction | `~/rg-c_farm/analyses/pi0/` |
| Result CSVs | `~/rg-c_farm/results/pi0_*/` |
| Plots | `~/rg-c_farm/plots/pi0_*/` |
| Production ROOT files | `/work/clas12b/users/ouillon/clas-framework/pi0_multiplicity_{LD2,CxC,Cu,Sn}_OB.root` (farm-only) |

## Important

This note documents the **kd-tree binning** (`clas-framework` @ `e8334b1`),
because that is what produced every number and figure in it. The newer
`clas-analysis-1` checkout **deleted** that code on 14 July 2026 (commit
`eac29a1`) and replaced it with factorized (Q²,x_B)×(z,pT²) grids. See
§2.4 (Code provenance) and §5.6 (What the refactor changes).

All quoted numbers are **statistical only** and provisional. See §10
(Corrections and systematic uncertainties) — no systematic has been
evaluated — and §12.2 for the list of publication blockers.

Numbers in the note were recomputed directly from the result CSVs rather
than transcribed.
