# `rgd-pi0` — the Python extraction stage

Reads **Stage B**'s binned ROOT file and produces the four observables. uproot
only: no ROOT, no dictionaries, no `TTree::Draw`.

```
Stage A (skim)  ->  Stage B (bin)  ->  binned.root  ->  THIS PACKAGE  ->  CSV + figures
```

## Running it

Nothing needs installing:

```sh
UV="uv run --with uproot --with numpy --with scipy --with matplotlib"
cd /Users/mathieuouillon/Documents/C++/RGD_pi0_analysis/python

# QA: the summed spectrum, the fit, alpha vs bin, and the abscissa-vs-centre plot
$UV python -m pi0.qa       --file binned.root --config ../config --outdir figures/

# R_A: needs TWO files -- a nuclear target and LD2. R_A is a ratio; one file cannot make it.
$UV python -m pi0.ratio    --target CxC.root --ld2 LD2.root --config ../config --out ra.csv

# Delta<pT2>, sideband-subtracted
$UV python -m pi0.broadening --target CxC.root --ld2 LD2.root --config ../config --out ptb.csv

# A_LU: REFUSES to run without an explicit beam polarization
$UV python -m pi0.bsa --file LD2.root --config ../config --out bsa.csv \
      --polarization 0.86 --polarization-err 0.03

# tests
$UV --with pytest python -m pytest tests/ -q
```

Two flags exist to say "yes, I know":

* `--allow-unpublishable` — proceed despite fatal provenance blockers. Every
  figure gets a `DIAGNOSTIC ONLY` watermark and every CSV a `DO NOT QUOTE`
  header. Without it, a file whose photons came from a fallback GBT model, or
  whose grids are placeholders, **raises**.
* `--shared-window` (`pi0.qa` only) — locate every bin's mass window with the
  fit to the *summed* spectrum instead of the bin's own. Not the production
  path; a fallback for files too thin to fit per bin. The output says so.

## The modules

| Module | What it owns |
|---|---|
| `config.py` | `cuts.json` + the two grid JSONs. **The only place that knows the index formulas.** A missing key raises; nothing is defaulted. |
| `io.py` | The four trees + both provenance blocks. **Enforces** the provenance gate. |
| `extract.py` | Sideband scale, subtraction, the core-restricted Gaussian fit, the ±3σ window **sum**, the quality gates, and the count-weighted sideband-subtracted **abscissa**. |
| `ratio.py` | `R_A` per 4D bin → a tidy CSV / numpy structured array. |
| `broadening.py` | `Δ⟨pT²⟩` per 3D bin, **sideband-subtracted**. |
| `bsa.py` | `A_LU`, dilution-corrected, with a mandatory explicit polarization and a producer-side fit-quality mask. |
| `qa.py` | Diagnostics, including the abscissa-vs-box-centre money plot. |

## Three invariants

**1. No cut value is hard-coded.** Every number comes from `config/cuts.json`
and the two grid JSONs. `config.require()` raises `ConfigError` on a missing
key rather than defaulting — a silently defaulted cut is the failure mode this
project exists to end. Writing this stage needed four keys `cuts.json` did not
yet have, and they were **added to the config**, not to the code:
`bsa/fit_bounds`, `bsa/quality`, `bsa/polarization` (null on purpose) and
`broadening/min_counts`.

**2. No result is ever reported at a geometric bin centre.** Every abscissa is
the count-weighted, sideband-subtracted mean. A bin whose abscissa cannot be
computed is **dropped and counted**, never back-filled with a box midpoint.
`Grid1D.centres` exists solely so the QA plot can show how far off the old
positions were.

**3. Provenance is enforced, not decorated.** `io.load` reads both blocks and
refuses to proceed when they say the data is not fit to quote.

## The yield, precisely

1. `alpha = Σ n_same(SB) / Σ n_mixed(SB)`, single **high**-mass sideband.
2. `S(m) = n_same(m) − alpha·n_mixed(m)`, `var(m) = n_same(m) + alpha²·n_mixed(m)`.
3. Pure Gaussian fit (no background term — it is already gone) over the fit range.
4. Iterative core restriction to `[µ ± 1.5σ]`, ≤2 iterations, stop at <1% change.
5. **`Y = Σ_{|m−µ|<3σ} S(m)`**, `σ_Y = sqrt(Σ var)`.

Step 5 is not a detail. The superseded script's own docstring claimed
`N = amp·σ·√(2π)/binW` — **that is not what produced its results**; the
analytic integral existed only in a demonstration plot. `gaussian_integral_yield`
is implemented here for exactly one purpose: so the test suite can inject a
non-Gaussian shoulder and prove the two differ by ~25%, pinning which one the
code uses.

## The abscissa, and its one assumption

Stage B's `sum_*` columns are **same-event only** — no mixed pair contributes,
because a mixed pair is not a π⁰. So the background's `⟨X⟩` is estimated from
the sideband and its share removed:

```
⟨X⟩_B  ≈ sum_X(SB) / n_same(SB)          (the sideband is ~pure background)
B_W     = alpha · Σ n_mixed(W)
⟨X⟩_S  ≈ (sum_X(W) − ⟨X⟩_B · B_W) / (n_same(W) − B_W)
```

This is algebraically exact **given** that `⟨X⟩_B` does not vary with `m_gg`
between the peak and the sideband. That assumption is the weak point, and it is
not obviously true for `z`: `m_gg` correlates with the pair's energy, and
`z = E_π⁰/ν`, so the background's `⟨z⟩` plausibly rises with `m_gg` and a flat
extrapolation down to the peak would over-subtract.

`qa.sideband_slope_diagnostic` measures that trend inside the sideband, where
the sample is essentially pure background, and prints it. It is reported, not
corrected — no error carries it.

**The exact fix belongs in Stage B**: write a mixed-event `sum_X` per
`(bin4d, imgg)` — `sum_q2_mixed`, `sum_xb_mixed`, `sum_z_mixed`,
`sum_pt2_mixed` — each mixed pair's kinematics computed against its *current*
event's DIS electron. Then

```
⟨X⟩_S = (sum_X_same(W) − alpha·sum_X_mixed(W)) / (n_same(W) − alpha·n_mixed(W))
```

with no flatness assumption at all — exactly the structure already used for the
counts. It is a small change and it retires the last approximation in the
abscissa.

## Known omissions, carried deliberately

Stated here because they are invisible in the numbers:

* The uncertainty on `alpha` is not propagated; the subtraction treats it as exact.
* The ±3σ window sum's error ignores the uncertainty on the fitted `µ` and `σ`
  that *define* the window.
* **`R_CxC`, `R_Cu`, `R_Sn` are correlated** — they share the LD2 yield and DIS
  count in every denominator — but `σ_R` treats them as independent. Harmless
  for any single `R_A`; **wrong for any fit of the A-dependence**, which will
  understate its own uncertainty. The shared LD2 columns are kept in every CSV
  row so a caller can build the covariance.
* `σ²_⟨pT²⟩ = (⟨pT⁴⟩ − ⟨pT²⟩²)/N` ignores the variance the subtraction injects.
* In `bsa`, `S/(S+B)` and `A_fit` come from the same events and are correlated;
  that correlation is ignored.
