"""Summary figures for the diagnostic first run of the rewritten pipeline.

Reads the Stage C result CSVs (``pi0.ratio``, ``pi0.broadening``) and draws the
three figures the note's results section shows:

* ``results_RA_vs_z.pdf``   -- R_A vs the sideband-subtracted ``<z>`` per target,
* ``results_dpt2_vs_A.pdf`` -- the A-dependence of the pT broadening,
* ``results_dpt2_vs_z.pdf`` -- the broadening vs ``<z>`` per target.

Each is stamped DIAGNOSTIC: the inputs carry the fallback + truncation
blockers, so nothing here is a measurement. Only matplotlib is used (no
mplhep), so this runs anywhere the extraction stack does.

    python -m pi0.plots_results --results <dir> --outdir note/figures

The R_A points are plotted at each z bin's true ``<z>`` (the ``z_mean`` column,
sideband-subtracted), NOT at the box midpoint -- the top z bin [0.37, 1.0]
therefore sits near its real centre, which is the whole point of the rewrite.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

#: (file/provenance target, mass number A, plot label). C, Cu, Sn -- ratioed to LD2.
TARGETS = [("CxC", 12.0, "C"), ("Cu", 63.5, "Cu"), ("Sn", 118.7, "Sn")]
N_Z = 5
COLORS = {"CxC": "#1f77b4", "Cu": "#d62728", "Sn": "#2ca02c"}


def _tofloat(x: str) -> float:
    try:
        return float(x)
    except ValueError:
        return np.nan


def _load(path: Path) -> dict[str, np.ndarray]:
    """A result CSV as name -> array, skipping the ``#`` provenance header."""
    with path.open() as fh:
        rows = [r for r in csv.reader(fh) if r and not r[0].lstrip().startswith("#")]
    hdr = rows[0]
    return {name: np.array([_tofloat(r[i]) for r in rows[1:]]) for i, name in enumerate(hdr)}


def _wmean(val: np.ndarray, err: np.ndarray) -> tuple[float, float]:
    """Inverse-variance weighted mean and its error over finite, err>0 entries."""
    m = np.isfinite(val) & np.isfinite(err) & (err > 0)
    if not m.any():
        return np.nan, np.nan
    w = 1.0 / err[m] ** 2
    return float(np.sum(w * val[m]) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w)))


def _by_z(d: dict[str, np.ndarray], val: str, err: str, zc: np.ndarray) -> np.ndarray:
    """Per z-bin (<z>, weighted mean, error), dropping empty bins."""
    out = []
    for k in range(N_Z):
        m = d["i_z"] == k
        mean, e = _wmean(d[val][m], d[err][m])
        if np.isfinite(mean) and np.isfinite(zc[k]):
            out.append((zc[k], mean, e))
    return np.array(out) if out else np.empty((0, 3))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Diagnostic result figures (R_A, Delta<pT2>).")
    p.add_argument("--results", type=Path, required=True,
                   help="directory holding ratio_<T>.csv and broadening_<T>.csv")
    p.add_argument("--outdir", type=Path, required=True, help="output directory for the PDFs")
    args = p.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({
        "savefig.dpi": 120, "savefig.bbox": "tight",
        "font.size": 10, "axes.grid": True, "grid.alpha": 0.25,
    })
    import matplotlib.pyplot as plt

    ratio = {t: _load(args.results / f"ratio_{t}.csv") for t, _, _ in TARGETS}
    broad = {t: _load(args.results / f"broadening_{t}.csv") for t, _, _ in TARGETS}

    # Representative <z> per z bin, from the ratio table's sideband-subtracted z_mean.
    with np.errstate(invalid="ignore"):
        zc = {t: np.array([np.nanmean(ratio[t]["z_mean"][ratio[t]["i_z"] == k]) for k in range(N_Z)])
              for t, _, _ in TARGETS}

    def stamp(fig: "plt.Figure") -> None:
        fig.text(0.5, 0.5, "DIAGNOSTIC — DO NOT QUOTE", ha="center", va="center",
                 fontsize=20, color="0.85", rotation=28, zorder=-10, alpha=0.7)

    # --- Fig 1: R_A vs <z> --------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    for t, _, lbl in TARGETS:
        xy = _by_z(ratio[t], "R", "sR", zc[t])
        ax.errorbar(xy[:, 0], xy[:, 1], yerr=xy[:, 2], marker="o", ms=4, capsize=2,
                    lw=1.2, color=COLORS[t], label=f"$R_{{{lbl}}}$")
    ax.axhline(1.0, color="0.4", lw=0.8, ls="--")
    ax.set_xlabel(r"$\langle z \rangle$  (sideband-subtracted)")
    ax.set_ylabel(r"$R_A$")
    ax.set_title(r"Multiplicity ratio vs $z$")
    ax.legend(frameon=False)
    stamp(fig)
    fig.savefig(args.outdir / "results_RA_vs_z.pdf")
    plt.close(fig)

    # --- Fig 2: Delta<pT2> vs A --------------------------------------------
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    A_arr, d_arr = [], []
    for t, A, lbl in TARGETS:
        mean, err = _wmean(broad[t]["delta"], broad[t]["sdelta"])
        A_arr.append(A); d_arr.append(mean)
        ax.errorbar(A, mean, yerr=err, marker="s", ms=7, capsize=3, color=COLORS[t], label=lbl)
    A_arr, d_arr = np.array(A_arr), np.array(d_arr)
    good = np.isfinite(d_arr) & (d_arr > 0)
    if good.sum() >= 2:
        n, lc = np.polyfit(np.log(A_arr[good]), np.log(d_arr[good]), 1)
        xs = np.linspace(A_arr.min(), A_arr.max(), 50)
        ax.plot(xs, np.exp(lc) * xs ** n, color="0.4", lw=0.9, label=fr"$\propto A^{{{n:.2f}}}$")
    ax.set_xlabel("$A$")
    ax.set_ylabel(r"$\Delta\langle p_T^2 \rangle$  [GeV$^2$], iv-weighted")
    ax.set_title(r"$p_T$ broadening vs $A$")
    ax.legend(frameon=False)
    stamp(fig)
    fig.savefig(args.outdir / "results_dpt2_vs_A.pdf")
    plt.close(fig)

    # --- Fig 3: Delta<pT2> vs <z> ------------------------------------------
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    for t, _, lbl in TARGETS:
        xy = _by_z(broad[t], "delta", "sdelta", zc[t])
        ax.errorbar(xy[:, 0], xy[:, 1], yerr=xy[:, 2], marker="o", ms=4, capsize=2,
                    lw=1.2, color=COLORS[t], label=lbl)
    ax.axhline(0.0, color="0.4", lw=0.8, ls="--")
    ax.set_xlabel(r"$\langle z \rangle$")
    ax.set_ylabel(r"$\Delta\langle p_T^2 \rangle$  [GeV$^2$]")
    ax.set_title(r"$p_T$ broadening vs $z$")
    ax.legend(frameon=False)
    stamp(fig)
    fig.savefig(args.outdir / "results_dpt2_vs_z.pdf")
    plt.close(fig)

    print(f"wrote 3 figures to {args.outdir}: results_RA_vs_z.pdf, "
          f"results_dpt2_vs_A.pdf, results_dpt2_vs_z.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
