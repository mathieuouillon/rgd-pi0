"""Diagnostics for the extraction stage.

Four plots, each answering a question the superseded analysis could not:

* :func:`plot_spectrum` -- same / scaled-mixed / subtracted with the fit
  overlaid. Does the subtraction work at all?
* :func:`plot_fit_vs_bin` -- fitted ``mu``, ``sigma`` and ``chi2/ndf`` per bin.
  Is the peak where it should be, everywhere?
* :func:`plot_alpha_vs_bin` -- the sideband scale per bin. Is the background
  normalisation stable across kinematics?
* :func:`plot_abscissa_vs_centre` -- **the money plot**. How far is each bin's
  count-weighted, sideband-subtracted abscissa from the geometric box centre
  the old analysis reported it at? This is what the rewrite bought.

Plus :func:`sideband_slope_diagnostic`, which measures the one assumption the
abscissa estimator rests on.

Every figure is stamped with the provenance of the file it came from.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from . import extract, io
from .config import Config, load_config
from .extract import BinYield, FitResult, gaussian, subtract

__all__ = [
    "plot_spectrum",
    "plot_fit_vs_bin",
    "plot_alpha_vs_bin",
    "plot_abscissa_vs_centre",
    "sideband_slope_diagnostic",
    "main",
]

_AXIS_LABEL = {
    "q2": r"$Q^2$ [GeV$^2$]",
    "xb": r"$x_B$",
    "z": r"$z$",
    "pt2": r"$p_T^2$ [GeV$^2$]",
}


def _stamp(fig, data: io.StageBData) -> None:
    """Print the file's provenance along the bottom of a figure."""
    txt = data.stamp()
    colour = "crimson" if data.blockers else "0.35"
    fig.text(0.005, 0.004, txt, fontsize=5.0, color=colour, ha="left", va="bottom")
    if any(x.fatal for x in data.blockers):
        fig.text(
            0.5, 0.5, "DIAGNOSTIC ONLY\nDO NOT QUOTE", fontsize=44, color="crimson",
            alpha=0.10, ha="center", va="center", rotation=30, zorder=100,
            transform=fig.transFigure,
        )


def plot_spectrum(
    n_same: NDArray,
    n_mixed: NDArray,
    cfg: Config,
    by: BinYield,
    path: Path,
    data: io.StageBData,
    title: str = "",
) -> None:
    """The mgg spectrum: same, scaled mixed, subtracted, with the fit overlaid.

    The shaded band is the +-``n_sigma`` window whose **sum** is the yield. The
    dashed curve is the fitted Gaussian, drawn to show that it located the
    window -- its integral is *not* the yield and is not drawn as one.

    Parameters
    ----------
    n_same, n_mixed : numpy.ndarray
        The spectra to draw.
    cfg : Config
        Configuration.
    by : BinYield
        The extraction whose alpha, fit and window are overlaid.
    path : pathlib.Path
        Output file.
    data : pi0.io.StageBData
        Only for the provenance stamp.
    title : str, optional
        Extra title text.
    """
    m, cuts = cfg.mgg.centres, cfg.extraction
    alpha = by.alpha
    s, var = subtract(n_same, n_mixed, alpha)

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True,
                                   gridspec_kw={"height_ratios": [1.15, 1]})
    ax0.step(m, n_same, where="mid", color="k", lw=1.0, label="same-event")
    ax0.step(m, alpha * n_mixed, where="mid", color="tab:blue", lw=1.0,
             label=rf"mixed $\times\ \alpha$ = {alpha:.4g}")
    ax0.axvspan(cuts.sideband_min, cuts.sideband_max, color="tab:orange", alpha=0.13,
                label="sideband (sets $\\alpha$)")
    ax0.set_ylabel("counts / 1.5 MeV")
    ax0.legend(fontsize=8, loc="upper right")
    ax0.set_title(title or "$m_{\\gamma\\gamma}$ spectrum", fontsize=10)

    ax1.errorbar(m, s, yerr=np.sqrt(var), fmt=".", ms=3, lw=0.7, color="k",
                 label="subtracted $S(m)$")
    ax1.axhline(0, color="0.6", lw=0.7)
    if by.fit is not None and by.fit.converged:
        f = by.fit
        fine = np.linspace(m[0], m[-1], 800)
        ax1.plot(fine, gaussian(fine, f.amp, f.mu, f.sigma), "--", color="tab:red", lw=1.2,
                 label=(rf"fit: $\mu$={f.mu*1000:.1f} MeV, $\sigma$={f.sigma*1000:.1f} MeV, "
                        rf"$\chi^2/$ndf={f.chi2_ndf:.2f}"))
        lo, hi = f.mu - cuts.n_sigma_window * f.sigma, f.mu + cuts.n_sigma_window * f.sigma
        ax1.axvspan(lo, hi, color="tab:green", alpha=0.13,
                    label=rf"$\pm{cuts.n_sigma_window:g}\sigma$ window (its SUM is the yield)")
        ax1.axvline(f.range_lo, color="tab:red", ls=":", lw=0.8)
        ax1.axvline(f.range_hi, color="tab:red", ls=":", lw=0.8,
                    label="final (core-restricted) fit range")
        ax1.annotate(
            f"$Y$ = {by.value:.1f} $\\pm$ {by.error:.1f}  (window sum, NOT the Gaussian integral)\n"
            f"$S/(S+B)$ = {by.signal_fraction:.3f}",
            xy=(0.02, 0.96), xycoords="axes fraction", va="top", fontsize=8,
        )
    ax1.set_xlabel(r"$m_{\gamma\gamma}$ [GeV]")
    ax1.set_ylabel("counts / 1.5 MeV")
    ax1.legend(fontsize=7.5, loc="upper right")
    _stamp(fig, data)
    fig.tight_layout(rect=(0, 0.015, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_fit_vs_bin(res: list[BinYield], path: Path, data: io.StageBData) -> None:
    """Fitted ``mu``, ``sigma`` and ``chi2/ndf`` against 4D bin index."""
    ok = [r for r in res if r.fit is not None and r.fit.converged]
    fig, axes = plt.subplots(3, 1, figsize=(7.6, 6.6), sharex=True)
    if not ok:
        for ax in axes:
            ax.text(0.5, 0.5, "no converged fits", ha="center", va="center",
                    transform=ax.transAxes, color="crimson")
    else:
        x = [r.bin4d for r in ok]
        for ax, vals, lab in (
            (axes[0], [r.fit.mu * 1000 for r in ok], r"fitted $\mu$ [MeV]"),
            (axes[1], [r.fit.sigma * 1000 for r in ok], r"fitted $\sigma$ [MeV]"),
            (axes[2], [r.fit.chi2_ndf for r in ok], r"$\chi^2/$ndf"),
        ):
            ax.plot(x, vals, ".", ms=3)
            ax.set_ylabel(lab, fontsize=9)
            ax.grid(alpha=0.25)
        axes[0].axhline(134.98, color="tab:green", lw=0.8, ls="--", label=r"PDG $m_{\pi^0}$")
        axes[0].legend(fontsize=8)
        axes[2].set_yscale("log")
    axes[2].set_xlabel("4D bin index")
    axes[0].set_title("Mass fit vs bin", fontsize=10)
    _stamp(fig, data)
    fig.tight_layout(rect=(0, 0.015, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_alpha_vs_bin(res: list[BinYield], path: Path, data: io.StageBData) -> None:
    """The sideband scale ``alpha`` against 4D bin index.

    A stable ``alpha`` means the mixed sample tracks the same-event background
    normalisation consistently across kinematics; structure here is the first
    place a subtraction problem shows up.
    """
    pts = [(r.bin4d, r.alpha) for r in res if np.isfinite(r.alpha)]
    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    if pts:
        x, y = zip(*pts)
        ax.plot(x, y, ".", ms=3)
        ax.axhline(float(np.median(y)), color="tab:red", lw=0.8, ls="--",
                   label=f"median = {np.median(y):.4g}")
        ax.set_yscale("log")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no bin has a defined alpha", ha="center", va="center",
                transform=ax.transAxes, color="crimson")
    ax.set_xlabel("4D bin index")
    ax.set_ylabel(r"$\alpha$ = $\Sigma$same$_{SB}$ / $\Sigma$mixed$_{SB}$")
    ax.set_title("Sideband scale vs bin", fontsize=10)
    ax.grid(alpha=0.25)
    _stamp(fig, data)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_abscissa_vs_centre(
    res: list[BinYield], cfg: Config, path: Path, data: io.StageBData, note: str = ""
) -> int:
    """**The money plot**: count-weighted abscissa against geometric box centre.

    Each point is one bin: the x axis is where the old analysis would have
    reported it (the box midpoint), the y axis is where the pi0 in it actually
    are. The diagonal is where the two agree. The horizontal bars are the bins'
    full widths -- the outermost bin in each dimension is enormously wide, and
    that is exactly where the two disagree most.

    Returns
    -------
    int
        How many bins are drawn.
    """
    centres = cfg.binning.centres_4d()
    widths = cfg.binning.widths_4d()
    ok = [r for r in res if r.abscissae and all(np.isfinite(v) for v in r.abscissae.values())]

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 8.0))
    for ax, name in zip(axes.ravel(), ("q2", "xb", "z", "pt2")):
        if ok:
            cx = np.array([centres[name][r.bin4d] for r in ok])
            wx = np.array([widths[name][r.bin4d] for r in ok])
            cy = np.array([r.abscissae[name] for r in ok])
            ax.errorbar(cx, cy, xerr=wx / 2, fmt="o", ms=3.4, lw=0.7, alpha=0.8,
                        color="tab:blue", ecolor="0.7")
            span = [min(cx.min() - wx.max() / 2, cy.min()), max(cx.max() + wx.max() / 2, cy.max())]
            pad = 0.05 * (span[1] - span[0] + 1e-9)
            ax.plot([span[0] - pad, span[1] + pad], [span[0] - pad, span[1] + pad], "--",
                    color="0.4", lw=0.9, label="$\\langle X\\rangle$ = box centre")
            med = float(np.median(cy - cx))
            ax.set_title(f"{_AXIS_LABEL[name]}   median shift = {med:+.3g}", fontsize=9)
            ax.legend(fontsize=7.5, loc="upper left")
        else:
            ax.text(0.5, 0.5, "no bin has a computable abscissa", ha="center", va="center",
                    transform=ax.transAxes, color="crimson", fontsize=9)
            ax.set_title(_AXIS_LABEL[name], fontsize=9)
        ax.set_xlabel(f"geometric box centre  (what the old analysis reported)", fontsize=8)
        ax.set_ylabel(f"count-weighted, sideband-subtracted mean", fontsize=8)
        ax.grid(alpha=0.25)

    fig.suptitle(
        "Where the bin actually is, vs where it was reported\n"
        + (note or "x error bars are the bins' FULL widths"),
        fontsize=10,
    )
    _stamp(fig, data)
    fig.tight_layout(rect=(0, 0.015, 1, 0.95))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return len(ok)


def sideband_slope_diagnostic(
    n_same: NDArray, sum_x: NDArray, cfg: Config, name: str
) -> tuple[float, float, float]:
    """Measure the assumption the abscissa estimator rests on.

    The estimator assumes the *background's* ``<X>`` is constant in ``m_gg``
    between the peak and the sideband, so that a single number measured in the
    sideband can be subtracted from under the peak
    (:func:`pi0.extract.abscissa`). That assumption is testable *inside* the
    sideband, where the sample is essentially pure background: split it and
    look for a trend.

    Parameters
    ----------
    n_same, sum_x : numpy.ndarray
        The same-event counts and ``sum_X`` per mgg bin.
    cfg : Config
        Configuration, for the sideband range.
    name : str
        Axis name, for the message.

    Returns
    -------
    slope : float
        d``<X>_bkg``/d``m_gg``, per GeV, from a counts-weighted linear fit
        across the sideband bins.
    slope_err : float
        Its uncertainty.
    extrapolated_shift : float
        ``slope * (m_peak - m_sideband_centre)`` -- roughly how far a flat
        extrapolation misplaces the background's ``<X>`` at the peak. This is
        an unpropagated systematic on every abscissa; it is reported, not
        corrected.

    Notes
    -----
    A slope consistent with zero means the assumption holds and the abscissa is
    trustworthy. A significant slope means it does not, and the honest fix is
    the mixed-event ``sum_X`` columns described in :func:`pi0.extract.abscissa`.
    """
    m = cfg.mgg.centres
    sb = cfg.mgg.mask(cfg.extraction.sideband_min, cfg.extraction.sideband_max)
    sel = sb & (n_same > 0)
    if int(sel.sum()) < 3:
        return float("nan"), float("nan"), float("nan")
    x = m[sel]
    y = sum_x[sel] / n_same[sel]
    # Weight by counts: <X> in a bin with n entries has variance ~ 1/n.
    w = n_same[sel].astype(float)
    coef, cov = np.polyfit(x, y, 1, w=np.sqrt(w), cov=True)
    slope, slope_err = float(coef[0]), float(np.sqrt(cov[0, 0]))
    m_peak = cfg.extraction.mu_seed
    m_sb = float(np.average(x, weights=w))
    return slope, slope_err, slope * (m_peak - m_sb)


def main(argv: list[str] | None = None) -> int:
    """CLI: run the QA suite over one Stage B file."""
    p = argparse.ArgumentParser(description="QA diagnostics for the pi0 extraction stage.")
    p.add_argument("--file", required=True, type=Path, help="Stage B file")
    p.add_argument("--config", required=True, type=Path, help="config/ directory")
    p.add_argument("--outdir", required=True, type=Path, help="directory for the figures")
    p.add_argument("--allow-unpublishable", action="store_true",
                   help="proceed despite fatal provenance blockers. Figures get a watermark.")
    p.add_argument(
        "--shared-window", action="store_true",
        help="locate every bin's mass window with the fit to the SUMMED spectrum instead of the "
             "bin's own. NOT the production path -- a fallback for files whose per-bin statistics "
             "cannot support a fit. Output is labelled accordingly.",
    )
    args = p.parse_args(argv)

    cfg = load_config(args.config)
    data = io.load(args.file, cfg, allow_unpublishable=args.allow_unpublishable)
    args.outdir.mkdir(parents=True, exist_ok=True)

    # --- the summed spectrum: always meaningful, even when no single bin is.
    n_same, n_mixed, sums = extract.summed_spectrum(data, cfg)
    total = extract.extract_bin(-1, n_same, n_mixed, sums, cfg)
    print(f"\nSUMMED SPECTRUM over all {cfg.binning.n_4d} bins")
    print(f"  same-event entries : {int(n_same.sum())}")
    print(f"  mixed-event pairs  : {int(n_mixed.sum())}")
    print(f"  alpha              : {total.alpha:.6g}")
    if total.fit is not None and total.fit.converged:
        f = total.fit
        print(f"  fit                : mu = {f.mu*1000:.2f} +- {f.mu_err*1000:.2f} MeV, "
              f"sigma = {f.sigma*1000:.2f} +- {f.sigma_err*1000:.2f} MeV, "
              f"amp = {f.amp:.3g}, chi2/ndf = {f.chi2:.1f}/{f.ndf} = {f.chi2_ndf:.3f}, "
              f"core iters = {f.n_iter}, final range = [{f.range_lo:.4f}, {f.range_hi:.4f}]")
        print(f"  YIELD (window sum) : {total.value:.2f} +- {total.error:.2f}")
        print(f"  Gaussian integral  : {extract.gaussian_integral_yield(f, cfg.mgg):.2f}  "
              f"<- NOT the yield; shown only to make the difference visible")
        print(f"  S/(S+B) in window  : {total.signal_fraction:.4f}")
        print(f"  abscissae          : " + ", ".join(
            f"<{k}> = {v:.4g}" for k, v in total.abscissae.items()))
    else:
        print(f"  fit                : FAILED ({total.reason})")

    print("\nABSCISSA ASSUMPTION -- background <X> vs m_gg inside the sideband:")
    for name, sx in sums.items():
        sl, sle, shift = sideband_slope_diagnostic(n_same, sx, cfg, name)
        flag = ""
        if np.isfinite(sl) and np.isfinite(sle) and sle > 0 and abs(sl) / sle > 3:
            flag = "  <-- SIGNIFICANT: the flat-background assumption is not supported here"
        print(f"  d<{name}>/dm_gg = {sl:+.4g} +- {sle:.4g} /GeV  =>  extrapolated shift at the "
              f"peak {shift:+.4g}{flag}")

    plot_spectrum(n_same, n_mixed, cfg, total, args.outdir / "qa_spectrum_summed.png", data,
                  title="Summed over all 4D bins (diagnostic: mixes kinematics)")

    # --- per-bin
    shared = total.fit if (args.shared_window and total.fit and total.fit.converged) else None
    res = extract.extract_all(data, cfg, shared_fit=shared)
    n_ok = sum(1 for r in res if r.ok)
    counts: dict[str, int] = {}
    for r in res:
        counts[r.reason] = counts.get(r.reason, 0) + 1
    print(f"\nPER-BIN EXTRACTION: {n_ok} of {cfg.binning.n_4d} bins passed every gate.")
    for code, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"    {n:6d}  {code}")

    plot_fit_vs_bin(res, args.outdir / "qa_fit_vs_bin.png", data)
    plot_alpha_vs_bin(res, args.outdir / "qa_alpha_vs_bin.png", data)
    note = (
        "WINDOW FROM THE SUMMED-SPECTRUM FIT (--shared-window), not per-bin fits: DIAGNOSTIC"
        if shared else "x error bars are the bins' FULL widths"
    )
    n_drawn = plot_abscissa_vs_centre(res, cfg, args.outdir / "qa_abscissa_vs_centre.png", data,
                                      note)
    print(f"\nmoney plot: {n_drawn} bins have a computable abscissa.")
    print(f"figures written to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
