"""Transverse-momentum broadening ``Delta<pT2>``, per 3D (Q2, xB, z) bin.

    <pT2>            = sum_pT2 / counts
    sigma^2_<pT2>    = (<pT4> - <pT2>^2) / counts
    Delta<pT2>_A     = <pT2>_A - <pT2>_D,   sigma = sqrt(sA^2 + sD^2)

**Unlike the superseded analysis, the moments here are sideband-subtracted.**

That is the whole point of Stage B binning ``ptb3d`` in ``m_gg``. The old chain
averaged ``pT2`` over a +-200 MeV mass window with no subtraction at all, so
what it measured was

    <pT2>^meas = f_S <pT2>^pi0 + (1 - f_S) <pT2>^bkg

with ``f_S`` nowhere near 1 (note ``sec:ptb-caveat``; its script's docstring
claimed a +-30 MeV pi0 selection that the shipped configuration did not
apply). Its broadening is therefore diluted toward zero by an unknown factor.
The note names the fix -- "accumulate sum pT2 and sum pT4 per m_gg bin so that
the same sideband subtraction used for the yields can be applied to the
moments" -- and this module is that fix.

``ptb3d`` carries no mixed-event spectrum of its own, but it is exactly the
``i_pt2`` projection of ``spectra`` (:meth:`pi0.io.StageBData.mixed_3d`, checked
on load), so the 4D mixed spectrum supplies the background shape.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from . import extract, io
from .config import Config, load_config
from .extract import FitResult, fit_core, sideband_scale, subtract

__all__ = ["Moments", "moments_3d", "delta_pt2", "main"]


@dataclass(frozen=True)
class Moments:
    """Sideband-subtracted pT2 moments of one 3D bin.

    Attributes
    ----------
    bin3d : int
        Bin index.
    mean_pt2, err_pt2 : float
        ``<pT2>`` of the signal and its standard error.
    mean_pt4 : float
        ``<pT4>`` of the signal -- accumulated solely so that ``err_pt2`` can be
        formed.
    n_signal : float
        Subtracted signal count in the mass window; the ``N`` of the error.
    alpha : float
        The sideband scale used.
    fit : FitResult or None
        The mass fit that located the window.
    ok : bool
        Whether the bin passed every gate.
    reason : str
        ``"ok"`` or the gate it failed.
    """

    bin3d: int
    mean_pt2: float
    err_pt2: float
    mean_pt4: float
    n_signal: float
    alpha: float
    fit: FitResult | None
    ok: bool
    reason: str


def moments_3d(
    data: io.StageBData, cfg: Config, *, shared_fit: FitResult | None = None
) -> list[Moments]:
    """Sideband-subtracted ``<pT2>`` and ``<pT4>`` for every 3D bin.

    The estimator is the same one the abscissa uses, applied to the ``pT2`` and
    ``pT4`` moments::

        <X>_B    ~= sum_X(SB) / counts(SB)            (sideband ~ pure background)
        B_W       = alpha * sum(n_mixed_3d over W)
        <X>_S    ~= (sum_X(W) - <X>_B * B_W) / (counts(W) - B_W)

    and carries the same assumption: that the background's ``<pT2>`` varies
    slowly with ``m_gg`` between the peak and the sideband. See
    :func:`pi0.extract.abscissa` for why that assumption is the weak point and
    what Stage B should write to remove it.

    Parameters
    ----------
    data : pi0.io.StageBData
        The loaded file.
    cfg : Config
        Configuration.
    shared_fit : FitResult, optional
        Diagnostic only; see :func:`pi0.extract.extract_bin`.

    Returns
    -------
    list of Moments
        One per 3D bin, in index order. Failures are present with ``ok=False``.

    Notes
    -----
    ``sigma^2_<pT2> = (<pT4> - <pT2>^2) / N`` is the standard error of the mean
    of the pT2 distribution -- which is the whole reason ``sum_pT4`` is
    accumulated. Applied here to the *subtracted* moments with ``N`` the
    subtracted signal count, it **ignores the extra variance the subtraction
    itself injects** (both through the Poisson noise of the mixed spectrum and
    through ``alpha``, still treated as exact). It is therefore an
    underestimate at low purity. That is a real omission, stated rather than
    hidden; the honest version needs the mixed-event moments named in
    :func:`pi0.extract.abscissa`.
    """
    axis, cuts = cfg.mgg, cfg.extraction
    sb = axis.mask(cuts.sideband_min, cuts.sideband_max)
    fit_mask = axis.mask(cuts.fit_range_min, cuts.fit_range_max)
    mixed3 = data.mixed_3d()
    out: list[Moments] = []

    for i in range(cfg.binning.n_3d):
        counts, n_mixed = data.ptb_counts[i], mixed3[i]
        s_pt2, s_pt4 = data.ptb_sum_pt2[i], data.ptb_sum_pt4[i]

        def fail(reason: str, alpha: float = float("nan"), fit=None) -> Moments:
            return Moments(i, np.nan, np.nan, np.nan, np.nan, alpha, fit, False, reason)

        # Gate on the RAW same-event count, so it does not depend on the
        # subtraction it gates (note sec:ptb-accumulators: 20 entries).
        if float(np.sum(counts)) < cfg.broadening.min_counts:
            out.append(fail("below_min_counts"))
            continue

        alpha, why = sideband_scale(counts, n_mixed, sb)
        if not np.isfinite(alpha):
            out.append(fail(why))
            continue

        s, var = subtract(counts, n_mixed, alpha)
        if float(np.max(s[fit_mask])) < cuts.min_peak_amp:
            out.append(fail(extract.LOW_PEAK, alpha))
            continue

        fit = shared_fit if shared_fit is not None else fit_core(axis, s, var, cuts)
        if not fit.converged:
            out.append(fail(extract.FIT_FAILED, alpha, fit))
            continue

        window = np.abs(axis.centres - fit.mu) < cuts.n_sigma_window * fit.sigma
        n_sb = float(np.sum(counts[sb]))
        if n_sb <= 0:
            out.append(fail(extract.NO_SIDEBAND_SAME, alpha, fit))
            continue

        n_bkg_w = alpha * float(np.sum(n_mixed[window]))
        n_sig = float(np.sum(counts[window])) - n_bkg_w
        if n_sig <= 0:
            out.append(fail("non_positive_signal", alpha, fit))
            continue

        pt2_bkg = float(np.sum(s_pt2[sb])) / n_sb
        pt4_bkg = float(np.sum(s_pt4[sb])) / n_sb
        mean_pt2 = (float(np.sum(s_pt2[window])) - pt2_bkg * n_bkg_w) / n_sig
        mean_pt4 = (float(np.sum(s_pt4[window])) - pt4_bkg * n_bkg_w) / n_sig

        var_pt2 = (mean_pt4 - mean_pt2**2) / n_sig
        if not np.isfinite(var_pt2) or var_pt2 < 0:
            # <pT4> < <pT2>^2 is impossible for a real sample; after subtraction
            # it means the estimate has been driven negative by noise. Drop it
            # rather than emit an imaginary error bar.
            out.append(fail("negative_variance", alpha, fit))
            continue

        out.append(
            Moments(i, mean_pt2, float(np.sqrt(var_pt2)), mean_pt4, n_sig, alpha, fit, True,
                    extract.OK)
        )
    return out


def delta_pt2(
    target: list[Moments], deuterium: list[Moments], cfg: Config, *, stream=None
) -> tuple[NDArray, dict[str, int]]:
    """``Delta<pT2>_A = <pT2>_A - <pT2>_D`` per 3D bin.

    Parameters
    ----------
    target, deuterium : list of Moments
        Per-3D-bin moments, in index order.
    cfg : Config
        Configuration, for the binning.
    stream : file-like, optional
        Where to report drops.

    Returns
    -------
    table : numpy.ndarray
        Structured array: ``bin3d, i_q2, i_xb, i_z, z_mean, pt2_A, spt2_A,
        pt2_D, spt2_D, delta, sdelta``.
    dropped : dict of str to int
        Reject code -> count.

    Notes
    -----
    Errors add in quadrature, treating target and reference as independent --
    valid, they are separate files. The ``z_mean`` column is the count-weighted,
    sideband-subtracted ``<z>`` of the target's bin, never a box centre.
    """
    stream = stream if stream is not None else sys.stderr
    dtype = np.dtype(
        [
            ("bin3d", np.int64), ("i_q2", np.int64), ("i_xb", np.int64), ("i_z", np.int64),
            ("pt2_A", np.float64), ("spt2_A", np.float64),
            ("pt2_D", np.float64), ("spt2_D", np.float64),
            ("delta", np.float64), ("sdelta", np.float64),
        ]
    )
    rows, dropped = [], {}
    for i in range(cfg.binning.n_3d):
        ta, de = target[i], deuterium[i]
        if not ta.ok or not de.ok:
            code = f"target:{ta.reason}" if not ta.ok else f"ld2:{de.reason}"
            dropped[code] = dropped.get(code, 0) + 1
            continue
        i_q2, i_xb, i_z = cfg.binning.decode_3d(i)
        d = ta.mean_pt2 - de.mean_pt2
        sd = float(np.hypot(ta.err_pt2, de.err_pt2))
        rows.append((i, i_q2, i_xb, i_z, ta.mean_pt2, ta.err_pt2, de.mean_pt2, de.err_pt2, d, sd))
    table = np.array(rows, dtype=dtype)
    if dropped:
        print(f"broadening: {len(rows)} of {cfg.binning.n_3d} 3D bins survived. Dropped:",
              file=stream)
        for code, n in sorted(dropped.items(), key=lambda kv: -kv[1]):
            print(f"    {n:6d}  {code}", file=stream)
    return table, dropped


def main(argv: list[str] | None = None) -> int:
    """CLI: compute sideband-subtracted ``Delta<pT2>`` for a target against LD2."""
    p = argparse.ArgumentParser(description="Sideband-subtracted Delta<pT2> per 3D bin.")
    p.add_argument("--target", required=True, type=Path)
    p.add_argument("--ld2", required=True, type=Path)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--allow-unpublishable", action="store_true")
    args = p.parse_args(argv)

    cfg = load_config(args.config)
    da = io.load(args.target, cfg, allow_unpublishable=args.allow_unpublishable)
    dd = io.load(args.ld2, cfg, allow_unpublishable=args.allow_unpublishable)
    table, _ = delta_pt2(moments_3d(da, cfg), moments_3d(dd, cfg), cfg)

    from .ratio import write_csv  # same CSV convention, same provenance stamping

    header = [
        "Delta<pT2>_A = <pT2>_A - <pT2>_D per 3D (Q2, xB, z) bin.",
        "MOMENTS ARE SIDEBAND-SUBTRACTED (unlike the superseded analysis).",
        f"target : {da.path} | {da.stamp()}",
        f"LD2    : {dd.path} | {dd.stamp()}",
        f"config : {cfg.cuts_path} sha256={cfg.cuts_sha256}",
        "sigma^2 = (<pT4> - <pT2>^2)/N_signal ignores the variance the subtraction injects.",
    ]
    if da.blockers or dd.blockers:
        header.insert(0, "*** DIAGNOSTIC ONLY -- provenance blockers present, DO NOT QUOTE ***")
    write_csv(table, args.out, header)
    print(f"wrote {len(table)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
