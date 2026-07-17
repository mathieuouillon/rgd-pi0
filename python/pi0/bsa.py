"""Beam-spin asymmetry ``A_LU``, per 4D bin.

    A_LU(phi_h) = (1/P) * (N+ - N-) / (N+ + N-)
    sigma_A_LU  = (1/P) * sqrt( (1 - (P * A_LU)^2) / (N+ + N-) )

fitted with the **single-parameter** model ``f(phi_h) = A sin(phi_h)``,
``ndf = n_phi - 1``.

Three things this module does that the superseded analysis did not:

**1. It refuses to run without a beam polarization.**
   ``cuts.json`` carries ``bsa/polarization/value: null`` -- null meaning *not
   measured*, because no RG-D Møller number exists for these run periods. The
   old code hard-coded ``BEAM_POLARIZATION = 0.85``, self-declared as a
   placeholder in its own source, so every ``A_LU`` it published is provisional
   by a factor ``0.85 / P_true``. There is **no default here**: supply
   ``--polarization``/``--polarization-err``, or fill the config key, or get
   nothing. ``sigma_P`` is propagated (the old ``BEAM_POLARIZATION_ERR = 0.03``
   was defined and never imported anywhere).

**2. It applies the dilution correction.**
   ``A_LU^meas = A_LU^true * S/(S+B)``. Stage B bins the helicity counts in
   ``m_gg``, so ``S/(S+B)`` is measurable per bin from the very same spectra fit
   the multiplicity analysis performs, over the very same window. The old
   analysis ignored this, so every ``A_LU`` it quoted was a **lower bound** on
   its true magnitude by a factor it never measured.

**3. It masks railed fits in the producer.**
   Two leaves in the old production railed at ``|A| = 0.5`` with
   ``sigma_A ~ 4e-6`` and ``chi2/ndf ~ 1e11``; inverse-variance weighting then
   returned ``A_LU(Cu) = -0.50`` instead of ``+0.012``. The mask lives here, in
   the producer, so that no consumer can be poisoned by a rejected fit it never
   knew to look for.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import OptimizeWarning, curve_fit

from . import extract, io
from .config import BsaCuts, Config, load_config

__all__ = ["PolarizationError", "BsaResult", "compute_bsa", "main"]


class PolarizationError(RuntimeError):
    """No beam polarization was supplied, so ``A_LU`` cannot be formed."""


@dataclass(frozen=True)
class BsaResult:
    """The asymmetry of one 4D bin.

    Attributes
    ----------
    bin4d : int
        Bin index.
    amp, amp_err : float
        ``A_LU^{sin phi_h}``, **dilution-corrected** if ``dilution_applied``,
        with ``sigma_P`` folded in.
    amp_raw, amp_raw_err : float
        The fit parameter before the dilution correction -- what the old
        analysis would have quoted.
    signal_fraction, signal_fraction_err : float
        ``S/(S+B)`` in the mass window, from the spectra fit.
    dilution_applied : bool
        If False, ``amp`` is a **lower bound** on the true magnitude.
    chi2, ndf : float, int
        Fit quality. ``ndf = n_phi_used - 1`` (one free parameter).
    n_plus, n_minus : numpy.ndarray
        Helicity-summed counts per phi_h bin, inside the mass window.
    phi_used : numpy.ndarray of bool
        Which phi_h bins entered the fit.
    ok : bool
        Whether every gate and the fit-quality mask passed.
    reason : str
        ``"ok"`` or the gate it failed.
    """

    bin4d: int
    amp: float
    amp_err: float
    amp_raw: float
    amp_raw_err: float
    signal_fraction: float
    signal_fraction_err: float
    dilution_applied: bool
    chi2: float
    ndf: int
    n_plus: NDArray[np.int64]
    n_minus: NDArray[np.int64]
    phi_used: NDArray[np.bool_]
    ok: bool
    reason: str

    @property
    def chi2_ndf(self) -> float:
        """float : chi2 per degree of freedom, or ``nan``."""
        return self.chi2 / self.ndf if self.ndf > 0 else float("nan")


def _resolve_polarization(
    cuts: BsaCuts, polarization: float | None, polarization_err: float | None
) -> tuple[float, float]:
    """Resolve the beam polarization from the CLI argument or the config key.

    Raises
    ------
    PolarizationError
        If neither source supplies a value. **This function has no default.**
    """
    p = polarization if polarization is not None else cuts.polarization
    if p is None:
        raise PolarizationError(
            "No beam polarization. A_LU = (1/P) * (N+ - N-)/(N+ + N-) cannot be formed without "
            "one, and this module will not invent it.\n"
            "  cuts.json /bsa/polarization/value is null, meaning NOT MEASURED: there is no "
            "RG-D Moller number for this run period yet.\n"
            "  Supply --polarization P --polarization-err sigma_P, or fill the config key once "
            "the measurement exists.\n"
            "  The superseded analysis hard-coded P = 0.85 -- self-declared as a placeholder in "
            "its own source -- so every A_LU it published scales by 0.85/P_true. Reproducing "
            "that number here would reproduce the defect. Obtaining the measured Moller "
            "polarization is the single blocking item for this observable."
        )
    pe = polarization_err if polarization_err is not None else cuts.polarization_err
    if pe is None:
        raise PolarizationError(
            f"A beam polarization P = {p} was supplied without its uncertainty. A_LU ~ 1/P, so "
            f"sigma_P is a fully-correlated scale uncertainty on every quoted A. The old code "
            f"defined BEAM_POLARIZATION_ERR and never imported it anywhere, which is precisely "
            f"why it is mandatory here. Pass --polarization-err (0.0 is accepted, but say it)."
        )
    if not (0.0 < p <= 1.0):
        raise PolarizationError(f"beam polarization P = {p} is not in (0, 1].")
    if pe < 0.0:
        raise PolarizationError(f"beam polarization error sigma_P = {pe} is negative.")
    return float(p), float(pe)


def _sin_model(phi: NDArray[np.float64], amp: float) -> NDArray[np.float64]:
    """The single-parameter fit model ``A sin(phi_h)``, phi in radians."""
    return amp * np.sin(phi)


def compute_bsa(
    data: io.StageBData,
    cfg: Config,
    *,
    polarization: float | None = None,
    polarization_err: float | None = None,
    require_dilution: bool = True,
    stream=None,
) -> list[BsaResult]:
    """Compute ``A_LU`` for every 4D bin.

    Parameters
    ----------
    data : pi0.io.StageBData
        The loaded file.
    cfg : Config
        Configuration.
    polarization, polarization_err : float, optional
        Beam polarization and its uncertainty. If ``None``, the config key is
        used; if that is null too, :class:`PolarizationError` is raised. There
        is deliberately **no default value**.
    require_dilution : bool, default True
        If True, a bin whose ``S/(S+B)`` cannot be measured (its spectra fit
        failed) is **dropped**, rather than reported as an uncorrected lower
        bound. Set False only for diagnostics; the results are then flagged
        ``dilution_applied=False`` and are lower bounds.
    stream : file-like, optional
        Where to print the summary. Defaults to stderr.

    Returns
    -------
    list of BsaResult
        One per 4D bin, in index order. Rejected bins are present with
        ``ok=False`` and a reason.

    Raises
    ------
    PolarizationError
        If no polarization is available.

    Notes
    -----
    The mass window is the fitted +-``n_sigma_window`` window of the bin's own
    spectra fit -- the same window as the yield -- so the purity that divides
    out the dilution describes exactly the sample the helicity counts were taken
    from. (The old analysis used a fixed 0.110-0.160 GeV window and measured no
    purity at all.)

    ``sigma_A`` combines, in quadrature, the fit error, the relative error on
    ``S/(S+B)``, and ``sigma_P/P``. Two omissions are stated rather than hidden:
    ``S/(S+B)`` and ``A_fit`` are drawn from the same events and are therefore
    correlated, which this ignores; and ``sigma_f`` itself neglects the positive
    correlation between the window's subtracted yield and its raw count, so it
    is conservative.
    """
    stream = stream if stream is not None else sys.stderr
    p, p_err = _resolve_polarization(cfg.bsa, polarization, polarization_err)
    bc = cfg.bsa
    phi_rad = bc.phi.centres_rad

    # The dilution factor comes from the very same extraction the yields use.
    yields = extract.extract_all(data, cfg, check_box=False)

    out: list[BsaResult] = []
    rejected: dict[str, int] = {}

    def fail(i: int, reason: str, **kw) -> BsaResult:
        rejected[reason] = rejected.get(reason, 0) + 1
        base = dict(
            amp=np.nan, amp_err=np.nan, amp_raw=np.nan, amp_raw_err=np.nan,
            signal_fraction=np.nan, signal_fraction_err=np.nan, dilution_applied=False,
            chi2=np.nan, ndf=0,
            n_plus=np.zeros(bc.phi.n, dtype=np.int64), n_minus=np.zeros(bc.phi.n, dtype=np.int64),
            phi_used=np.zeros(bc.phi.n, dtype=bool),
        )
        base.update(kw)
        return BsaResult(bin4d=i, ok=False, reason=reason, **base)

    for i in range(cfg.binning.n_4d):
        by = yields[i]

        # --- the dilution correction, and what to do when it is unavailable
        f_s, f_s_err, dilution_applied = np.nan, np.nan, False
        if by.window is None or not np.isfinite(by.value):
            if require_dilution:
                out.append(fail(i, "no_dilution_available"))
                continue
            window = cfg.mgg.mask(cfg.extraction.fit_range_min, cfg.extraction.fit_range_max)
        else:
            window = by.window
            f_s = by.signal_fraction
            if np.isfinite(f_s) and f_s > 0:
                # sigma_f/f neglects the (positive) correlation between the
                # subtracted yield and the raw window count, so it is conservative.
                f_s_err = f_s * float(
                    np.sqrt((by.error / by.value) ** 2 + 1.0 / by.n_same_window)
                )
                dilution_applied = True
            elif require_dilution:
                out.append(fail(i, "non_positive_signal_fraction"))
                continue

        n_plus = data.bsa[i, window, :, 0].sum(axis=0).astype(np.int64)
        n_minus = data.bsa[i, window, :, 1].sum(axis=0).astype(np.int64)
        n_tot = n_plus + n_minus

        if float(n_tot.sum()) < bc.min_counts_per_bin:
            out.append(fail(i, "below_min_counts_per_bin", n_plus=n_plus, n_minus=n_minus))
            continue
        used = n_tot >= bc.min_counts_per_phi_bin
        if int(used.sum()) < bc.min_phi_bins:
            out.append(fail(i, "too_few_phi_bins", n_plus=n_plus, n_minus=n_minus))
            continue

        r = (n_plus[used] - n_minus[used]) / n_tot[used]  # r = P * A_LU
        a_lu = r / p
        # sigma = (1/P) sqrt((1 - r^2)/N): the exact binomial variance of r,
        # well defined for all |r| <= 1.
        sig = np.sqrt(np.clip(1.0 - r**2, 0.0, None) / n_tot[used]) / p
        if np.any(sig <= 0):
            # A bin where |r| == 1 has zero estimated variance and would carry
            # infinite weight. That is a statistics artefact, not a measurement.
            out.append(fail(i, "zero_variance_phi_bin", n_plus=n_plus, n_minus=n_minus))
            continue

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", OptimizeWarning)
                popt, pcov = curve_fit(
                    _sin_model, phi_rad[used], a_lu, p0=[0.0], sigma=sig,
                    absolute_sigma=True, bounds=([bc.amp_min], [bc.amp_max]), maxfev=10000,
                )
        except (RuntimeError, ValueError, OptimizeWarning, TypeError):
            out.append(fail(i, "fit_did_not_converge", n_plus=n_plus, n_minus=n_minus))
            continue
        if not np.all(np.isfinite(pcov)):
            out.append(fail(i, "fit_covariance_not_finite", n_plus=n_plus, n_minus=n_minus))
            continue

        amp_raw = float(popt[0])
        amp_raw_err = float(np.sqrt(pcov[0, 0]))
        chi2 = float(np.sum(((a_lu - _sin_model(phi_rad[used], amp_raw)) / sig) ** 2))
        ndf = int(used.sum()) - 1  # ONE free parameter

        # ---- the fit-quality mask, in the producer ----
        bound = max(abs(bc.amp_min), abs(bc.amp_max))
        common = dict(
            amp_raw=amp_raw, amp_raw_err=amp_raw_err, chi2=chi2, ndf=ndf,
            n_plus=n_plus, n_minus=n_minus, phi_used=used,
            signal_fraction=f_s, signal_fraction_err=f_s_err,
        )
        if abs(amp_raw) >= bound * (1.0 - bc.amp_at_bound_tol_frac):
            # The minimiser walked into the wall: its error is meaningless and
            # inverse-variance weighting downstream would be dominated by it.
            out.append(fail(i, "amp_at_bound", **common))
            continue
        if amp_raw_err < bc.sigma_amp_min:
            out.append(fail(i, "sigma_amp_below_floor", **common))
            continue
        if ndf > 0 and chi2 / ndf > bc.max_chi2_ndf:
            out.append(fail(i, "chi2_ndf_too_large", **common))
            continue

        if dilution_applied:
            amp = amp_raw / f_s
            rel = (amp_raw_err / abs(amp_raw)) ** 2 if amp_raw else np.inf
            rel += (f_s_err / f_s) ** 2 + (p_err / p) ** 2
            amp_err = abs(amp) * float(np.sqrt(rel))
        else:
            amp, amp_err = amp_raw, float(np.hypot(amp_raw_err, abs(amp_raw) * p_err / p))

        out.append(
            BsaResult(
                bin4d=i, amp=amp, amp_err=amp_err, dilution_applied=dilution_applied,
                ok=True, reason=extract.OK, **common
            )
        )

    n_ok = sum(1 for x in out if x.ok)
    print(f"bsa: {n_ok} of {cfg.binning.n_4d} bins survived (P = {p} +- {p_err}).", file=stream)
    if rejected:
        for code, n in sorted(rejected.items(), key=lambda kv: -kv[1]):
            print(f"    {n:6d}  {code}", file=stream)
    n_undiluted = sum(1 for x in out if x.ok and not x.dilution_applied)
    if n_undiluted:
        print(
            f"\n    *** {n_undiluted} surviving bins have NO DILUTION CORRECTION. Their A_LU is a\n"
            f"    *** LOWER BOUND on the true magnitude, by an unmeasured and not necessarily\n"
            f"    *** target-independent factor. Do not average them with corrected bins.\n",
            file=stream,
        )
    return out


def main(argv: list[str] | None = None) -> int:
    """CLI: compute ``A_LU`` per 4D bin. Refuses to run without a polarization."""
    p = argparse.ArgumentParser(
        description="Beam-spin asymmetry A_LU per 4D bin. Requires an explicit beam polarization."
    )
    p.add_argument("--file", required=True, type=Path, help="Stage B file")
    p.add_argument("--config", required=True, type=Path, help="config/ directory")
    p.add_argument("--out", required=True, type=Path, help="output CSV")
    p.add_argument(
        "--polarization", type=float, default=None,
        help="measured beam polarization P. MANDATORY unless cuts.json /bsa/polarization/value "
             "is set. There is no default -- 0.85 was the old code's self-declared placeholder.",
    )
    p.add_argument(
        "--polarization-err", type=float, default=None,
        help="uncertainty sigma_P on the beam polarization. Mandatory alongside --polarization; "
             "it is a fully-correlated scale uncertainty on every A.",
    )
    p.add_argument(
        "--allow-uncorrected-dilution", action="store_true",
        help="report bins whose S/(S+B) could not be measured. Their A_LU is a LOWER BOUND.",
    )
    p.add_argument("--allow-unpublishable", action="store_true")
    args = p.parse_args(argv)

    cfg = load_config(args.config)
    data = io.load(args.file, cfg, allow_unpublishable=args.allow_unpublishable)
    try:
        res = compute_bsa(
            data, cfg,
            polarization=args.polarization,
            polarization_err=args.polarization_err,
            require_dilution=not args.allow_uncorrected_dilution,
        )
    except PolarizationError as exc:
        print(f"\nREFUSING TO COMPUTE A_LU:\n{exc}\n", file=sys.stderr)
        return 2

    dtype = np.dtype(
        [
            ("bin4d", np.int64), ("i_q2", np.int64), ("i_xb", np.int64),
            ("i_z", np.int64), ("i_pt2", np.int64),
            ("A_LU", np.float64), ("sA_LU", np.float64),
            ("A_raw", np.float64), ("sA_raw", np.float64),
            ("signal_fraction", np.float64), ("dilution_applied", np.int64),
            ("chi2", np.float64), ("ndf", np.int64), ("n_events", np.int64),
        ]
    )
    rows = []
    for r in res:
        if not r.ok:
            continue
        i_q2, i_xb, i_z, i_pt2 = cfg.binning.decode_4d(r.bin4d)
        rows.append(
            (r.bin4d, i_q2, i_xb, i_z, i_pt2, r.amp, r.amp_err, r.amp_raw, r.amp_raw_err,
             r.signal_fraction, int(r.dilution_applied), r.chi2, r.ndf,
             int(r.n_plus.sum() + r.n_minus.sum()))
        )
    table = np.array(rows, dtype=dtype)

    from .ratio import write_csv

    header = [
        "A_LU^{sin phi_h} per 4D bin, from f(phi_h) = A sin(phi_h), ndf = n_phi - 1.",
        f"file   : {data.path} | {data.stamp()}",
        f"config : {cfg.cuts_path} sha256={cfg.cuts_sha256}",
        f"beam polarization P = {args.polarization or cfg.bsa.polarization} "
        f"+- {args.polarization_err if args.polarization_err is not None else cfg.bsa.polarization_err}",
        "A_LU IS DILUTION-CORRECTED (A_meas = A_true * S/(S+B)) where dilution_applied = 1.",
        "A_raw is the uncorrected fit parameter -- what the superseded analysis quoted.",
        "Railed fits are rejected in the producer (|A| at bound, sigma_A < floor, chi2/ndf > max).",
    ]
    if data.blockers:
        header.insert(0, "*** DIAGNOSTIC ONLY -- provenance blockers present, DO NOT QUOTE ***")
    write_csv(table, args.out, header)
    print(f"wrote {len(table)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
