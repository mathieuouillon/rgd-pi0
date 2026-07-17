"""Yield extraction and count-weighted abscissae, per 4D bin.

The procedure is the note's ``sec:yield-extraction``, reproduced exactly:

1. Scale the mixed spectrum to the same-event one in a single **high-mass**
   sideband, ``alpha = sum(n_same in SB) / sum(n_mixed in SB)``.
2. Subtract bin by bin: ``S(m) = n_same(m) - alpha * n_mixed(m)``, with
   ``var(m) = n_same(m) + alpha^2 * n_mixed(m)``.
3. Fit ``S(m)`` with a **pure Gaussian** -- no background term, the background
   is already gone -- over the configured fit range.
4. Iteratively restrict to the core: refit over ``[mu - k*sigma, mu + k*sigma]``.
5. **The yield is the +-n_sigma window SUM of the subtracted histogram**, not
   the Gaussian integral. The fit exists only to locate the window.

.. important::
   Step 5 is not a detail. The superseded script's own docstring claimed
   ``N = amp * sigma * sqrt(2 pi) / binW`` with a covariance-propagated error.
   **That is not what produced its results** -- the analytic integral existed
   only in a pedagogical demonstration plot. This module reproduces the *code*,
   not the docstring, and :func:`gaussian_integral_yield` exists solely so the
   test suite can prove the two differ and pin which one is used.

Known omissions in the error model, carried deliberately (note
``sec:ra-uncertainties``, ``sec:subtraction``):

* the uncertainty on ``alpha`` is **not** propagated -- the subtraction treats
  it as exact;
* the window sum's error ignores the uncertainty on the fitted ``mu`` and
  ``sigma`` that *define* the window.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import OptimizeWarning, curve_fit

from .config import Config, ExtractionCuts, MggAxis

__all__ = [
    "FitResult",
    "BinYield",
    "gaussian",
    "sideband_scale",
    "subtract",
    "fit_core",
    "extract_bin",
    "extract_all",
    "abscissa",
    "gaussian_integral_yield",
]

# Reject codes. Every dropped bin carries one; nothing is dropped silently.
OK = "ok"
NO_SIDEBAND_MIXED = "no_mixed_in_sideband"
NO_SIDEBAND_SAME = "no_same_in_sideband"
LOW_STATS = "below_min_stats"
LOW_PEAK = "below_min_peak_amp"
FIT_FAILED = "fit_did_not_converge"
EMPTY_WINDOW = "empty_3sigma_window"
NO_ABSCISSA = "abscissa_not_computable"
ABSCISSA_OUT_OF_BOX = "abscissa_outside_own_box"


def gaussian(m: NDArray[np.float64], amp: float, mu: float, sigma: float) -> NDArray[np.float64]:
    """A pure Gaussian ``amp * exp(-(m - mu)^2 / (2 sigma^2))``.

    No background term: by the time this is fitted, the background has already
    been subtracted.
    """
    return amp * np.exp(-0.5 * ((m - mu) / sigma) ** 2)


@dataclass(frozen=True)
class FitResult:
    """Outcome of the iterative core-restricted Gaussian fit.

    Attributes
    ----------
    amp, mu, sigma : float
        Fitted parameters. ``amp`` never enters the yield.
    amp_err, mu_err, sigma_err : float
        Parameter errors from ``curve_fit(absolute_sigma=True)``, so they are
        *not* rescaled by chi2/ndf.
    chi2, ndf : float, int
        Fit quality over the final (core-restricted) range.
    n_iter : int
        How many core-restriction iterations actually ran.
    converged : bool
        False if any iteration raised or the covariance was not finite.
    range_lo, range_hi : float
        The final fit range in GeV.
    reason : str
        ``"ok"``, or why it failed.
    """

    amp: float
    mu: float
    sigma: float
    amp_err: float
    mu_err: float
    sigma_err: float
    chi2: float
    ndf: int
    n_iter: int
    converged: bool
    range_lo: float
    range_hi: float
    reason: str = OK

    @property
    def chi2_ndf(self) -> float:
        """float : chi2 per degree of freedom, or ``nan`` if ndf <= 0."""
        return self.chi2 / self.ndf if self.ndf > 0 else float("nan")


@dataclass(frozen=True)
class BinYield:
    """Everything the extraction knows about one 4D bin.

    Attributes
    ----------
    bin4d : int
        The bin index.
    alpha : float
        The sideband scale.
    fit : FitResult or None
        The mass fit. ``None`` if the bin failed before fitting.
    value, error : float
        The yield ``Y`` and ``sigma_Y``: the +-n_sigma window **sum** of the
        subtracted histogram, and the root of the summed per-bin variance.
    n_same_window, n_bkg_window : float
        Raw same-event counts in the window, and the estimated background under
        it (``alpha * sum(n_mixed)``). Their difference is ``value``.
    signal_fraction : float
        ``value / n_same_window`` -- the purity ``S / (S + B)`` in the window.
        This is what :mod:`pi0.bsa` needs for its dilution correction.
    window : numpy.ndarray or None
        Boolean mask over mgg bins selecting the window.
    abscissae : dict of str to float
        Count-weighted, sideband-subtracted ``<q2>``, ``<xb>``, ``<z>``,
        ``<pt2>``. Empty if they could not be computed.
    ok : bool
        Whether the bin passed every gate.
    reason : str
        ``"ok"``, or the first gate it failed.
    """

    bin4d: int
    alpha: float
    fit: FitResult | None
    value: float
    error: float
    n_same_window: float
    n_bkg_window: float
    window: NDArray[np.bool_] | None
    abscissae: dict[str, float]
    ok: bool
    reason: str

    @property
    def signal_fraction(self) -> float:
        """float : Purity ``S / (S + B)`` in the mass window, or ``nan``."""
        if self.n_same_window <= 0:
            return float("nan")
        return self.value / self.n_same_window


def sideband_scale(
    n_same: NDArray, n_mixed: NDArray, sb: NDArray[np.bool_]
) -> tuple[float, str]:
    """Compute the sideband scale ``alpha``.

    ``alpha = sum(n_same in SB) / sum(n_mixed in SB)`` over a single high-mass
    sideband (note ``eq:alpha``). The low-mass sideband is deliberately unused:
    there the spectrum is distorted by the opening-angle cut, which removes
    precisely the small-angle pairs that populate it, so same and mixed are not
    expected to share a shape.

    Parameters
    ----------
    n_same, n_mixed : numpy.ndarray
        Same- and mixed-event counts per mgg bin.
    sb : numpy.ndarray of bool
        Mask selecting the sideband bins.

    Returns
    -------
    alpha : float
        The scale, or ``nan`` if it is undefined.
    reason : str
        ``"ok"``, or why alpha is undefined.

    Notes
    -----
    The uncertainty on ``alpha`` is not returned and not propagated anywhere
    downstream -- a real omission, inherited deliberately from the note's
    ``sec:subtraction``, and trivially fixable if it ever matters.
    """
    den = float(np.sum(n_mixed[sb]))
    num = float(np.sum(n_same[sb]))
    if den <= 0:
        return float("nan"), NO_SIDEBAND_MIXED
    if num <= 0:
        return float("nan"), NO_SIDEBAND_SAME
    return num / den, OK


def subtract(
    n_same: NDArray, n_mixed: NDArray, alpha: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Sideband-scaled subtraction, bin by bin (note ``eq:subtraction``).

    Returns
    -------
    S : numpy.ndarray
        ``n_same - alpha * n_mixed``.
    var : numpy.ndarray
        ``n_same + alpha^2 * n_mixed``. Poisson on both terms; ``alpha`` is
        treated as exact.
    """
    s = n_same.astype(np.float64) - alpha * n_mixed.astype(np.float64)
    var = n_same.astype(np.float64) + alpha**2 * n_mixed.astype(np.float64)
    return s, var


def _fit_once(
    m: NDArray[np.float64],
    s: NDArray[np.float64],
    var: NDArray[np.float64],
    lo: float,
    hi: float,
    cuts: ExtractionCuts,
    seed: tuple[float, float, float],
) -> tuple[tuple[float, float, float], tuple[float, float, float], float, int, bool]:
    """One pure-Gaussian fit over ``[lo, hi]``. Returns (pars, errs, chi2, ndf, ok)."""
    # Bins with var == 0 carry no information but would get infinite weight.
    # They are empty bins of an empty spectrum; drop them rather than let them
    # dominate the chi2.
    sel = (m >= lo) & (m <= hi) & (var > 0)
    n = int(sel.sum())
    if n < 4:  # 3 free parameters; a fit with ndf <= 0 is not a measurement
        return seed, (np.nan,) * 3, np.nan, n - 3, False
    x, y, sy = m[sel], s[sel], np.sqrt(var[sel])
    bounds = (
        [cuts.amp_min, cuts.mu_min, cuts.sigma_min],
        [cuts.amp_max, cuts.mu_max, cuts.sigma_max],
    )
    p0 = [
        float(np.clip(seed[0], cuts.amp_min, cuts.amp_max)),
        float(np.clip(seed[1], cuts.mu_min, cuts.mu_max)),
        float(np.clip(seed[2], cuts.sigma_min, cuts.sigma_max)),
    ]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            popt, pcov = curve_fit(
                gaussian, x, y, p0=p0, sigma=sy, absolute_sigma=True, bounds=bounds, maxfev=20000
            )
    except (RuntimeError, ValueError, OptimizeWarning, TypeError):
        return seed, (np.nan,) * 3, np.nan, n - 3, False
    if not np.all(np.isfinite(pcov)):
        return tuple(popt), (np.nan,) * 3, np.nan, n - 3, False
    resid = (y - gaussian(x, *popt)) / sy
    chi2 = float(np.sum(resid**2))
    errs = tuple(np.sqrt(np.diag(pcov)))
    return tuple(popt), errs, chi2, n - 3, True


def fit_core(
    axis: MggAxis,
    s: NDArray[np.float64],
    var: NDArray[np.float64],
    cuts: ExtractionCuts,
) -> FitResult:
    """Fit the subtracted spectrum with a pure Gaussian, restricting to the core.

    Seeded at ``(max S over the fit range, mu_seed, sigma_seed)``, fitted over
    ``[fit_range_min, fit_range_max]``, then refitted over
    ``[mu - core_k*sigma, mu + core_k*sigma]`` for up to ``n_iter`` iterations,
    stopping as soon as ``mu`` **and** ``sigma`` both change by less than
    ``convergence_tol_frac``.

    Parameters
    ----------
    axis : MggAxis
        The mass axis, supplying bin centres.
    s, var : numpy.ndarray
        The subtracted spectrum and its per-bin variance.
    cuts : ExtractionCuts
        Fit range, bounds, seeds, ``core_k``, ``n_iter``, tolerance -- all from
        ``cuts.json``.

    Returns
    -------
    FitResult
    """
    m = axis.centres
    lo, hi = cuts.fit_range_min, cuts.fit_range_max
    in_range = (m >= lo) & (m <= hi)
    amp_seed = float(np.max(s[in_range])) if in_range.any() else 0.0
    seed = (amp_seed, cuts.mu_seed, cuts.sigma_seed)

    pars, errs, chi2, ndf, ok = _fit_once(m, s, var, lo, hi, cuts, seed)
    if not ok:
        return FitResult(*pars, *errs, chi2, ndf, 0, False, lo, hi, FIT_FAILED)

    n_done = 0
    for _ in range(cuts.n_iter):
        mu, sigma = pars[1], pars[2]
        new_lo, new_hi = mu - cuts.core_k * sigma, mu + cuts.core_k * sigma
        p2, e2, c2, n2, ok2 = _fit_once(m, s, var, new_lo, new_hi, cuts, pars)
        if not ok2:
            # The core refit failed; keep the wider fit, which did converge.
            break
        n_done += 1
        pars, errs, chi2, ndf, lo, hi = p2, e2, c2, n2, new_lo, new_hi
        d_mu = abs(pars[1] - mu) / abs(mu) if mu else np.inf
        d_sig = abs(pars[2] - sigma) / abs(sigma) if sigma else np.inf
        if d_mu < cuts.convergence_tol_frac and d_sig < cuts.convergence_tol_frac:
            break

    return FitResult(*pars, *errs, chi2, ndf, n_done, True, lo, hi, OK)


def abscissa(
    sum_x: NDArray[np.float64],
    n_same: NDArray,
    n_mixed: NDArray,
    alpha: float,
    window: NDArray[np.bool_],
    sb: NDArray[np.bool_],
) -> float:
    """Count-weighted, sideband-subtracted ``<X>`` for one bin. **The point of the rewrite.**

    Never returns a geometric bin centre. If the estimate is not computable the
    caller must drop the bin: a bin reported at the wrong place is worse than a
    missing bin (note ``sec:binning-caveat``).

    The estimator
    -------------
    Stage B's ``sum_*`` columns are **same-event only** -- no mixed pair
    contributes to any of them, because a mixed pair is not a pi0 and its
    kinematics describe nothing. So there is no mixed-event ``<X>`` to subtract
    directly, and the background's share must be estimated from the sideband.

    Decompose the same-event pairs in mass bin ``m`` into true pi0 and
    combinatorial background::

        sum_X(m) = S(m) <X>_S(m) + B(m) <X>_B(m),    B(m) = alpha n_mixed(m)

    Summed over the peak window ``W``, and **assuming <X>_B varies slowly with
    m_gg** so that it may be taken constant across the window and the
    sideband::

        <X>_B  ~=  sum_X(SB) / n_same(SB)          (the sideband is ~pure background)
        B_W     =  alpha * sum(n_mixed over W)
        <X>_S  ~=  (sum_X(W) - <X>_B * B_W) / (n_same(W) - B_W)

    which is what this function computes.

    .. warning::
       **This is an approximation, and its assumption is the weak point.** It is
       exact only if ``<X>_B`` is independent of ``m_gg`` between the peak
       (~0.135 GeV) and the sideband (0.17-0.28 GeV). For ``z`` in particular
       that is not obviously true: ``m_gg`` correlates with the pair's energy at
       fixed opening angle, and ``z = E_pi0 / nu``, so the background's ``<z>``
       plausibly rises with ``m_gg`` and a flat extrapolation down to the peak
       would then over-subtract. The bias is *not* propagated into any error
       here. :func:`pi0.qa.sideband_slope_diagnostic` measures the trend inside
       the sideband so that the size of the assumption is at least visible.

       **The exactly-correct version needs Stage B to write a mixed-event**
       ``sum_X`` **per (bin4d, imgg)**, i.e. ``sum_q2_mixed``, ``sum_xb_mixed``,
       ``sum_z_mixed``, ``sum_pt2_mixed`` alongside ``n_mixed``, with each mixed
       pair's kinematics computed against its *current* event's DIS electron
       (which is what the same-event background it models is made of). Then::

           <X>_S = (sum_X_same(W) - alpha * sum_X_mixed(W)) / (n_same(W) - alpha * n_mixed(W))

       with no flatness assumption at all -- exactly the structure already used
       for the counts. That is a small Stage B change and it is the right fix.

    Parameters
    ----------
    sum_x : numpy.ndarray
        Same-event ``sum_X`` per mgg bin.
    n_same, n_mixed : numpy.ndarray
        Counts per mgg bin.
    alpha : float
        Sideband scale.
    window : numpy.ndarray of bool
        The +-n_sigma peak window.
    sb : numpy.ndarray of bool
        The sideband mask.

    Returns
    -------
    float
        ``<X>`` of the signal, or ``nan`` if it is not computable (empty
        sideband, or a non-positive signal count under the peak).
    """
    n_same_sb = float(np.sum(n_same[sb]))
    if n_same_sb <= 0:
        return float("nan")
    x_bkg = float(np.sum(sum_x[sb])) / n_same_sb

    n_bkg_w = alpha * float(np.sum(n_mixed[window]))
    n_sig_w = float(np.sum(n_same[window])) - n_bkg_w
    if n_sig_w <= 0:
        return float("nan")

    num = float(np.sum(sum_x[window])) - x_bkg * n_bkg_w
    val = num / n_sig_w
    return val if np.isfinite(val) else float("nan")


def gaussian_integral_yield(fit: FitResult, axis: MggAxis) -> float:
    """``amp * sigma * sqrt(2 pi) / bin_width`` -- **NOT** the yield.

    This is the formula the superseded script's docstring claimed produced its
    results. It did not: that path existed only in a demonstration plot, while
    ``R_A`` came from the window sum. It is implemented here for exactly one
    reason -- so ``tests/test_extract.py`` can inject a non-Gaussian tail and
    assert that :func:`extract_bin` returns the window sum and *not* this.

    Do not call it from an analysis path.
    """
    return fit.amp * fit.sigma * np.sqrt(2.0 * np.pi) / axis.width


def extract_bin(
    bin4d: int,
    n_same: NDArray,
    n_mixed: NDArray,
    sums: dict[str, NDArray[np.float64]],
    cfg: Config,
    *,
    shared_fit: FitResult | None = None,
    box: dict[str, tuple[float, float]] | None = None,
) -> BinYield:
    """Run the full extraction on one 4D bin's spectra.

    Parameters
    ----------
    bin4d : int
        Bin index, carried through to the result.
    n_same, n_mixed : numpy.ndarray
        This bin's same- and mixed-event mgg spectra.
    sums : dict of str to numpy.ndarray
        This bin's same-event ``sum_X`` spectra, keyed by axis name.
    cfg : Config
        Configuration. Every cut value comes from here.
    shared_fit : FitResult, optional
        If given, skip this bin's own fit and take ``mu``/``sigma`` from the
        supplied fit. **This is not the production path**: it exists for
        low-statistics diagnostics, where a per-bin fit is meaningless but the
        abscissa machinery can still be exercised against a window located on
        the summed spectrum. Callers must label any output produced this way.
    box : dict, optional
        This bin's ``(lo, hi)`` box per axis. If given, an abscissa falling
        outside its own box is treated as not computable and the bin is
        dropped -- the true mean of a sample confined to a box lies inside it,
        so a value outside is a broken estimate, and a bin at the wrong place
        is worse than a missing bin.

    Returns
    -------
    BinYield
    """
    axis, cuts = cfg.mgg, cfg.extraction
    sb = axis.mask(cuts.sideband_min, cuts.sideband_max)
    fit_mask = axis.mask(cuts.fit_range_min, cuts.fit_range_max)

    def fail(reason: str, alpha: float = float("nan")) -> BinYield:
        return BinYield(bin4d, alpha, None, np.nan, np.nan, np.nan, np.nan, None, {}, False, reason)

    # Gate 1: raw same-event statistics, over the fit range (the region the
    # extraction actually uses). Applied to the RAW counts, so the gate never
    # depends on the subtraction it is gating.
    if float(np.sum(n_same[fit_mask])) < cuts.min_stats:
        return fail(LOW_STATS)

    alpha, why = sideband_scale(n_same, n_mixed, sb)
    if not np.isfinite(alpha):
        return fail(why)

    s, var = subtract(n_same, n_mixed, alpha)

    # Gate 2: subtracted peak amplitude.
    if float(np.max(s[fit_mask])) < cuts.min_peak_amp:
        return fail(LOW_PEAK, alpha)

    if shared_fit is not None:
        fit = shared_fit
    else:
        fit = fit_core(axis, s, var, cuts)
        # Gate 3: the fit must converge.
        if not fit.converged:
            return BinYield(bin4d, alpha, fit, np.nan, np.nan, np.nan, np.nan, None, {}, False,
                            FIT_FAILED)

    # THE YIELD: the +-n_sigma window SUM of the SUBTRACTED histogram.
    # Not the Gaussian integral. The fit supplied mu and sigma and nothing else.
    half = cuts.n_sigma_window * fit.sigma
    window = np.abs(axis.centres - fit.mu) < half
    if not window.any():
        return BinYield(bin4d, alpha, fit, np.nan, np.nan, np.nan, np.nan, window, {}, False,
                        EMPTY_WINDOW)

    value = float(np.sum(s[window]))
    error = float(np.sqrt(np.sum(var[window])))
    n_same_w = float(np.sum(n_same[window]))
    n_bkg_w = alpha * float(np.sum(n_mixed[window]))

    absc: dict[str, float] = {}
    for name, sum_x in sums.items():
        absc[name] = abscissa(sum_x, n_same, n_mixed, alpha, window, sb)

    reason, ok = OK, True
    if not all(np.isfinite(v) for v in absc.values()):
        reason, ok = NO_ABSCISSA, False
    elif box is not None:
        for name, val in absc.items():
            lo, hi = box[name]
            if not (lo <= val <= hi):
                reason, ok = ABSCISSA_OUT_OF_BOX, False
                break

    return BinYield(bin4d, alpha, fit, value, error, n_same_w, n_bkg_w, window, absc, ok, reason)


def extract_all(
    data,
    cfg: Config,
    *,
    shared_fit: FitResult | None = None,
    check_box: bool = True,
) -> list[BinYield]:
    """Run :func:`extract_bin` over every 4D bin of a Stage B file.

    Parameters
    ----------
    data : pi0.io.StageBData
        The loaded file.
    cfg : Config
        Configuration.
    shared_fit : FitResult, optional
        See :func:`extract_bin`. Diagnostic use only.
    check_box : bool, default True
        Whether to drop bins whose abscissa falls outside their own box.

    Returns
    -------
    list of BinYield
        One entry per 4D bin, in index order. Failed bins are present with
        ``ok=False`` and a reason -- nothing is dropped silently.
    """
    b = cfg.binning
    sums_all = data.sums
    out: list[BinYield] = []
    for i in range(b.n_4d):
        box = _box_of(b, i) if check_box else None
        out.append(
            extract_bin(
                i,
                data.n_same[i],
                data.n_mixed[i],
                {k: v[i] for k, v in sums_all.items()},
                cfg,
                shared_fit=shared_fit,
                box=box,
            )
        )
    return out


def _box_of(binning, bin4d: int) -> dict[str, tuple[float, float]]:
    """The ``(lo, hi)`` edges of a 4D bin on each of its four axes."""
    i_q2, i_xb, i_z, i_pt2 = binning.decode_4d(bin4d)
    q2, xb = binning.grid_a.axes
    z, pt2 = binning.grid_b.axes
    return {
        "q2": (float(q2.edges[i_q2]), float(q2.edges[i_q2 + 1])),
        "xb": (float(xb.edges[i_xb]), float(xb.edges[i_xb + 1])),
        "z": (float(z.edges[i_z]), float(z.edges[i_z + 1])),
        "pt2": (float(pt2.edges[i_pt2]), float(pt2.edges[i_pt2 + 1])),
    }


def summed_spectrum(data, cfg: Config, bins: Sequence[int] | None = None):
    """Sum every 4D bin's spectra into one, for a global diagnostic fit.

    Useful when per-bin statistics cannot support a fit (a smoke file), and for
    the QA overview plot. The result is **not** a physics measurement: summing
    over kinematics mixes bins with different peak positions and resolutions.

    Returns
    -------
    tuple
        ``(n_same, n_mixed, sums)`` with the same shapes as one bin's.
    """
    sel = slice(None) if bins is None else np.asarray(bins)
    n_same = data.n_same[sel].sum(axis=0)
    n_mixed = data.n_mixed[sel].sum(axis=0)
    sums = {k: v[sel].sum(axis=0) for k, v in data.sums.items()}
    return n_same, n_mixed, sums
