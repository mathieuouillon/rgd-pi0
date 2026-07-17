"""Configuration and binning geometry for the RG-D pi0 extraction stage.

This module is the **only** place in the Python package that knows

* where a cut value comes from (always ``config/cuts.json``), and
* how a bin index is encoded or decoded.

Nothing else may hard-code either. Every lookup goes through :func:`require`,
which raises :class:`ConfigError` on a missing key rather than defaulting --
a silently defaulted cut is the failure mode this project exists to end.

Index formulas
--------------
Reproduced here verbatim from Stage B (``provenance/binning.index_formula_*``),
because the superseded analysis's leaf formula was written down *nowhere* and
had to be reverse-engineered out of its output files::

    cell_a = i_q2 * n_xb  + i_xb          index into Grid A (Q2 x xB)
    cell_b = i_z  * n_pt2 + i_pt2         index into Grid B (z x pT2)
    bin4d  = cell_a * (n_z * n_pt2) + cell_b
    bin3d  = cell_a * n_z + i_z           NOT bin4d / n_pt2
    mgg bin centre = (imgg + 0.5) * (mgg_max - mgg_min) / n_mgg + mgg_min

Axis semantics (``provenance/axis.semantics``): every axis is half-open
``[lo, hi)``, a value exactly on the top edge lands in the last bin, and
out-of-range or NaN gives ``-1``. There are no flow bins.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "ConfigError",
    "Grid1D",
    "Grid2D",
    "Binning",
    "MggAxis",
    "PhiAxis",
    "ExtractionCuts",
    "BsaCuts",
    "BroadeningCuts",
    "Config",
    "load_config",
    "require",
]


class ConfigError(KeyError):
    """A required configuration key is missing, null, or of the wrong type."""


def require(node: Any, path: str, *, expect: type | tuple[type, ...] | None = None) -> Any:
    """Fetch a nested key from a parsed JSON document, or fail loudly.

    Parameters
    ----------
    node : Any
        The parsed JSON document (or sub-document) to look in.
    path : str
        Slash- or dot-separated key path, e.g. ``"extraction/min_stats"``.
    expect : type or tuple of type, optional
        If given, the value's type is checked against it.

    Returns
    -------
    Any
        The value at ``path``.

    Raises
    ------
    ConfigError
        If any component of ``path`` is absent, or the value is ``None``, or
        it fails the ``expect`` type check. Never returns a default.
    """
    parts = path.replace(".", "/").split("/")
    cur = node
    for i, key in enumerate(parts):
        if not isinstance(cur, dict) or key not in cur:
            so_far = "/".join(parts[:i]) or "<root>"
            avail = sorted(k for k in cur) if isinstance(cur, dict) else "<not a mapping>"
            raise ConfigError(
                f"required config key {path!r} is missing: no {key!r} under {so_far}. "
                f"Available there: {avail}"
            )
        cur = cur[key]
    if cur is None:
        raise ConfigError(
            f"required config key {path!r} is present but null. It must be given an "
            f"explicit value; this stage does not substitute a default for it."
        )
    if expect is not None and not isinstance(cur, expect):
        raise ConfigError(f"config key {path!r} has type {type(cur).__name__}, expected {expect}")
    return cur


def _sha256(path: Path) -> str:
    """Return the hex SHA-256 of a file's contents."""
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


@dataclass(frozen=True)
class Grid1D:
    """A one-dimensional binning axis with explicit edges.

    Mirrors ``pi0::Grid1D`` in the C++: half-open ``[lo, hi)``, top edge folded
    into the last bin, ``-1`` for out-of-range or NaN, no flow bins.

    Attributes
    ----------
    name : str
        Axis name as given in the grid JSON (``q2``, ``xb``, ``z``, ``pt2``).
    edges : numpy.ndarray
        Strictly increasing edge array of length ``n + 1``.
    """

    name: str
    edges: NDArray[np.float64]

    @property
    def n(self) -> int:
        """int : Number of bins on this axis."""
        return len(self.edges) - 1

    @property
    def centres(self) -> NDArray[np.float64]:
        """numpy.ndarray : Geometric midpoints of each bin.

        Provided for the QA comparison plot **only**. No result may be reported
        at these positions -- see :mod:`pi0.extract` and the note's
        ``sec:binning-caveat``.
        """
        return 0.5 * (self.edges[:-1] + self.edges[1:])

    @property
    def widths(self) -> NDArray[np.float64]:
        """numpy.ndarray : Width of each bin."""
        return np.diff(self.edges)

    def find(self, x: float) -> int:
        """Return the bin index containing ``x``, or ``-1`` if out of range."""
        if not np.isfinite(x) or x < self.edges[0] or x > self.edges[-1]:
            return -1
        return int(min(np.searchsorted(self.edges, x, side="right") - 1, self.n - 1))


@dataclass(frozen=True)
class Grid2D:
    """A product of two :class:`Grid1D` axes, as stored in one grid JSON file.

    Attributes
    ----------
    path : pathlib.Path
        File the grid was read from.
    axes : tuple of Grid1D
        The two axes, in the file's declared order (slow axis first).
    source : str
        ``provenance.source`` from the file. ``"placeholder"`` means the edges
        are hand-authored, not fitted to data, and **no result may be quoted
        from them**.
    n_events : int
        ``provenance.n_events`` -- how many events the quantiles were computed
        from. ``0`` for a placeholder grid.
    """

    path: Path
    axes: tuple[Grid1D, Grid1D]
    source: str
    n_events: int

    @property
    def is_placeholder(self) -> bool:
        """bool : True if these edges were never fitted to data."""
        return self.source.strip().lower() == "placeholder" or self.n_events == 0

    @property
    def n(self) -> int:
        """int : Total number of cells (product of the two axes)."""
        return self.axes[0].n * self.axes[1].n


def _load_grid(path: Path, expect_axes: Sequence[str]) -> Grid2D:
    """Read a grid JSON file, validating axis names and edge monotonicity."""
    doc = json.loads(path.read_text())
    axes_doc = require(doc, "axes")
    if len(axes_doc) != 2:
        raise ConfigError(f"{path}: expected exactly 2 axes, found {len(axes_doc)}")
    axes: list[Grid1D] = []
    for ax_doc, want in zip(axes_doc, expect_axes):
        name = require(ax_doc, "name", expect=str)
        if name != want:
            raise ConfigError(f"{path}: axis order is {name!r}, expected {want!r}")
        edges = np.asarray(require(ax_doc, "edges", expect=list), dtype=np.float64)
        if edges.ndim != 1 or len(edges) < 2 or not np.all(np.diff(edges) > 0):
            raise ConfigError(f"{path}: axis {name!r} edges are not strictly increasing: {edges}")
        axes.append(Grid1D(name=name, edges=edges))
    return Grid2D(
        path=path,
        axes=(axes[0], axes[1]),
        source=str(require(doc, "provenance/source")),
        n_events=int(require(doc, "provenance/n_events")),
    )


@dataclass(frozen=True)
class Binning:
    """The 4D binning: Grid A (Q2 x xB) crossed with Grid B (z x pT2).

    This class owns the index formulas. Nothing else may reimplement them.
    """

    grid_a: Grid2D
    grid_b: Grid2D

    @property
    def n_q2(self) -> int:
        """int : Number of Q2 bins."""
        return self.grid_a.axes[0].n

    @property
    def n_xb(self) -> int:
        """int : Number of xB bins."""
        return self.grid_a.axes[1].n

    @property
    def n_z(self) -> int:
        """int : Number of z bins."""
        return self.grid_b.axes[0].n

    @property
    def n_pt2(self) -> int:
        """int : Number of pT2 bins."""
        return self.grid_b.axes[1].n

    @property
    def n_a(self) -> int:
        """int : Number of Grid A cells."""
        return self.n_q2 * self.n_xb

    @property
    def n_b(self) -> int:
        """int : Number of Grid B cells."""
        return self.n_z * self.n_pt2

    @property
    def n_4d(self) -> int:
        """int : Number of 4D bins."""
        return self.n_a * self.n_b

    @property
    def n_3d(self) -> int:
        """int : Number of 3D (Q2, xB, z) bins used for pT broadening."""
        return self.n_a * self.n_z

    def encode_4d(self, i_q2: int, i_xb: int, i_z: int, i_pt2: int) -> int:
        """Encode a 4D bin index.

        ``bin4d = (i_q2 * n_xb + i_xb) * (n_z * n_pt2) + i_z * n_pt2 + i_pt2``
        """
        return (i_q2 * self.n_xb + i_xb) * self.n_b + i_z * self.n_pt2 + i_pt2

    def decode_4d(self, bin4d: int) -> tuple[int, int, int, int]:
        """Decode a 4D bin index into ``(i_q2, i_xb, i_z, i_pt2)``."""
        a, b = divmod(int(bin4d), self.n_b)
        i_q2, i_xb = divmod(a, self.n_xb)
        i_z, i_pt2 = divmod(b, self.n_pt2)
        return i_q2, i_xb, i_z, i_pt2

    def cell_a_of_4d(self, bin4d: int | NDArray[np.int64]) -> NDArray[np.int64] | int:
        """Return the Grid A cell of a 4D bin (works elementwise on arrays).

        All Grid B cells of one Grid A cell share the same ``N_DIS``; this is
        the map that expresses that.
        """
        return np.asarray(bin4d) // self.n_b if np.ndim(bin4d) else int(bin4d) // self.n_b

    def encode_3d(self, i_q2: int, i_xb: int, i_z: int) -> int:
        """Encode a 3D bin index: ``bin3d = (i_q2 * n_xb + i_xb) * n_z + i_z``."""
        return (i_q2 * self.n_xb + i_xb) * self.n_z + i_z

    def decode_3d(self, bin3d: int) -> tuple[int, int, int]:
        """Decode a 3D bin index into ``(i_q2, i_xb, i_z)``."""
        a, i_z = divmod(int(bin3d), self.n_z)
        i_q2, i_xb = divmod(a, self.n_xb)
        return i_q2, i_xb, i_z

    def bin3d_of_4d(self, bin4d: int) -> int:
        """Return the 3D bin containing a given 4D bin (drops ``i_pt2``)."""
        i_q2, i_xb, i_z, _ = self.decode_4d(bin4d)
        return self.encode_3d(i_q2, i_xb, i_z)

    def centres_4d(self) -> dict[str, NDArray[np.float64]]:
        """Geometric box centres of every 4D bin, keyed by axis name.

        .. warning::
           These are **wrong positions to report a measurement at** and are
           computed here for exactly one purpose: the QA plot that shows how
           far they sit from the count-weighted means (``sec:binning-caveat``).
           In the superseded production 52.6% of bins sat in a top box, 15.5%
           implied ``y > 0.85`` -- violating the DIS cut their own events
           passed -- and 6.9% implied ``nu > E_beam``.
        """
        q2c, xbc = self.grid_a.axes[0].centres, self.grid_a.axes[1].centres
        zc, pt2c = self.grid_b.axes[0].centres, self.grid_b.axes[1].centres
        idx = np.arange(self.n_4d)
        a, b = np.divmod(idx, self.n_b)
        i_q2, i_xb = np.divmod(a, self.n_xb)
        i_z, i_pt2 = np.divmod(b, self.n_pt2)
        return {"q2": q2c[i_q2], "xb": xbc[i_xb], "z": zc[i_z], "pt2": pt2c[i_pt2]}

    def widths_4d(self) -> dict[str, NDArray[np.float64]]:
        """Full box width of every 4D bin, keyed by axis name.

        Plotted as the x error bar so that a wide outer bin can never be
        mistaken for a point measurement.
        """
        q2w, xbw = self.grid_a.axes[0].widths, self.grid_a.axes[1].widths
        zw, pt2w = self.grid_b.axes[0].widths, self.grid_b.axes[1].widths
        idx = np.arange(self.n_4d)
        a, b = np.divmod(idx, self.n_b)
        i_q2, i_xb = np.divmod(a, self.n_xb)
        i_z, i_pt2 = np.divmod(b, self.n_pt2)
        return {"q2": q2w[i_q2], "xb": xbw[i_xb], "z": zw[i_z], "pt2": pt2w[i_pt2]}


@dataclass(frozen=True)
class MggAxis:
    """The gamma-gamma invariant-mass axis, from ``cuts.json /mgg_histogram``."""

    lo: float
    hi: float
    n: int

    @property
    def width(self) -> float:
        """float : Bin width in GeV."""
        return (self.hi - self.lo) / self.n

    @property
    def centres(self) -> NDArray[np.float64]:
        """numpy.ndarray : ``(imgg + 0.5) * width + lo`` for every bin."""
        return self.lo + (np.arange(self.n) + 0.5) * self.width

    @property
    def edges(self) -> NDArray[np.float64]:
        """numpy.ndarray : The ``n + 1`` bin edges."""
        return np.linspace(self.lo, self.hi, self.n + 1)

    def mask(self, lo: float, hi: float) -> NDArray[np.bool_]:
        """Boolean mask of bins whose *centre* lies in ``[lo, hi)``."""
        c = self.centres
        return (c >= lo) & (c < hi)


@dataclass(frozen=True)
class PhiAxis:
    """The ``phi_h`` axis, from ``cuts.json /bsa``. Degrees, Trento convention."""

    lo: float
    hi: float
    n: int

    @property
    def centres_deg(self) -> NDArray[np.float64]:
        """numpy.ndarray : Bin centres in degrees."""
        w = (self.hi - self.lo) / self.n
        return self.lo + (np.arange(self.n) + 0.5) * w

    @property
    def centres_rad(self) -> NDArray[np.float64]:
        """numpy.ndarray : Bin centres in radians."""
        return np.deg2rad(self.centres_deg)


@dataclass(frozen=True)
class ExtractionCuts:
    """The ``extraction`` block of ``cuts.json``. Every field is required."""

    sideband_min: float
    sideband_max: float
    fit_range_min: float
    fit_range_max: float
    n_sigma_window: float
    core_k: float
    n_iter: int
    convergence_tol_frac: float
    amp_min: float
    amp_max: float
    mu_min: float
    mu_max: float
    sigma_min: float
    sigma_max: float
    mu_seed: float
    sigma_seed: float
    min_stats: float
    min_peak_amp: float


@dataclass(frozen=True)
class BsaCuts:
    """The ``bsa`` block of ``cuts.json``.

    ``polarization`` is deliberately allowed to be ``None`` here so that
    :mod:`pi0.bsa` can raise its own explanatory error. It is **never**
    defaulted to a number -- see the note's ``sec:bsa-polarization``.
    """

    phi: PhiAxis
    amp_min: float
    amp_max: float
    min_counts_per_phi_bin: float
    min_phi_bins: int
    min_counts_per_bin: float
    amp_at_bound_tol_frac: float
    sigma_amp_min: float
    max_chi2_ndf: float
    polarization: float | None
    polarization_err: float | None


@dataclass(frozen=True)
class BroadeningCuts:
    """The ``broadening`` block of ``cuts.json``."""

    min_counts: float


@dataclass(frozen=True)
class Config:
    """Everything the extraction stage needs to know, loaded from disk.

    Attributes
    ----------
    config_dir : pathlib.Path
        Directory the configuration was loaded from.
    cuts_path : pathlib.Path
        Path to ``cuts.json``.
    cuts_sha256 : str
        SHA-256 of ``cuts.json`` as loaded. Compared against the hash Stage B
        recorded in its provenance -- see :func:`pi0.io.check_provenance`.
    """

    config_dir: Path
    cuts_path: Path
    cuts_sha256: str
    raw: dict[str, Any]
    binning: Binning
    mgg: MggAxis
    extraction: ExtractionCuts
    bsa: BsaCuts
    broadening: BroadeningCuts


def load_config(config_dir: str | Path) -> Config:
    """Load ``cuts.json`` and the two grid JSONs, validating them against each other.

    Parameters
    ----------
    config_dir : str or pathlib.Path
        The repository's ``config/`` directory.

    Returns
    -------
    Config
        The fully validated configuration.

    Raises
    ------
    ConfigError
        On any missing key, non-monotonic edge array, or disagreement between
        ``cuts.json``'s declared grid shape and the grid files' actual shape.
        Nothing is defaulted and nothing is repaired.
    """
    config_dir = Path(config_dir)
    cuts_path = config_dir / "cuts.json"
    if not cuts_path.is_file():
        raise ConfigError(f"no cuts.json at {cuts_path}")
    cuts = json.loads(cuts_path.read_text())

    grid_a = _load_grid(config_dir / "binning" / "grid_A_q2_xb.json", ("q2", "xb"))
    grid_b = _load_grid(config_dir / "binning" / "grid_B_z_pt2.json", ("z", "pt2"))
    binning = Binning(grid_a=grid_a, grid_b=grid_b)

    # cuts.json declares the grid shape independently of the grid files. If the
    # two ever disagree, one of them is stale -- refuse rather than pick one.
    for key, got in (
        ("binning/grid_a/n_q2", binning.n_q2),
        ("binning/grid_a/n_xb", binning.n_xb),
        ("binning/grid_b/n_z", binning.n_z),
        ("binning/grid_b/n_pt2", binning.n_pt2),
    ):
        want = int(require(cuts, key))
        if want != got:
            raise ConfigError(
                f"cuts.json {key} = {want} but the grid file gives {got}. "
                f"The configuration is internally inconsistent; fix it, do not guess."
            )

    mgg = MggAxis(
        lo=float(require(cuts, "mgg_histogram/min_gev")),
        hi=float(require(cuts, "mgg_histogram/max_gev")),
        n=int(require(cuts, "mgg_histogram/bins")),
    )

    amp_max_raw = require(cuts, "extraction/fit_bounds").get("amp_max", "__missing__")
    if amp_max_raw == "__missing__":
        raise ConfigError("required config key 'extraction/fit_bounds/amp_max' is missing")
    extraction = ExtractionCuts(
        sideband_min=float(require(cuts, "extraction/sideband_min_gev")),
        sideband_max=float(require(cuts, "extraction/sideband_max_gev")),
        fit_range_min=float(require(cuts, "extraction/fit_range_min_gev")),
        fit_range_max=float(require(cuts, "extraction/fit_range_max_gev")),
        n_sigma_window=float(require(cuts, "extraction/n_sigma_window")),
        core_k=float(require(cuts, "extraction/core_k")),
        n_iter=int(require(cuts, "extraction/n_iter")),
        convergence_tol_frac=float(require(cuts, "extraction/convergence_tol_frac")),
        amp_min=float(require(cuts, "extraction/fit_bounds/amp_min")),
        # null amp_max means "unbounded above" and is the one documented null in
        # this block -- see the _amp_max_comment in cuts.json.
        amp_max=np.inf if amp_max_raw is None else float(amp_max_raw),
        mu_min=float(require(cuts, "extraction/fit_bounds/mu_min_gev")),
        mu_max=float(require(cuts, "extraction/fit_bounds/mu_max_gev")),
        sigma_min=float(require(cuts, "extraction/fit_bounds/sigma_min_gev")),
        sigma_max=float(require(cuts, "extraction/fit_bounds/sigma_max_gev")),
        mu_seed=float(require(cuts, "extraction/fit_seed/mu_seed_gev")),
        sigma_seed=float(require(cuts, "extraction/fit_seed/sigma_seed_gev")),
        min_stats=float(require(cuts, "extraction/min_stats")),
        min_peak_amp=float(require(cuts, "extraction/min_peak_amp")),
    )
    amp_seed_rule = str(require(cuts, "extraction/fit_seed/amp_seed"))
    if amp_seed_rule != "max_of_subtracted_spectrum":
        raise ConfigError(
            f"extraction/fit_seed/amp_seed = {amp_seed_rule!r} is not a rule this code "
            f"implements (only 'max_of_subtracted_spectrum')."
        )

    bsa_doc = require(cuts, "bsa")
    bsa = BsaCuts(
        phi=PhiAxis(
            lo=float(require(cuts, "bsa/phi_min_deg")),
            hi=float(require(cuts, "bsa/phi_max_deg")),
            n=int(require(cuts, "bsa/n_phi_bins")),
        ),
        amp_min=float(require(cuts, "bsa/fit_bounds/amp_min")),
        amp_max=float(require(cuts, "bsa/fit_bounds/amp_max")),
        min_counts_per_phi_bin=float(require(cuts, "bsa/quality/min_counts_per_phi_bin")),
        min_phi_bins=int(require(cuts, "bsa/quality/min_phi_bins")),
        min_counts_per_bin=float(require(cuts, "bsa/quality/min_counts_per_bin")),
        amp_at_bound_tol_frac=float(require(cuts, "bsa/quality/amp_at_bound_tol_frac")),
        sigma_amp_min=float(require(cuts, "bsa/quality/sigma_amp_min")),
        max_chi2_ndf=float(require(cuts, "bsa/quality/max_chi2_ndf")),
        # These two are REQUIRED KEYS but MAY BE NULL. Null means "not measured",
        # which pi0.bsa turns into a refusal to run, not into 0.85.
        polarization=_optional_float(bsa_doc, "polarization/value"),
        polarization_err=_optional_float(bsa_doc, "polarization/error"),
    )

    broadening = BroadeningCuts(min_counts=float(require(cuts, "broadening/min_counts")))

    return Config(
        config_dir=config_dir,
        cuts_path=cuts_path,
        cuts_sha256=_sha256(cuts_path),
        raw=cuts,
        binning=binning,
        mgg=mgg,
        extraction=extraction,
        bsa=bsa,
        broadening=broadening,
    )


def _optional_float(node: Any, path: str) -> float | None:
    """Fetch a key that must EXIST but is allowed to be explicitly null.

    The distinction matters: a missing key is a broken config, whereas an
    explicit ``null`` is a deliberate statement that the quantity has not been
    measured. Only the beam polarization uses this.
    """
    parts = path.replace(".", "/").split("/")
    cur = node
    for i, key in enumerate(parts):
        if not isinstance(cur, dict) or key not in cur:
            raise ConfigError(
                f"required config key {path!r} is missing (it may be null, but it must be present)"
            )
        cur = cur[key]
    return None if cur is None else float(cur)
