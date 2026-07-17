"""The nuclear multiplicity ratio ``R_A``, per 4D bin.

    R_A = (Y_A / N_DIS_A) / (Y_D / N_DIS_D)

with ``N_DIS`` integrated over the bin's **(Q2, xB) cell only** -- not over
``z`` or ``pT2``. All 25 Grid B cells of one Grid A cell therefore share the
same ``N_DIS``, which is directly verifiable in the output table. LD2 is always
the denominator.

Every bin is reported at its count-weighted, sideband-subtracted abscissa
(:func:`pi0.extract.abscissa`). A bin whose abscissa is not computable is
**dropped and counted**, never reported at a box centre.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from . import extract, io
from .config import Config, load_config
from .extract import BinYield

__all__ = ["ratio_table", "write_csv", "TABLE_DTYPE", "main"]

TABLE_DTYPE = np.dtype(
    [
        ("bin4d", np.int64),
        ("i_q2", np.int64),
        ("i_xb", np.int64),
        ("i_z", np.int64),
        ("i_pt2", np.int64),
        ("cell_a", np.int64),
        ("q2_mean", np.float64),
        ("xb_mean", np.float64),
        ("z_mean", np.float64),
        ("pt2_mean", np.float64),
        ("Y_A", np.float64),
        ("sY_A", np.float64),
        ("Y_D", np.float64),
        ("sY_D", np.float64),
        ("N_DIS_A", np.int64),
        ("N_DIS_D", np.int64),
        ("R", np.float64),
        ("sR", np.float64),
    ]
)


def ratio_table(
    target: list[BinYield],
    deuterium: list[BinYield],
    n_dis_target: NDArray[np.int64],
    n_dis_deuterium: NDArray[np.int64],
    cfg: Config,
    *,
    stream=None,
) -> tuple[NDArray, dict[str, int]]:
    """Build the tidy ``R_A`` table from two sets of per-bin extractions.

    Parameters
    ----------
    target : list of BinYield
        The nuclear target's extractions, in 4D bin order.
    deuterium : list of BinYield
        LD2's extractions, in the same order. Same binning, enforced.
    n_dis_target, n_dis_deuterium : numpy.ndarray
        Inclusive DIS counts per Grid A cell.
    cfg : Config
        Configuration, for the binning.
    stream : file-like, optional
        Where to report dropped bins. Defaults to stderr.

    Returns
    -------
    table : numpy.ndarray
        A structured array with dtype :data:`TABLE_DTYPE`, one row per
        surviving bin.
    dropped : dict of str to int
        Reject code -> how many bins it dropped. Nothing vanishes unaccounted.

    Notes
    -----
    The uncertainty is the note's ``eq:ra-error``::

        sigma_R / R = sqrt( (sY_D/Y_D)^2 + (sY_A/Y_A)^2 + 1/N_DIS_D + 1/N_DIS_A )

    treating the DIS counts as Poisson. **Three omissions are carried
    deliberately**, exactly as the note lists them (``sec:ra-uncertainties``):

    1. The uncertainty on the sideband scale ``alpha`` is not propagated; the
       subtraction treats it as exact.
    2. The +-3 sigma window sum's error ignores the uncertainty on the fitted
       ``mu`` and ``sigma`` that *define* the window.
    3. **R_CxC, R_Cu and R_Sn are correlated** -- they share the LD2 yield and
       the LD2 DIS count in every denominator -- but this formula treats them as
       independent. That is harmless for any single ``R_A`` and *wrong for any
       fit of the A-dependence*, which will understate its own uncertainty.
       Anyone fitting A-dependence off this table must build the covariance
       between targets themselves; the information needed (the shared LD2
       columns ``Y_D``, ``sY_D``, ``N_DIS_D``) is kept in every row for exactly
       that reason.
    """
    stream = stream if stream is not None else sys.stderr
    b = cfg.binning
    if len(target) != b.n_4d or len(deuterium) != b.n_4d:
        raise ValueError(
            f"expected {b.n_4d} bins per target, got {len(target)} and {len(deuterium)}"
        )

    rows: list[tuple] = []
    dropped: dict[str, int] = {}

    def drop(code: str) -> None:
        dropped[code] = dropped.get(code, 0) + 1

    for i in range(b.n_4d):
        ta, de = target[i], deuterium[i]
        if not ta.ok:
            drop(f"target:{ta.reason}")
            continue
        if not de.ok:
            drop(f"ld2:{de.reason}")
            continue
        cell = int(b.cell_a_of_4d(i))
        nda, ndd = int(n_dis_target[cell]), int(n_dis_deuterium[cell])
        if nda <= 0 or ndd <= 0:
            drop("no_dis_normalisation")
            continue
        if ta.value <= 0 or de.value <= 0:
            # A non-positive subtracted yield is a statistical fluctuation, not
            # a measurement; the ratio and its relative error are undefined.
            drop("non_positive_yield")
            continue

        r = (ta.value / nda) / (de.value / ndd)
        rel = np.sqrt(
            (de.error / de.value) ** 2
            + (ta.error / ta.value) ** 2
            + 1.0 / ndd
            + 1.0 / nda
        )
        i_q2, i_xb, i_z, i_pt2 = b.decode_4d(i)
        # The abscissa is the TARGET's. It differs slightly from LD2's, since
        # the two samples populate the box differently -- another reason a shared
        # box centre was never the right label.
        rows.append(
            (
                i, i_q2, i_xb, i_z, i_pt2, cell,
                ta.abscissae["q2"], ta.abscissae["xb"], ta.abscissae["z"], ta.abscissae["pt2"],
                ta.value, ta.error, de.value, de.error, nda, ndd, r, r * rel,
            )
        )

    table = np.array(rows, dtype=TABLE_DTYPE)
    if dropped:
        print(f"ratio: {len(rows)} of {b.n_4d} bins survived. Dropped:", file=stream)
        for code, n in sorted(dropped.items(), key=lambda kv: -kv[1]):
            print(f"    {n:6d}  {code}", file=stream)
    return table, dropped


def write_csv(table: NDArray, path: str | Path, header_lines: list[str]) -> None:
    """Write the table to CSV, stamped with provenance.

    Parameters
    ----------
    table : numpy.ndarray
        A structured array with dtype :data:`TABLE_DTYPE`.
    path : str or pathlib.Path
        Output file.
    header_lines : list of str
        Provenance and caveats, written as leading ``#`` comments. Every result
        file must carry the provenance of the data it came from -- an
        unstamped CSV is how a placeholder-grid number becomes a plot in a talk.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for line in header_lines:
            fh.write(f"# {line}\n")
        fh.write(",".join(table.dtype.names) + "\n")
        for row in table:
            fh.write(",".join(_fmt(v) for v in row) + "\n")


def _fmt(v) -> str:
    """Format one cell: integers plain, floats at repr precision."""
    if isinstance(v, (np.integer, int)):
        return str(int(v))
    return repr(float(v))


def main(argv: list[str] | None = None) -> int:
    """CLI: compute ``R_A`` from a nuclear-target file and an LD2 file."""
    p = argparse.ArgumentParser(description="Nuclear multiplicity ratio R_A per 4D bin.")
    p.add_argument("--target", required=True, type=Path, help="Stage B file for the nuclear target")
    p.add_argument("--ld2", required=True, type=Path, help="Stage B file for LD2 (the denominator)")
    p.add_argument("--config", required=True, type=Path, help="config/ directory")
    p.add_argument("--out", required=True, type=Path, help="output CSV")
    p.add_argument(
        "--allow-unpublishable",
        action="store_true",
        help="proceed despite fatal provenance blockers. Output is DIAGNOSTIC ONLY.",
    )
    args = p.parse_args(argv)

    cfg = load_config(args.config)
    da = io.load(args.target, cfg, allow_unpublishable=args.allow_unpublishable)
    dd = io.load(args.ld2, cfg, allow_unpublishable=args.allow_unpublishable)

    if da.target.upper() == dd.target.upper():
        print(
            f"ERROR: both files carry target={da.target!r}. R_A is a ratio of two DIFFERENT "
            f"targets to LD2; dividing a target by itself measures nothing.",
            file=sys.stderr,
        )
        return 2
    if dd.target.upper() not in {"LD2", "D", "D2"}:
        print(f"ERROR: --ld2 file carries target={dd.target!r}, not LD2.", file=sys.stderr)
        return 2

    ya = extract.extract_all(da, cfg)
    yd = extract.extract_all(dd, cfg)
    table, _ = ratio_table(ya, yd, da.n_dis, dd.n_dis, cfg)

    header = [
        f"R_A = (Y_A/N_DIS_A) / (Y_D/N_DIS_D), per 4D bin.",
        f"target file : {da.path}",
        f"target prov : {da.stamp()}",
        f"LD2 file    : {dd.path}",
        f"LD2 prov    : {dd.stamp()}",
        f"config      : {cfg.cuts_path} sha256={cfg.cuts_sha256}",
        "q2_mean/xb_mean/z_mean/pt2_mean are COUNT-WEIGHTED, SIDEBAND-SUBTRACTED means,",
        "not box centres. Bins whose abscissa was not computable are absent, not defaulted.",
        "sR omits: the alpha uncertainty; the mu/sigma fit error on the window; and the",
        "correlation between targets through the shared LD2 columns (matters for A-dependence).",
    ]
    if da.blockers or dd.blockers:
        header.insert(0, "*** DIAGNOSTIC ONLY -- provenance blockers present, DO NOT QUOTE ***")
    write_csv(table, args.out, header)
    print(f"wrote {len(table)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
