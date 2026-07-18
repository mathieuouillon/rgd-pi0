"""Reading Stage B's output, and refusing to use it when it says not to.

Stage B writes four flat TTrees plus two blocks of ``TNamed`` provenance
(its own under ``/provenance``, Stage A's propagated verbatim under
``/provenance_stageA``). Everything here is uproot-only: no ROOT, no
dictionary.

The provenance is not decoration. :func:`check_provenance` reads it and turns
it into blockers -- a file skimmed with an RGA-trained GBT fallback, or binned
on placeholder grids, cannot produce a quotable physics number, and this module
will not let one out silently. That refusal is the entire reason Stage B writes
provenance at all.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import uproot
from numpy.typing import NDArray

from .config import Binning, Config

__all__ = [
    "ProvenanceError",
    "SchemaError",
    "Provenance",
    "StageBData",
    "Blocker",
    "check_provenance",
    "load",
]


class ProvenanceError(RuntimeError):
    """The file's provenance says its contents are not fit to quote."""


class SchemaError(RuntimeError):
    """The file's trees do not match the configured binning."""


@dataclass(frozen=True)
class Provenance:
    """One block of ``TNamed`` provenance, read as name -> title.

    Attributes
    ----------
    directory : str
        The ROOT directory the block was read from.
    entries : dict of str to str
        ``fName -> fTitle`` for every ``TNamed`` in the directory.
    """

    directory: str
    entries: dict[str, str]

    def get(self, key: str, default: str | None = None) -> str | None:
        """Return one provenance value, or ``default`` if absent."""
        return self.entries.get(key, default)

    def require(self, key: str) -> str:
        """Return one provenance value, raising if it is absent.

        Raises
        ------
        ProvenanceError
            If ``key`` is not in the block. A provenance key that the analysis
            depends on and that the file does not carry is a hard error: the
            file was written by a Stage B that did not make the promise this
            code is checking.
        """
        if key not in self.entries:
            raise ProvenanceError(
                f"{self.directory}/{key} is absent. This file was not written by a Stage B "
                f"that records it, so the check that depends on it cannot be performed."
            )
        return self.entries[key]

    def starts_true(self, key: str) -> bool:
        """True if a provenance value begins with ``TRUE``/``true``/``yes``.

        Stage B writes booleans as prose (``"TRUE -- photons scored by a model
        trained on OTHER data"``), so the leading token is the machine-readable
        part and the rest is for the human reading the file.
        """
        val = self.require(key).strip().lower()
        return val.startswith(("true", "yes"))

    def summary(self, keys: tuple[str, ...] = ()) -> str:
        """A compact one-line stamp for plots and CSV headers."""
        picked = keys or tuple(self.entries)
        return "; ".join(f"{k}={self.entries[k]}" for k in picked if k in self.entries)


@dataclass
class Blocker:
    """One reason the file's contents must not be quoted as physics.

    Attributes
    ----------
    key : str
        The provenance key (or config comparison) that raised it.
    message : str
        Human-readable explanation, quoting the file.
    fatal : bool
        ``True`` if this alone disqualifies every physics number in the file.
        ``False`` marks a serious caveat that does not by itself invalidate the
        extraction machinery.
    """

    key: str
    message: str
    fatal: bool = True

    def __str__(self) -> str:
        return f"[{'BLOCKER' if self.fatal else 'WARNING'}] {self.key}: {self.message}"


@dataclass
class StageBData:
    """Stage B's four trees, densified onto the configured binning.

    Attributes
    ----------
    path : pathlib.Path
        The file read.
    binning : Binning
        The binning the arrays are indexed on.
    n_same, n_mixed : numpy.ndarray
        Shape ``(n_4d, n_mgg)``, integer. Same- and mixed-event gamma-gamma
        pair counts.
    sum_q2, sum_xb, sum_z, sum_pt2 : numpy.ndarray
        Shape ``(n_4d, n_mgg)``, float. **Same-event only** kinematic sums,
        weighted by ``n_same`` in the same cell. No mixed pair contributes to
        any of these -- a mixed pair is not a pi0 and its kinematics describe
        nothing (``provenance/abscissa.mixed_excluded``).
    ptb_counts, ptb_sum_pt2, ptb_sum_pt4 : numpy.ndarray
        Shape ``(n_3d, n_mgg)``. Same-event only, as above.
    n_dis : numpy.ndarray
        Shape ``(n_a,)``, integer. Inclusive DIS events per Grid A cell, filled
        once per event. The normalisation denominator.
    bsa : numpy.ndarray
        Shape ``(n_4d, n_mgg, n_phi, 2)``, integer, with the last axis
        ``[helicity +1, helicity -1]``. Densified from the sparse tree **via
        its index columns** -- a positional read of that tree is silently wrong.
    prov, prov_stage_a : Provenance
        Stage B's own provenance and Stage A's, propagated verbatim.
    blockers : list of Blocker
        What :func:`check_provenance` found.
    """

    path: Path
    binning: Binning
    n_same: NDArray[np.int64]
    n_mixed: NDArray[np.int64]
    sum_q2: NDArray[np.float64]
    sum_xb: NDArray[np.float64]
    sum_z: NDArray[np.float64]
    sum_pt2: NDArray[np.float64]
    ptb_counts: NDArray[np.int64]
    ptb_sum_pt2: NDArray[np.float64]
    ptb_sum_pt4: NDArray[np.float64]
    n_dis: NDArray[np.int64]
    bsa: NDArray[np.int64]
    prov: Provenance
    prov_stage_a: Provenance
    blockers: list[Blocker] = field(default_factory=list)

    @property
    def target(self) -> str:
        """str : The target this file was skimmed for (from Stage A's provenance)."""
        return self.prov_stage_a.require("target").strip()

    @property
    def sums(self) -> dict[str, NDArray[np.float64]]:
        """dict : The four same-event kinematic sums, keyed by axis name."""
        return {"q2": self.sum_q2, "xb": self.sum_xb, "z": self.sum_z, "pt2": self.sum_pt2}

    def stamp(self) -> str:
        """A provenance stamp to print on every plot and CSV this data produces."""
        a, b = self.prov_stage_a, self.prov
        bits = [
            f"target={a.get('target', '?')}",
            f"run={a.get('run', '?')}",
            f"pol={a.get('polarity', '?')}",
            f"events={b.get('events.read', '?')}",
            f"binning_hash={b.get('binning.provenance_hash', '?')}",
            f"stageA={a.get('stageA_skim.created_utc', '?')}",
            f"stageB={b.get('stageB_bin.created_utc', '?')}",
        ]
        if self.blockers:
            n_fatal = sum(1 for x in self.blockers if x.fatal)
            bits.append(f"NOT PUBLISHABLE ({n_fatal} blockers)")
        return " | ".join(bits)

    def mixed_3d(self) -> NDArray[np.int64]:
        """Project the 4D mixed spectrum onto the 3D (Q2, xB, z) bins.

        Stage B writes no mixed-event spectrum for ``ptb3d``, but ``ptb3d`` is
        exactly the ``i_pt2`` projection of ``spectra`` (verified: its
        ``counts`` equal ``n_same`` summed over ``i_pt2``, and its ``sum_pt2``
        likewise). So the mixed background for a 3D bin is the 4D mixed
        spectrum summed over ``i_pt2``, and the moments can be sideband-
        subtracted like the yield.

        Returns
        -------
        numpy.ndarray
            Shape ``(n_3d, n_mgg)``.
        """
        b = self.binning
        shaped = self.n_mixed.reshape(b.n_a, b.n_z, b.n_pt2, -1)
        return shaped.sum(axis=2).reshape(b.n_3d, -1)


def _read_dir_entries(f: uproot.ReadOnlyDirectory, name: str) -> dict[str, str]:
    """The ``fName -> fTitle`` map of every ``TNamed`` in directory ``name``."""
    entries: dict[str, str] = {}
    for key in f[name].keys(recursive=False):
        obj = f[name][key]
        members = getattr(obj, "all_members", {})
        if "fName" in members and "fTitle" in members:
            entries[str(members["fName"])] = str(members["fTitle"])
    return entries


def _read_provenance(f: uproot.ReadOnlyDirectory, name: str) -> Provenance:
    """Read a directory of ``TNamed`` provenance into name -> title."""
    if name not in {k.rstrip(";1234567890").split(";")[0] for k in f.keys(recursive=False)}:
        raise ProvenanceError(
            f"{f.file_path} has no /{name} block. Every file this stage reads must carry "
            f"its provenance; one that does not cannot be checked and is not usable."
        )
    return Provenance(directory=name, entries=_read_dir_entries(f, name))


#: Stage B writes one ``provenance_stageA_NNN`` directory per input file, in
#: canonical order, VERBATIM and UNMERGED (stageB_bin main.cxx: "ONE
#: /provenance_stageA_NNN PER INPUT ... they are not merged either"). There is
#: no single ``provenance_stageA`` -- reading one is what :func:`_read_stage_a`
#: fixes.
_STAGE_A_RE = re.compile(r"provenance_stageA_\d+")


def _tripped(block: Provenance, key: str) -> bool:
    """``starts_true`` without raising on an absent key (a missing flag is not set)."""
    val = block.get(key)
    return bool(val) and val.strip().lower().startswith(("true", "yes"))


def _read_stage_a(f: uproot.ReadOnlyDirectory) -> Provenance:
    """Aggregate the per-input ``provenance_stageA_NNN`` blocks into one.

    Stage B keeps them unmerged, one per input, so a file holding several runs
    cannot pass one run's provenance off as the whole. For the publish guard,
    though, the file is only as clean as its dirtiest input: a fallback or a
    truncation in ANY input taints every number. So ``gbt.fallback_used`` and
    the ``events.max_events_requested`` truncation marker are carried over from
    any input that trips them; the rest is taken from the first input as
    representative (all inputs of one Stage B run share a skim config).
    """
    names = sorted(
        n for n in {k.split(";")[0] for k in f.keys(recursive=False)} if _STAGE_A_RE.fullmatch(n)
    )
    if not names:
        raise ProvenanceError(
            f"{f.file_path} has no /provenance_stageA_NNN block. Every file this stage reads "
            f"must carry its Stage A provenance (one per input); one that does not cannot be "
            f"checked and is not usable."
        )
    blocks = [Provenance(directory=n, entries=_read_dir_entries(f, n)) for n in names]

    merged: dict[str, str] = dict(blocks[0].entries)
    fallback = next((b for b in blocks if _tripped(b, "gbt.fallback_used")), None)
    if fallback is not None:
        merged["gbt.fallback_used"] = fallback.require("gbt.fallback_used")
        model = fallback.get("gbt.model")
        if model:
            merged["gbt.model"] = model
    truncated = next(
        (b for b in blocks if "TRUNCATED" in (b.get("events.max_events_requested") or "").upper()),
        None,
    )
    if truncated is not None:
        merged["events.max_events_requested"] = truncated.require("events.max_events_requested")

    return Provenance(directory=f"provenance_stageA_[{len(names)} input(s)]", entries=merged)


def check_provenance(prov: Provenance, prov_a: Provenance, cfg: Config) -> list[Blocker]:
    """Decide whether this file's numbers may be quoted as physics.

    Parameters
    ----------
    prov : Provenance
        Stage B's provenance block.
    prov_a : Provenance
        Stage A's propagated provenance block.
    cfg : Config
        The configuration this stage loaded, for the hash comparison.

    Returns
    -------
    list of Blocker
        Empty if the file is clean. Callers pass this to :func:`load` via
        ``allow_unpublishable`` to decide between refusing and proceeding.

    Notes
    -----
    The checks, and why each one exists:

    * ``gbt.fallback_used`` -- the photons were scored by a GBT model trained
      on *other* data (RGA rather than RG-D). The photon sample is then not the
      one the cuts describe.
    * placeholder grids -- the bin edges are hand-authored product edges, not
      equal-statistics edges fitted to data. Stage B says it directly: "No
      result may be quoted from a placeholder grid."
    * ``config.sha256_matches_stageA`` -- the skim and the binning ran against
      different ``cuts.json`` files, so the electron and photon selections are
      not this run's.
    * config hash drift -- the ``cuts.json`` on disk now is not the one Stage B
      recorded. Not fatal (a comment edit changes the hash), but it means the
      cut values in play here are not provably the ones that made the file.
    * truncated run -- ``--max-events`` was given, so no yield or normalisation
      in the file is complete.
    """
    blockers: list[Blocker] = []

    if prov_a.starts_true("gbt.fallback_used"):
        blockers.append(
            Blocker(
                "provenance_stageA/gbt.fallback_used",
                f"the photon sample was built with a FALLBACK GBT model "
                f"(gbt.model={prov_a.get('gbt.model', '?')!r}): "
                f"{prov_a.require('gbt.fallback_used')}. "
                f"cuts.json sets photon.allow_rga_fallback = "
                f"{cfg.raw.get('photon', {}).get('allow_rga_fallback')!r}. "
                f"The photon selection is not the one the cuts describe.",
            )
        )

    for label, grid in (("grid_a", cfg.binning.grid_a), ("grid_b", cfg.binning.grid_b)):
        if grid.is_placeholder:
            blockers.append(
                Blocker(
                    f"{label}.provenance.source",
                    f"{grid.path.name} carries source={grid.source!r} with n_events="
                    f"{grid.n_events}: these are hand-authored product edges, NOT "
                    f"equal-statistics edges fitted to data. Stage B's own provenance says: "
                    f"'No result may be quoted from a placeholder grid.'",
                )
            )

    matches = prov.get("config.sha256_matches_stageA", "")
    if matches and matches.strip().lower().startswith(("false", "no")):
        blockers.append(Blocker("provenance/config.sha256_matches_stageA", matches))

    file_cfg_hash = prov.get("config.sha256", "")
    if file_cfg_hash and file_cfg_hash.strip() != cfg.cuts_sha256:
        blockers.append(
            Blocker(
                "config.sha256",
                f"the cuts.json loaded here ({cfg.cuts_path}, sha256 "
                f"{cfg.cuts_sha256[:16]}...) is NOT the one Stage B ran against "
                f"({file_cfg_hash[:16]}...). The cut values applied to this data are not "
                f"provably the ones in play now.",
                fatal=False,
            )
        )

    max_ev = prov_a.get("events.max_events_requested", "")
    if "TRUNCATED" in max_ev.upper():
        blockers.append(Blocker("provenance_stageA/events.max_events_requested", max_ev))

    return blockers


def load(
    path: str | Path,
    cfg: Config,
    *,
    allow_unpublishable: bool = False,
    stream: object = None,
) -> StageBData:
    """Read a Stage B file, checking its provenance before returning it.

    Parameters
    ----------
    path : str or pathlib.Path
        The Stage B ROOT file.
    cfg : Config
        Configuration, which fixes the binning the trees are densified onto.
    allow_unpublishable : bool, default False
        If ``False`` and the provenance carries any fatal blocker, raise
        :class:`ProvenanceError`. If ``True``, print every blocker loudly to
        stderr and proceed anyway -- for QA and debugging, never for a number
        that leaves the building.
    stream : file-like, optional
        Where to print warnings. Defaults to ``sys.stderr``.

    Returns
    -------
    StageBData
        The densified trees plus both provenance blocks.

    Raises
    ------
    ProvenanceError
        If the file is not fit to quote and ``allow_unpublishable`` is False.
    SchemaError
        If a tree's length or index column disagrees with the configured
        binning.
    """
    stream = stream if stream is not None else sys.stderr
    path = Path(path)
    b = cfg.binning
    n_mgg = cfg.mgg.n

    with uproot.open(path) as f:
        prov = _read_provenance(f, "provenance")
        prov_a = _read_stage_a(f)
        blockers = check_provenance(prov, prov_a, cfg)

        fatal = [x for x in blockers if x.fatal]
        if blockers:
            print("\n" + "=" * 78, file=stream)
            print(f"PROVENANCE CHECK FAILED for {path}", file=stream)
            print("=" * 78, file=stream)
            for x in blockers:
                print(f"  {x}\n", file=stream)
            if fatal and not allow_unpublishable:
                raise ProvenanceError(
                    f"{path} carries {len(fatal)} fatal provenance blocker(s); refusing to "
                    f"produce physics from it. Pass allow_unpublishable=True (or "
                    f"--allow-unpublishable) to proceed anyway for QA -- the results are then "
                    f"DIAGNOSTIC ONLY and must not be quoted."
                )
            if fatal:
                print(
                    "PROCEEDING ANYWAY because allow_unpublishable=True.\n"
                    "EVERY NUMBER BELOW IS DIAGNOSTIC ONLY AND MUST NOT BE QUOTED.\n"
                    + "=" * 78 + "\n",
                    file=stream,
                )

        spectra = f["spectra"].arrays(library="np")
        ptb = f["ptb3d"].arrays(library="np")
        ndis_t = f["n_dis"].arrays(library="np")
        bsa_t = f["bsa"].arrays(library="np")

    # --- spectra: dense, n_4d * n_mgg rows in index order. Verify, do not assume.
    want = b.n_4d * n_mgg
    if len(spectra["bin4d"]) != want:
        raise SchemaError(
            f"{path}: spectra has {len(spectra['bin4d'])} rows, but the configured binning "
            f"({b.n_4d} 4D bins x {n_mgg} mgg bins) implies {want}. The file and the config "
            f"describe different binnings."
        )
    if not (
        np.array_equal(spectra["bin4d"], np.repeat(np.arange(b.n_4d), n_mgg))
        and np.array_equal(spectra["imgg"], np.tile(np.arange(n_mgg), b.n_4d))
    ):
        raise SchemaError(
            f"{path}: spectra's index columns are not in dense (bin4d, imgg) order. "
            f"Stage B promises they are; refusing to reshape positionally."
        )

    def _sq(name: str, dtype: type) -> NDArray:
        return spectra[name].astype(dtype).reshape(b.n_4d, n_mgg)

    want3 = b.n_3d * n_mgg
    if len(ptb["bin3d"]) != want3:
        raise SchemaError(
            f"{path}: ptb3d has {len(ptb['bin3d'])} rows, expected {want3} "
            f"({b.n_3d} 3D bins x {n_mgg} mgg bins)."
        )
    if not (
        np.array_equal(ptb["bin3d"], np.repeat(np.arange(b.n_3d), n_mgg))
        and np.array_equal(ptb["imgg"], np.tile(np.arange(n_mgg), b.n_3d))
    ):
        raise SchemaError(f"{path}: ptb3d's index columns are not in dense order.")

    # --- n_dis: dense over Grid A cells, but scatter by cell_a rather than trust order.
    n_dis = np.zeros(b.n_a, dtype=np.int64)
    cell_a = ndis_t["cell_a"].astype(np.int64)
    if np.any((cell_a < 0) | (cell_a >= b.n_a)):
        raise SchemaError(f"{path}: n_dis has cell_a outside [0, {b.n_a}).")
    np.add.at(n_dis, cell_a, ndis_t["n_dis"].astype(np.int64))

    # --- bsa: SPARSE. Decode via index columns; a positional read is silently wrong.
    bsa = np.zeros((b.n_4d, n_mgg, cfg.bsa.phi.n, 2), dtype=np.int64)
    if len(bsa_t["counts"]):
        i4 = bsa_t["bin4d"].astype(np.int64)
        im = bsa_t["imgg"].astype(np.int64)
        ip = bsa_t["iphi"].astype(np.int64)
        hel = bsa_t["helicity"].astype(np.int64)
        bad = np.unique(hel[(hel != 1) & (hel != -1)])
        if len(bad):
            raise SchemaError(
                f"{path}: bsa carries helicity values {bad.tolist()}; only +1/-1 are allowed "
                f"(helicity == 0 is undefined and must be dropped by Stage B)."
            )
        for arr, hi, nm in ((i4, b.n_4d, "bin4d"), (im, n_mgg, "imgg"), (ip, cfg.bsa.phi.n, "iphi")):
            if np.any((arr < 0) | (arr >= hi)):
                raise SchemaError(f"{path}: bsa has {nm} outside [0, {hi}).")
        # ihel = 0 for helicity +1, 1 for helicity -1 (provenance/binning.index_formula_bsa)
        ihel = np.where(hel == 1, 0, 1)
        np.add.at(bsa, (i4, im, ip, ihel), bsa_t["counts"].astype(np.int64))

    data = StageBData(
        path=path,
        binning=b,
        n_same=_sq("n_same", np.int64),
        n_mixed=_sq("n_mixed", np.int64),
        sum_q2=_sq("sum_q2", np.float64),
        sum_xb=_sq("sum_xb", np.float64),
        sum_z=_sq("sum_z", np.float64),
        sum_pt2=_sq("sum_pt2", np.float64),
        ptb_counts=ptb["counts"].astype(np.int64).reshape(b.n_3d, n_mgg),
        ptb_sum_pt2=ptb["sum_pt2"].astype(np.float64).reshape(b.n_3d, n_mgg),
        ptb_sum_pt4=ptb["sum_pt4"].astype(np.float64).reshape(b.n_3d, n_mgg),
        n_dis=n_dis,
        bsa=bsa,
        prov=prov,
        prov_stage_a=prov_a,
        blockers=blockers,
    )

    # ptb3d must be the i_pt2 projection of spectra. broadening() relies on it to
    # borrow the mixed spectrum; if it is ever untrue, that borrowing is invalid.
    proj = data.n_same.reshape(b.n_a, b.n_z, b.n_pt2, n_mgg).sum(axis=2).reshape(b.n_3d, n_mgg)
    if not np.array_equal(proj, data.ptb_counts):
        raise SchemaError(
            f"{path}: ptb3d/counts is NOT the i_pt2 projection of spectra/n_same. The 3D and 4D "
            f"trees describe different samples, so the 3D bins cannot borrow the 4D mixed "
            f"spectrum and the pT2 moments cannot be sideband-subtracted."
        )

    return data
