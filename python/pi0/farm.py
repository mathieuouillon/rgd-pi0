"""Farm plumbing: run lists, input discovery, chunking, and SWIF2 synthesis.

This module is the library half of the batch tooling. :mod:`pi0.batch` (JLab
farm, SWIF2) and :mod:`pi0.local_batch` (one machine, N subprocesses) are the
two CLIs on top of it; everything they share lives here.

It is deliberately dependency-free -- no uproot, no numpy, no ROOT. It runs on
an ifarm login node against the system python with nothing installed, because
that is where you submit from.

What this module refuses to do
------------------------------
The whole point of putting a layer between "a /mss directory" and "28,000 jobs"
is that the layer can say no. It refuses to:

* build a production out of files whose run is not in the requested
  ``(polarity, target)`` run list -- see :func:`discover_inputs`. Inbending and
  outbending runs share one directory tree, so a plain recursive scan silently
  mixes torus polarities and *nothing downstream would notice*: the polarity is
  not in the slim schema.
* submit Stage A when ``photon.allow_rga_fallback`` is false, because every one
  of those jobs would exit 3 on RG-D (no GBT photon model covers runs
  18305-19131). Finding that out from 28,000 identical failures is a bad evening.
* submit Stage A when it *is* true without saying so, loudly, in the summary.
* accept a malformed ``time``. clas-framework's ``format_swif2_time`` returns
  malformed input unchanged and lets SWIF2 complain at submit time; this
  refuses at parse time instead. Same conversion, different failure mode.

Fidelity to clas-framework
--------------------------
The emitted SWIF2 is modelled on ``clas-framework``'s ``tools/batch.cpp`` +
``framework/include/clas12/framework/core/BatchSystem.hpp``, not on its
``docs/mkdocs/user-guide/batch-system.md``, which is stale: the doc shows a
``-shell '<cmd>'`` argument, but the code emits the command as a **trailing
positional** after the option flags, and clas-framework's own
``test/test_batch_synthesis.cpp`` asserts ``-shell`` is absent. Ported from the
code.

Divergences from clas-framework, and why
----------------------------------------
* **Config is JSON, not TOML** -- this project's config is JSON throughout.
* **``time`` refuses when malformed** rather than passing through (above).
* **Stage ordering is by separate workflow, not by antecedent.** clas-framework
  has no job-dependency support at all (its ``JobSpec`` has four fields and no
  ``-antecedent`` is emitted anywhere). This pipeline needs Stage B to follow
  every Stage A. Expressing that as antecedents would mean one Stage B job with
  ~28,000 antecedents; instead ``--stage a`` and ``--stage b`` are two workflows
  and :mod:`pi0.batch` refuses to build B until A's slims are on disk. That is
  the honest shape: B genuinely cannot start until *all* of A is done, so there
  is nothing for a scheduler to overlap.
* **Globs in ``inputs`` work here; they do not in clas-framework.** It has a
  ``FileUtils::is_glob`` its batch tool never calls, so ``inputs =
  ".../*.hipo"`` there falls through and is pushed verbatim as a fake path --
  no error, and the job gets a nonexistent ``-input``. Supported and tested
  here rather than reproduced.
"""

from __future__ import annotations

import glob
import hashlib
import json
import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

__all__ = [
    "FarmConfigError",
    "RunListError",
    "FarmConfig",
    "InputFile",
    "Rejected",
    "Chunk",
    "TARGET_TO_RUNLIST",
    "load_runs",
    "load_farm_config",
    "discover_inputs",
    "chunk_inputs",
    "format_swif2_time",
    "swif2_script",
    "stagea_wrapper",
    "job_inputs",
    "must_stage",
]


class FarmConfigError(RuntimeError):
    """The farm configuration is malformed, inconsistent, or unsafe to submit."""


class RunListError(RuntimeError):
    """An input file's run cannot be reconciled with the requested run list."""


# ---------------------------------------------------------------------------
# Targets vs run lists
# ---------------------------------------------------------------------------

#: There are FOUR analysis targets but THREE run lists. Cu and Sn are two foils
#: in one assembly exposed in the same CuSn runs, separated offline by the
#: electron vz window (note sec:vertex) -- so ``stageA_skim --target Cu`` and
#: ``--target Sn`` read the *same* input files and differ only in that window.
#: This mapping is the only place that fact is encoded.
TARGET_TO_RUNLIST: dict[str, str] = {
    "LD2": "LD2",
    "CxC": "CxC",
    "Cu": "CuSn",
    "Sn": "CuSn",
}

#: The run number in a CLAS12 HIPO filename. CLAS12 has more than one naming
#: convention and the production uses both:
#:
#:   dst/recon/018434/rec_clas_018434.evio.00000-00004.hipo   per-run DSTs
#:   dst/train/SIDIS/SIDIS_018431.hipo                        train skims
#:
#: so this matches ``_<run>.`` rather than a fixed prefix. Applied to the
#: BASENAME only -- a directory called ``run_12345.old`` must not be mistaken for
#: the run of a file inside it.
_RUN_IN_NAME = re.compile(r"_(\d{4,6})\.")

#: Fallback: a six-digit directory, as the per-run DST tree uses. Train skims
#: live in a directory named for the train (``SIDIS``), not the run, so they never
#: reach this.
_RUN_IN_DIR = re.compile(r"/(\d{6})/")


def _require(d: dict[str, Any], key: str, where: str) -> Any:
    """Fetch ``key`` or raise. Never defaults -- see :mod:`pi0.config`."""
    if key not in d:
        raise FarmConfigError(
            f"{where}: missing required key {key!r}. This file is the only place that value may "
            f"live; there is no default and inventing one here is how a production acquires a "
            f"parameter nobody chose."
        )
    return d[key]


# ---------------------------------------------------------------------------
# Run lists
# ---------------------------------------------------------------------------


def load_runs(config_dir: Path) -> dict[str, Any]:
    """Load ``config/runs.json``.

    The RG-D run lists, transcribed from clas-framework's ``clas12/Runs.hpp``.
    The note (sec:runs) says of the originals: *"These lists are not applied by
    any binary"*. Loading them here is what changes that.
    """
    path = Path(config_dir) / "runs.json"
    if not path.is_file():
        raise FarmConfigError(
            f"no run list at {path}. --config takes the config DIRECTORY (which must contain "
            f"runs.json and cuts.json), not a single file."
        )
    with path.open() as f:
        return json.load(f)


def runs_for(runs: dict[str, Any], polarity: str, target: str) -> set[int]:
    """The set of run numbers for one ``(polarity, target)``.

    Maps the target through :data:`TARGET_TO_RUNLIST` first, so Cu and Sn both
    resolve to the CuSn list.
    """
    if target not in TARGET_TO_RUNLIST:
        raise FarmConfigError(
            f"unknown target {target!r}. stageA_skim takes one of: "
            f"{', '.join(sorted(TARGET_TO_RUNLIST))}."
        )
    if polarity not in ("outbending", "inbending"):
        raise FarmConfigError(f"unknown polarity {polarity!r}; want 'outbending' or 'inbending'.")
    runlist = TARGET_TO_RUNLIST[target]
    if polarity not in runs:
        raise FarmConfigError(f"runs.json has no {polarity!r} block.")
    if runlist not in runs[polarity]:
        raise FarmConfigError(f"runs.json {polarity!r} has no {runlist!r} list.")
    return set(int(r) for r in runs[polarity][runlist])


def run_of(path: str | Path) -> int | None:
    """The run number a DST path claims, or None.

    Reads the filename first -- ``rec_clas_018434.evio...`` and
    ``SIDIS_018431.hipo`` are both production naming -- and falls back to a
    six-digit parent directory. Returns None rather than guessing: a caller that
    cannot identify a file's run must refuse it, not assume, because the run is
    what decides the target and the polarity.
    """
    m = _RUN_IN_NAME.search(Path(path).name)
    if m:
        return int(m.group(1))
    m = _RUN_IN_DIR.search(str(path))
    if m:
        return int(m.group(1))
    return None


# ---------------------------------------------------------------------------
# Farm configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Swif2Spec:
    workflow: str
    account: str
    partition: str
    cores: int
    ram_gb: int
    disk_gb: int
    time: str
    constraint: str = ""
    modules: tuple[str, ...] = ()
    tags: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class FarmConfig:
    target: str
    polarity: str
    inputs: tuple[str, ...]
    output_dir: str
    log_dir: str
    files_per_job: int
    exclude_runs: frozenset[int]
    swif2: Swif2Spec
    env: dict[str, str]
    source: Path


def load_farm_config(path: str | Path) -> FarmConfig:
    """Load and validate a farm JSON.

    Mirrors clas-framework's farm TOML schema, in JSON, minus the keys this
    pipeline resolves elsewhere. Note two nesting details ported from
    clas-framework's *code* (its struct layout suggests otherwise, and copying
    the struct gets them backwards):

    * ``env`` is a child of ``farm``, not of ``farm.swif2``.
    * ``modules`` is an object with a ``load`` array, not a bare array.
    """
    path = Path(path)
    if not path.is_file():
        raise FarmConfigError(f"no farm config at {path}")
    with path.open() as f:
        doc = json.load(f)

    farm = _require(doc, "farm", str(path))
    where = f"{path}:/farm"

    target = _require(farm, "target", where)
    if target not in TARGET_TO_RUNLIST:
        raise FarmConfigError(
            f"{where}/target: {target!r} is not a stageA_skim target "
            f"({', '.join(sorted(TARGET_TO_RUNLIST))})."
        )
    polarity = farm.get("polarity", "outbending")

    inputs = _require(farm, "inputs", where)
    if isinstance(inputs, str):
        inputs = [inputs]
    if not isinstance(inputs, list) or not inputs or not all(isinstance(i, str) for i in inputs):
        raise FarmConfigError(f"{where}/inputs: want a non-empty string or array of strings.")

    fpj = int(farm.get("files_per_job", 5))
    if fpj < 1:
        raise FarmConfigError(f"{where}/files_per_job: must be >= 1, got {fpj}.")

    sw = _require(farm, "swif2", where)
    swhere = f"{where}/swif2"
    modules_block = sw.get("modules", {})
    if isinstance(modules_block, list):
        raise FarmConfigError(
            f"{swhere}/modules: want an object with a 'load' array, not a bare array. "
            f'Write: "modules": {{"load": ["gcc/13.2.0", "root/6.36.04"]}}. '
            f"A bare array parses fine and silently loads nothing, and a job whose env lacks ROOT "
            f"fails on the node rather than at submit."
        )
    modules = tuple(modules_block.get("load", []))

    swif2 = Swif2Spec(
        workflow=_require(sw, "workflow", swhere),
        account=sw.get("account", "clas12"),
        partition=sw.get("partition", "production"),
        cores=int(sw.get("cores", 1)),
        ram_gb=int(sw.get("ram_gb", 4)),
        disk_gb=int(sw.get("disk_gb", 20)),
        time=str(_require(sw, "time", swhere)),
        constraint=str(sw.get("constraint", "")),
        modules=modules,
        tags={str(k): str(v) for k, v in sw.get("tags", {}).items()},
    )
    format_swif2_time(swif2.time, where=f"{swhere}/time")  # validate now, not on the node

    env = farm.get("env", {})
    if not isinstance(env, dict):
        raise FarmConfigError(f"{where}/env: want an object of NAME -> value.")

    return FarmConfig(
        target=target,
        polarity=polarity,
        inputs=tuple(inputs),
        output_dir=str(_require(farm, "output_dir", where)),
        log_dir=str(_require(farm, "log_dir", where)),
        files_per_job=fpj,
        exclude_runs=frozenset(int(r) for r in farm.get("exclude_runs", [])),
        swif2=swif2,
        env={str(k): str(v) for k, v in env.items()},
        source=path,
    )


def format_swif2_time(t: str, where: str = "time") -> str:
    """``"04:00:00"`` -> ``"14400s"``, the form SWIF2's ``-time`` wants.

    Same conversion as clas-framework's ``format_swif2_time``
    (BatchSystem.hpp), with one deliberate difference: that one returns
    malformed input unchanged ("fall through; let SWIF2 complain"), so a typo
    surfaces as a scheduler error after you have already submitted. This raises
    instead. ``HH:MM`` is read as ``HH:MM:00``, matching clas-framework.

    A bare ``<int>s`` is passed through, so an already-converted value is
    idempotent.
    """
    t = t.strip()
    if re.fullmatch(r"\d+s", t):
        return t
    m = re.fullmatch(r"(\d+):(\d{1,2})(?::(\d{1,2}))?", t)
    if not m:
        raise FarmConfigError(
            f"{where}: {t!r} is not a time. Want HH:MM:SS (e.g. \"04:00:00\"), HH:MM, or "
            f"<seconds>s. clas-framework would pass this through and let SWIF2 reject it after "
            f"submission; refusing here instead."
        )
    h, mi, s = int(m.group(1)), int(m.group(2)), int(m.group(3) or 0)
    if mi > 59 or s > 59:
        raise FarmConfigError(f"{where}: {t!r} has a minutes/seconds field above 59.")
    total = h * 3600 + mi * 60 + s
    if total <= 0:
        raise FarmConfigError(f"{where}: {t!r} is zero. A job with no time limit is not a job.")
    return f"{total}s"


# ---------------------------------------------------------------------------
# Input discovery
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class InputFile:
    path: str
    run: int


@dataclass(frozen=True)
class Rejected:
    path: str
    reason: str


def _expand_one(spec: str) -> list[str]:
    """A directory (scanned recursively for ``*.hipo``), a glob, or a literal file.

    Globs go through :mod:`glob`, not hand-rolled path arithmetic: a relative
    pattern like ``*.hipo`` has no directory part to split on, and getting that
    wrong silently returns nothing rather than erroring.
    """
    p = Path(spec)
    if p.is_dir():
        return sorted(str(m) for m in p.rglob("*.hipo") if m.is_file())
    if p.is_file():
        return [str(p)]
    if any(ch in spec for ch in "*?["):
        return sorted(m for m in glob.glob(spec, recursive=True) if Path(m).is_file())
    return []


def discover_inputs(
    cfg: FarmConfig,
    runs: dict[str, Any],
    *,
    listing: Sequence[str] | None = None,
) -> tuple[list[InputFile], list[Rejected]]:
    """Resolve ``cfg.inputs`` to files, then filter by the run list.

    Returns ``(accepted, rejected)``. Rejection is never silent: every dropped
    file comes back with a reason, and :mod:`pi0.batch` prints the tally.

    ``listing`` lets a caller supply the file list directly (a ``.dat`` file, or
    a test fixture) instead of touching the filesystem. That matters on a login
    node: clas-framework's ``recursive_directory_iterator`` **throws** on a
    missing directory and the throw is not caught, so a bad /mss path kills the
    tool with an uncaught exception. Here a missing path is a reported reason.

    THE FILTER IS THE POINT. ``/mss/.../LD2/dst/recon/`` holds LD2 *inbending*
    runs 18305-18336 alongside the outbending production. A scan without this
    filter mixes torus polarities into one dataset, and the slim schema does not
    record polarity, so nothing downstream can notice.
    """
    if listing is not None:
        raw = list(listing)
        missing_specs: list[Rejected] = []
    else:
        raw = []
        missing_specs = []
        for spec in cfg.inputs:
            found = _expand_one(spec)
            if not found:
                missing_specs.append(
                    Rejected(spec, "resolved to no files (missing path, empty directory, or a glob that matched nothing)")
                )
            raw.extend(found)

    wanted = runs_for(runs, cfg.polarity, cfg.target)
    runlist = TARGET_TO_RUNLIST[cfg.target]

    accepted: list[InputFile] = []
    rejected: list[Rejected] = list(missing_specs)
    seen: set[str] = set()

    for path in raw:
        if path in seen:
            rejected.append(Rejected(path, "duplicate: listed by more than one inputs entry"))
            continue
        seen.add(path)
        run = run_of(path)
        if run is None:
            rejected.append(
                Rejected(path, "no run number in the filename or a six-digit parent directory")
            )
            continue
        if run in cfg.exclude_runs:
            rejected.append(Rejected(path, f"run {run} is in exclude_runs"))
            continue
        if run not in wanted:
            rejected.append(
                Rejected(path, f"run {run} is not in the {cfg.polarity} {runlist} run list")
            )
            continue
        accepted.append(InputFile(path=path, run=run))

    # Canonical order. The pool downstream is order-dependent (reservoir
    # sampling is a function of the sequence of offers), so the file order must
    # not be whatever the filesystem happened to yield.
    accepted.sort(key=lambda f: (f.run, f.path))
    return accepted, rejected


# ---------------------------------------------------------------------------
# Chunking
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Chunk:
    index: int
    files: tuple[InputFile, ...]

    @property
    def name_suffix(self) -> str:
        return f"job_{self.index:05d}"


def chunk_inputs(files: Sequence[InputFile], files_per_job: int) -> list[Chunk]:
    """Slice the accepted files into per-job chunks.

    Chunks never span runs. stageA_skim takes ONE HIPO file per invocation and
    exits 5 on a multi-run file, so a job runs the skim once per file in its
    chunk; keeping a chunk within one run means its outputs share a run
    directory and a failed job retries a coherent unit.
    """
    if files_per_job < 1:
        raise FarmConfigError(f"files_per_job must be >= 1, got {files_per_job}")
    by_run: dict[int, list[InputFile]] = {}
    for f in files:
        by_run.setdefault(f.run, []).append(f)

    chunks: list[Chunk] = []
    for run in sorted(by_run):
        group = by_run[run]
        for i in range(0, len(group), files_per_job):
            chunks.append(Chunk(index=len(chunks), files=tuple(group[i : i + files_per_job])))
    return chunks


# ---------------------------------------------------------------------------
# SWIF2 synthesis
# ---------------------------------------------------------------------------


def _q(s: str) -> str:
    return shlex.quote(str(s))


#: Only /mss is tape. Everything else in the farm's namespace is a mounted
#: cluster filesystem that a worker node reads directly, so SWIF2 must NOT be
#: asked to stage it -- see :func:`must_stage`.
_TAPE = ("/mss/", "/w/mss/")


def must_stage(path: str) -> bool:
    """Does SWIF2 have to copy this input to the node, or can the node read it?

    ``-input`` COPIES. For a /mss path that is the whole point: the file is on
    tape and must be brought to disk. For /cache, /work and /volatile it is
    waste, and for RG-D it is fatal waste -- a train skim is ~200 GB
    (SIDIS_018431.hipo, measured), so staging one puts a 202 GB copy into a job
    that asked for 20 GB of disk and reads 2000 events. Observed: both jobs of
    the first submission sat in `preparing` (= staging) with
    site_job_disk_bytes = 20000000000 and local = input_0.hipo.

    A node-visible input is passed to the skim as its own absolute path instead,
    read straight off Lustre. No copy, no disk request, no staging wait.
    """
    return any(str(path).startswith(t) for t in _TAPE)


def job_inputs(chunk: "Chunk") -> list[str]:
    """What the skim opens, per file of a chunk, in job order.

    THE ONE PLACE that decides staged-name vs absolute-path. swif2_script emits
    ``-input`` for exactly the entries this returns as staged names, and
    stagea_wrapper opens exactly what this returns -- they cannot disagree
    because they call the same function.
    """
    return [f"input_{j}.hipo" if must_stage(f.path) else f.path
            for j, f in enumerate(chunk.files)]


def stagea_wrapper(
    cfg: FarmConfig,
    inputs: Sequence[str],
    exe: str,
    max_events: int | None = None,
) -> str:
    """The per-job script SWIF2 runs on the node.

    ``inputs`` is what the skim should open, in job order: either a staged name
    (``input_0.hipo``, relative to the job's scratch dir) for a tape input, or an
    absolute path for one the node reads directly. :func:`must_stage` decides.

    ``max_events`` caps each skim. It is for smoke tests -- an RG-D train skim is
    ~200 GB, so an unbounded job against one is a long, expensive way to discover
    a typo. Every output it produces is a PREFIX of its input and says so in its
    own provenance, and Stage B refuses such a file unless told otherwise.

    stageA_skim takes ONE input per invocation, so this loops. It does NOT use
    ``set -e``: a single bad file must not silently abandon the rest of the
    chunk with a zero exit. Every skim's status is recorded and the worst is
    returned, so SWIF2 marks the job failed and ``swif2 retry-jobs -problems``
    can find it.

    Exit 3 (no GBT photon model covers this run) is called out by name because
    on RG-D it is the expected failure when ``photon.allow_rga_fallback`` is
    false, and a wall of bare "exit 3" is otherwise a mystery.
    """
    lines = [
        "#!/usr/bin/env bash",
        "set -uo pipefail",
        "",
        "# Generated by pi0.batch. Do not edit -- regenerate.",
        "",
        "# The body runs inside run(), whose output is tee'd to log.txt. swif2_script",
        "# declares `-output log.txt`, and -output copies a FILE out of the job dir --",
        "# it does not capture stdout. Without this tee the file never exists, SWIF2",
        "# has nothing to reap, and the job's whole log is lost even though it ran and",
        "# exited 0. Observed on workflow rgd_pi0_stageA_LD2_TRAINTEST.",
        "run() {",
    ]
    for mod in cfg.swif2.modules:
        lines.append(f"module load {_q(mod)}")
    for k, v in sorted(cfg.env.items()):
        lines.append(f"export {k}={_q(v)}")
    # The input list is written out rather than derived from a counter: a job may
    # mix a staged name with an absolute Lustre path, and only the generator
    # knows which is which.
    lines += ["", "inputs=("]
    lines += [f"    {_q(i)}" for i in inputs]
    lines += [
        ")",
        "",
        "rc=0",
        'for i in "${!inputs[@]}"; do',
        '    in="${inputs[$i]}"',
        '    out="slim_${i}.root"',
        '    if [[ ! -f "$in" ]]; then',
        '        echo "INPUT NOT READABLE: $in" >&2',
        '        echo "  A staged name means SWIF2 did not deliver it; an absolute path means the" >&2',
        '        echo "  node cannot see that filesystem." >&2',
        "        rc=4; continue",
        "    fi",
        f"    {_q(exe)} --input \"$in\" --output \"$out\" \\",
        f"        --config cuts.json --target {_q(cfg.target)}"
        + (f" --max-events {int(max_events)}" if max_events else ""),
        "    s=$?",
        '    if [[ $s -ne 0 ]]; then',
        '        echo "stageA_skim FAILED on $in with exit $s" >&2',
        '        if [[ $s -eq 3 ]]; then',
        '            echo "  exit 3 = no GBT photon model covers this run. On RG-D this is EXPECTED" >&2',
        '            echo "  unless photon.allow_rga_fallback is true in cuts.json." >&2',
        "        fi",
        "        [[ $rc -eq 0 ]] && rc=$s",
        "    fi",
        "done",
        "return $rc",
        "}",
        "",
        "# PIPESTATUS[0], not $?: `run | tee` reports tee's status, which is always 0.",
        "run 2>&1 | tee log.txt",
        "exit ${PIPESTATUS[0]}",
    ]
    # Indent the body so run() reads as a function; the header lines up to and
    # including "run() {" and the trailing invocation are left alone.
    out = []
    inside = False
    for ln in lines:
        if ln == "run() {":
            out.append(ln); inside = True; continue
        if ln == "}":
            inside = False; out.append(ln); continue
        out.append(("    " + ln) if (inside and ln) else ln)
    return "\n".join(out) + "\n"


def swif2_script(
    cfg: FarmConfig,
    chunks: Sequence[Chunk],
    *,
    scripts_dir: Path,
    exe: str,
    cuts_path: Path,
) -> str:
    """The ``swif2 create`` / ``add-job`` xN / ``run`` script.

    Modelled on clas-framework's BatchSystem::synthesize: the command is a
    **trailing positional** (``bash wrapper.sh``), not a ``-shell`` argument.
    Tags are emitted in sorted order, matching clas-framework's ``std::map``.

    Every job gets an explicit ``mkdir -p`` of its output directory, because
    SWIF2's ``-output`` does not create the parent and the copy is simply lost
    if it is absent.
    """
    sw = cfg.swif2
    wf = sw.workflow
    out: list[str] = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f"# rgd-pi0 Stage A -- {len(chunks)} jobs, target {cfg.target}, {cfg.polarity}",
        f"# Generated by pi0.batch from {cfg.source}. Regenerate rather than edit.",
        "",
        f"swif2 create -workflow {_q(wf)}",
        "",
    ]
    for ch in chunks:
        job = f"{wf}_{ch.name_suffix}"
        run = ch.files[0].run
        job_out = f"{cfg.output_dir}/{run:06d}/{job}"
        out.append(f"mkdir -p {_q(job_out)} {_q(cfg.log_dir)}")
        parts = [
            f"swif2 add-job -workflow {_q(wf)}",
            f"    -name {_q(job)}",
            f"    -account {_q(sw.account)} -partition {_q(sw.partition)}",
            f"    -cores {sw.cores} -ram {sw.ram_gb}G -disk {sw.disk_gb}G"
            f" -time {_q(format_swif2_time(sw.time))}",
        ]
        if sw.constraint:
            parts.append(f"    -constraint {_q(sw.constraint)}")
        for k in sorted(sw.tags):
            parts.append(f"    -tag {_q(k)} {_q(sw.tags[k])}")
        parts.append(f"    -tag 'target' {_q(cfg.target)} -tag 'run' {_q(str(run))}")
        # Stage ONLY what the node cannot read. -input copies, and for a ~200 GB
        # train skim on /cache that copy is both pointless and larger than the
        # job's whole disk request.
        for j, f in enumerate(ch.files):
            if must_stage(f.path):
                parts.append(f"    -input {_q(f'input_{j}.hipo')} {_q(f.path)}")
        parts.append(f"    -input 'cuts.json' {_q(str(cuts_path.resolve()))}")
        parts.append(f"    -input 'wrapper.sh' {_q(str((scripts_dir / f'{job}.wrapper.sh').resolve()))}")
        for j in range(len(ch.files)):
            parts.append(f"    -output {_q(f'slim_{j}.root')} {_q(f'{job_out}/slim_{j}.root')}")
        parts.append(f"    -output 'log.txt' {_q(f'{cfg.log_dir}/{job}.log')}")
        parts.append("    bash wrapper.sh")
        out.append(" \\\n".join(parts))
        out.append("")

    out.append(f"swif2 run -workflow {_q(wf)}")
    out.append("")
    return "\n".join(out)


def config_sha256(path: Path) -> str:
    """The same digest Stage A stamps, so a drifted config is visible pre-submit."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
