"""Run the pi0 skim across many HIPO files locally, N processes at a time.

Usage::

    python -m pi0.local_batch --farm config/farm/LD2.farm.json --concurrent 8
    python -m pi0.local_batch --input /path/to/hipo/dir --target LD2 --concurrent 8
    python -m pi0.local_batch --farm config/farm/LD2.farm.json --concurrent 8 --stage-b

This is the farm's local sibling: same selection, same configs, same refusals,
no scheduler. Use it to prove a configuration works on a handful of files before
handing ~28,000 jobs to SWIF2.

Why this parallelises over FILES, not record ranges
---------------------------------------------------
clas-framework's ``local_batch`` slices ONE big HIPO file into ``--slices``
record ranges and hands each to a subprocess, because its analysis binary takes
one input and honours ``--record-range``. Neither is true here: ``stageA_skim``
takes one HIPO file per invocation *by contract* (it exits 5 on a multi-run
file, because its provenance header records a single run and a single model, and
recording one for a file holding several would be a lie), and it has no
``--record-range`` at all. RG-D production is ~30,866 files for LD2 alone, so
the file *is* the natural unit and there is nothing to gain from sub-slicing it.

Consequently ``--slices`` does not exist here; ``--concurrent`` is the only
parallelism knob, and it means what it says.

Two deliberate improvements over the tool this mirrors
------------------------------------------------------
* **A failure tells you which file, which exit code, and which log.**
  clas-framework's local_batch prints ``Slices: 100 (99 OK, 1 failed)`` and
  stops, even though it holds the index, exit code and duration of every slice.
  Finding the failure means grepping 100 logs by hand.
* **Ctrl-C kills the children.** clas-framework's has no signal handling: the
  children die only because the terminal happens to deliver SIGINT to the whole
  foreground process group, which stops being true the moment anything puts them
  in their own group. Here they are terminated explicitly.
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from .progress import Progress
from .farm import (
    FarmConfigError,
    InputFile,
    Rejected,
    TARGET_TO_RUNLIST,
    discover_inputs,
    load_farm_config,
    load_runs,
    run_of,
    runs_for,
)

__all__ = ["main"]

#: stageA_skim's exit codes, from src/stageA_skim/main.cxx. A number alone is a
#: mystery at 2am; this is the difference between "3" and "of course, RG-D".
EXIT_MEANING: dict[int, str] = {
    1: "usage error",
    2: "no run number found in the first events (RUN::config never filled)",
    3: "NO GBT PHOTON MODEL COVERS THIS RUN -- expected on RG-D unless "
       "photon.allow_rga_fallback is true",
    4: "general error (see the log)",
    5: "the file holds more than one run (stageA_skim takes one run per file)",
}


@dataclass
class Job:
    index: int
    src: InputFile
    out: Path
    log: Path


@dataclass
class Result:
    job: Job
    rc: int
    seconds: float


_ABORT = False


def _install_sigint(procs: dict[int, subprocess.Popen]) -> None:
    def handler(signum, frame):  # noqa: ARG001
        global _ABORT
        if _ABORT:
            return
        _ABORT = True
        print("\n\ninterrupted -- terminating running skims ...", file=sys.stderr)
        for p in list(procs.values()):
            try:
                p.terminate()
            except (OSError, ProcessLookupError):
                pass
    signal.signal(signal.SIGINT, handler)


def _run_one(job: Job, exe: str, cuts: Path, target: str, max_events: int | None,
             procs: dict[int, subprocess.Popen]) -> Result:
    if _ABORT:
        return Result(job, 130, 0.0)
    job.out.parent.mkdir(parents=True, exist_ok=True)
    job.log.parent.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "--input", job.src.path, "--output", str(job.out),
           "--config", str(cuts), "--target", target]
    if max_events:
        cmd += ["--max-events", str(max_events)]
    t0 = time.monotonic()
    try:
        with job.log.open("wb") as lf:
            p = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT)
            procs[job.index] = p
            rc = p.wait()
            procs.pop(job.index, None)
    except OSError as e:
        # A spawn failure is one job's problem, not the run's. Record it like any
        # other failure so the summary can name the file, rather than letting a
        # traceback out of the pool and losing every other result.
        job.log.parent.mkdir(parents=True, exist_ok=True)
        job.log.write_text(f"could not launch {' '.join(cmd)}\n{e}\n")
        return Result(job, 127, time.monotonic() - t0)
    return Result(job, rc, time.monotonic() - t0)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="pi0.local_batch",
        description="Run the pi0 Stage A skim across many HIPO files locally, N at a time.",
    )
    src = p.add_argument_group("input (one of --farm or --input)")
    src.add_argument("--farm", type=Path, default=None,
                     help="a farm JSON. Gives inputs, target and the run-list filter -- the same "
                          "selection the farm would make.")
    src.add_argument("--input", type=Path, default=None,
                     help="a directory of HIPO files, or one file. Needs --target.")
    src.add_argument("--target", default=None, choices=sorted(TARGET_TO_RUNLIST),
                     help="required with --input")
    src.add_argument("--polarity", default="outbending", choices=("outbending", "inbending"))
    src.add_argument("--no-run-filter", action="store_true",
                     help="skip the run-list filter. ONLY for a file that is not RG-D production "
                          "(a test file, a single run you are debugging). A production run without "
                          "the filter silently mixes torus polarities.")

    p.add_argument("--config", type=Path, default=Path("config"), help="config DIRECTORY (default: config)")
    p.add_argument("--exe", default="./build/src/stageA_skim/stageA_skim")
    p.add_argument("--stageb-exe", default="./build/src/stageB_bin/stageB_bin")
    p.add_argument("--outdir", type=Path, default=Path("slim"), help="where slims land (default: slim/)")
    p.add_argument("--concurrent", type=int, default=max(1, (os.cpu_count() or 4) // 2),
                   help="how many skims run at once (default: half your cores)")
    p.add_argument("--max-events", type=int, default=None,
                   help="passed to each skim. A PREFIX of each file, not a sample -- stamped into "
                        "every output's provenance. For smoke tests only.")
    p.add_argument("--max-files", type=int, default=None, help="process at most N files (for testing)")
    p.add_argument("--dry-run", action="store_true", help="plan and print; run nothing")
    p.add_argument("--stage-b", action="store_true",
                   help="after Stage A, run the single Stage B over every slim produced")
    args = p.parse_args(argv)

    if bool(args.farm) == bool(args.input):
        print("error: pass exactly one of --farm or --input.", file=sys.stderr)
        return 1
    cuts = args.config / "cuts.json"
    if not cuts.is_file():
        print(f"error: no cuts.json at {cuts}. --config takes the config DIRECTORY.", file=sys.stderr)
        return 1

    # Check the executable BEFORE planning. clas-framework's batch downgrades a
    # missing executable to a warning ("assuming it exists on the worker node"),
    # which is defensible for a farm job and useless here: locally there is no
    # other node, so a bad --exe means every job fails identically. The default
    # is repo-root-relative, and running from python/ is the obvious way to get
    # this wrong.
    exe = Path(args.exe)
    if not (exe.is_file() and os.access(exe, os.X_OK)):
        print(f"error: no executable at {args.exe}\n"
              f"       cwd is {Path.cwd()}. The default --exe is relative to the REPOSITORY ROOT;\n"
              f"       run from there, or pass an absolute --exe.", file=sys.stderr)
        return 1

    # ---- resolve the file list -----------------------------------------
    try:
        if args.farm:
            cfg = load_farm_config(args.farm)
            runs = load_runs(args.config)
            accepted, rejected = discover_inputs(cfg, runs)
            target = cfg.target
        else:
            if not args.target:
                print("error: --input needs --target.", file=sys.stderr)
                return 1
            target = args.target
            root = args.input
            paths = (sorted(str(x) for x in root.rglob("*.hipo")) if root.is_dir()
                     else [str(root)] if root.is_file() else [])
            if not paths:
                print(f"error: no .hipo files at {root}", file=sys.stderr)
                return 1
            accepted, rejected = [], []
            wanted = None if args.no_run_filter else runs_for(load_runs(args.config), args.polarity, target)
            for q in paths:
                r = run_of(q)
                if r is None and not args.no_run_filter:
                    rejected.append(Rejected(q, "no run number in the path"))
                    continue
                if wanted is not None and r not in wanted:
                    rejected.append(Rejected(q, f"run {r} not in the {args.polarity} "
                                                f"{TARGET_TO_RUNLIST[target]} run list"))
                    continue
                accepted.append(InputFile(path=q, run=r if r is not None else 0))
            accepted.sort(key=lambda f: (f.run, f.path))
    except (FarmConfigError, OSError, json.JSONDecodeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    if rejected:
        print(f"rejected {len(rejected)} file(s) before starting:")
        shown: dict[str, int] = {}
        for r in rejected:
            shown[r.reason] = shown.get(r.reason, 0) + 1
        for reason, n in sorted(shown.items(), key=lambda kv: -kv[1])[:6]:
            print(f"  {n:5d}  {reason}")
    if not accepted:
        print("error: no input files survived the run-list filter. Nothing to do.\n"
              "       If this is not RG-D production data, pass --no-run-filter.", file=sys.stderr)
        return 1
    if args.max_files:
        accepted = accepted[: args.max_files]

    jobs = [
        Job(index=i, src=f,
            out=args.outdir / f"{f.run:06d}" / f"slim_{Path(f.path).stem}.root",
            log=args.outdir / "logs" / f"{Path(f.path).stem}.log")
        for i, f in enumerate(accepted)
    ]

    print(f"target      : {target} ({args.polarity})")
    print(f"files       : {len(jobs)}")
    print(f"runs        : {len({j.src.run for j in jobs})}")
    print(f"concurrent  : {args.concurrent}")
    print(f"slims       : {args.outdir}/<run>/slim_<stem>.root")
    print(f"logs        : {args.outdir}/logs/<stem>.log")
    if args.max_events:
        print(f"max-events  : {args.max_events} PER FILE -- every output is a truncated prefix")

    if args.dry_run:
        print("\nDRY RUN. First 5 commands:")
        for j in jobs[:5]:
            print(f"  {args.exe} --input {j.src.path} --output {j.out} "
                  f"--config {cuts} --target {target}")
        if len(jobs) > 5:
            print(f"  ... and {len(jobs) - 5} more")
        return 0

    # ---- run ------------------------------------------------------------
    procs: dict[int, subprocess.Popen] = {}
    _install_sigint(procs)
    results: list[Result] = []
    t0 = time.monotonic()
    # A live bar on a terminal, a periodic line in a batch log. The old code
    # wrote a raw "\r ... " every file, which is right on a TTY and unreadable in
    # a captured log -- exactly the split Progress exists to handle.
    bar = Progress(f"skim {target}", len(jobs))
    with ThreadPoolExecutor(max_workers=args.concurrent) as ex:
        futs = {ex.submit(_run_one, j, args.exe, cuts, target, args.max_events, procs): j for j in jobs}
        for fut in as_completed(futs):
            results.append(fut.result())
            bar.add()
    bar.finish()

    results.sort(key=lambda r: r.job.index)
    ok = [r for r in results if r.rc == 0]
    bad = [r for r in results if r.rc != 0]
    elapsed = time.monotonic() - t0

    print(f"\n{'=' * 76}")
    print(f"Stage A: {len(jobs)} file(s) -- {len(ok)} ok, {len(bad)} failed, {elapsed:.1f}s")
    print("=" * 76)

    if bad:
        # The thing clas-framework's local_batch does not tell you.
        by_rc: dict[int, list[Result]] = {}
        for r in bad:
            by_rc.setdefault(r.rc, []).append(r)
        print("\nFAILURES, by exit code:\n")
        for rc in sorted(by_rc):
            rs = by_rc[rc]
            print(f"  exit {rc} -- {EXIT_MEANING.get(rc, 'unknown')}   ({len(rs)} file(s))")
            for r in rs[:5]:
                print(f"      {r.job.src.path}")
                print(f"        log: {r.job.log}")
            if len(rs) > 5:
                print(f"      ... and {len(rs) - 5} more with exit {rc}")
            print()
        if 3 in by_rc:
            print("  exit 3 is stageA_skim REFUSING to score RG-D photons with an RG-A inbending")
            print("  model. To take that fallback anyway, set photon.allow_rga_fallback = true in")
            print("  config/cuts.json -- every output is then stamped gbt.fallback_used = TRUE and")
            print("  the Python stage refuses to publish from it without --allow-unpublishable.\n")

    if _ABORT:
        print("interrupted; Stage B skipped.")
        return 130
    if not ok:
        print("no slims produced.")
        return 1

    if not args.stage_b:
        print(f"{len(ok)} slim(s) written. Next: make_grid over them, then --stage-b.")
        return 0 if not bad else 1

    if bad:
        print("refusing to run Stage B: some Stage A jobs failed, so the slim set is incomplete.")
        print("N_DIS -- the R_A normalisation denominator -- would be short by those files, with")
        print("nothing in the output recording it. Fix the failures, or re-run without --stage-b")
        print("and chain the slims you accept by hand.")
        return 1

    # ---- Stage B --------------------------------------------------------
    slims = [str(r.job.out) for r in ok]
    out = args.outdir / f"binned_{target}.root"
    cmd = [args.stageb_exe]
    for s in slims:
        cmd += ["--input", s]
    cmd += ["--output", str(out), "--config", str(cuts),
            "--grid-a", str(args.config / "binning" / "grid_A_q2_xb.json"),
            "--grid-b", str(args.config / "binning" / "grid_B_z_pt2.json")]
    print(f"\nStage B over {len(slims)} slim(s) -> {out}")
    print("One pool across every file, single-threaded by design.\n")
    return subprocess.run(cmd).returncode


if __name__ == "__main__":
    raise SystemExit(main())
