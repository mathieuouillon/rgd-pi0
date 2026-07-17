"""Submit the RG-D pi0 Stage A skim to the JLab farm via SWIF2.

Usage::

    python -m pi0.batch config/farm/LD2.farm.json               # dry-run
    python -m pi0.batch config/farm/LD2.farm.json --submit

Dry-run is the default, as in clas-framework. Unlike clas-framework's, this
dry-run is **side-effect-free apart from the scripts directory**: it writes the
wrapper scripts and the swif2 script (you asked to see them) and nothing else.
clas-framework's C++ tool and its Python skim both perform their full snapshot
copy -- hundreds of MB, clobbering a previous workflow's frozen program dir --
before ever checking ``--submit``.

The pipeline this drives
------------------------
::

    Stage A   one HIPO -> one slim         <- THIS TOOL, fanned out over ~28k jobs
    make_grid once, over the slims         <- one job, needs the slims
    Stage B   all slims -> one binned      <- ONE per target, needs ALL of A
    pi0.*     extraction                   <- local, cheap

Stage B is deliberately not submitted by this tool. It cannot start until every
Stage A job is done, so there is nothing for a scheduler to overlap, and
expressing "wait for all of A" as SWIF2 antecedents would mean one job with ~28,000
of them. Run ``--stage b`` once A is complete; it prints the exact
``stageB_bin`` command with every slim chained.

The pre-flight
--------------
The most valuable thing here happens before a single job is submitted. On RG-D,
``stageA_skim`` exits 3 for every run in 18305-19131 because no GBT photon model
covers them, so a production submitted with ``photon.allow_rga_fallback = false``
is ~28,000 identical failures and a wasted evening. That is checked, named, and
refused up front.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

from .farm import (
    FarmConfigError,
    RunListError,
    job_inputs,
    chunk_inputs,
    config_sha256,
    discover_inputs,
    format_swif2_time,
    load_farm_config,
    load_runs,
    stagea_wrapper,
    swif2_script,
)

__all__ = ["main"]


def _generalise(reason: str) -> str:
    """Collapse per-file reasons so the tally counts kinds, not runs.

    ``run 18306 is not in the outbending LD2 run list`` and the same for 18307
    are one fact reported twice, not two facts.
    """
    return re.sub(r"\brun \d+\b", "run <N>", reason)


def _read_listing(path: Path) -> list[str]:
    with path.open() as f:
        return [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith("#")]


def _preflight(cuts_path: Path, n_jobs: int) -> list[str]:
    """Everything that would make every job fail, checked once instead of N times.

    Returns a list of blocker strings; empty means go.
    """
    blockers: list[str] = []
    if not cuts_path.is_file():
        return [f"no cut config at {cuts_path}"]
    with cuts_path.open() as f:
        cuts = json.load(f)

    photon = cuts.get("photon", {})
    if not photon.get("allow_rga_fallback", False):
        blockers.append(
            "photon.allow_rga_fallback is false, and NO GBT PHOTON MODEL COVERS RG-D.\n"
            "      The model map stops at run 16772; RG-D is 18305-19131. stageA_skim will exit 3\n"
            f"      on every one of these {n_jobs} jobs. This is stageA_skim working as designed --\n"
            "      it refuses rather than scoring RG-D photons with an RG-A inbending model.\n"
            "      To take the fallback anyway, set photon.allow_rga_fallback = true in\n"
            "      config/cuts.json. Every output then carries gbt.fallback_used = TRUE and the\n"
            "      Python stage will refuse to publish from it without --allow-unpublishable.\n"
            "      Do that only for a study you intend to report as such."
        )
    return blockers


#: Filesystems an ifarm worker node mounts. A job's executable must live on one
#: of these: the job's CWD is a scratch directory holding only what SWIF2 staged,
#: and nothing else of yours is there.
#: BOTH the friendly names and the physical mounts behind them, and the physical
#: ones are not optional. On ifarm the friendly names are symlinks:
#:
#:     /work     -> /w/work    (and /work/clas12b -> /w/ceph24/hallb/clas12)
#:     /volatile -> /lustre24/expphy/volatile
#:     /cache    -> /lustre24/expphy/cache
#:     /mss      -> /w/mss
#:
#: and Python's os.getcwd() always reports the PHYSICAL path, so a relative
#: --config given from a /work checkout arrives here as /w/ceph24/... . Listing
#: only the friendly names made this reject the farm's own paths -- and ifarm is
#: the only place the check ever runs, since it gates --submit and swif2 exists
#: nowhere else.
_NODE_VISIBLE = (
    # what a human types
    "/work/", "/volatile/", "/cache/", "/mss/", "/home/", "/group/", "/scigroup/", "/u/",
    # what those resolve to, and what getcwd() reports
    "/w/", "/lustre24/",
)


def _resolve_exe(exe: str) -> tuple[str, list[str], list[str]]:
    """Make the executable path node-usable, or explain why it is not.

    Returns ``(absolute_path, problems, notes)``. ``problems`` block submission;
    ``notes`` are things worth saying that this machine cannot decide.

    THE JOB'S CWD IS NOT YOUR CHECKOUT. SWIF2 runs the wrapper in a fresh scratch
    directory containing only the staged inputs, so a relative ``--exe`` such as
    ``./build/src/stageA_skim/stageA_skim`` resolves to nothing there and every
    job exits 127, "No such file or directory". Found by running a generated
    wrapper in a simulated job directory -- a syntax check cannot see it, and
    neither can a dry run.

    Absolute is necessary but not sufficient: the path must also be on a
    filesystem the nodes mount. A build under /tmp, or on a laptop, is not.

    No libraries are staged, deliberately: meson bakes an rpath into the build
    tree, so the binary runs with no LD_LIBRARY_PATH provided the tree itself is
    node-visible. (clas-framework's swif2_pi0_skim.py instead snapshots the exe
    and every *.so into a frozen directory -- necessary there because it copies
    by basename into a flat dir, which loses the rpath.)
    """
    p = Path(exe)
    problems: list[str] = []
    notes: list[str] = []
    absolute = str(p if p.is_absolute() else (Path.cwd() / p).resolve())

    if not p.is_absolute():
        problems.append(
            f"--exe {exe!r} is a RELATIVE path.\n"
            f"      A SWIF2 job runs in a scratch directory holding only the staged inputs, so a\n"
            f"      relative executable path resolves to nothing there and EVERY job exits 127.\n"
            f"      Resolved against this cwd it would be:\n"
            f"        {absolute}\n"
            f"      Pass that with --exe if it is what you mean, and make sure it is on a\n"
            f"      filesystem the worker nodes mount ({', '.join(v.rstrip('/') for v in _NODE_VISIBLE)})."
        )
    elif not any(absolute.startswith(v) for v in _NODE_VISIBLE):
        problems.append(
            f"--exe {absolute!r} is absolute, but is not on a filesystem the worker nodes\n"
            f"      mount ({', '.join(v.rstrip('/') for v in _NODE_VISIBLE)}). Every job would exit 127.\n"
            f"      Build on /work (or another shared path) and point --exe there."
        )
    else:
        # The path names a node-visible filesystem. Whether it EXISTS is only
        # knowable if that filesystem is mounted here -- on ifarm it is, on a
        # laptop generating a script to run later it is not. Refuse what we know
        # is wrong; say what we cannot check rather than guessing either way.
        root = "/" + Path(absolute).parts[1]
        if not Path(root).is_dir():
            notes.append(
                f"{root} is not mounted on this machine, so whether {absolute} exists could not be\n"
                f"      checked. The path is the right SHAPE for a worker node; that is all this can\n"
                f"      tell you. Generating a script here to submit from ifarm is fine -- run the\n"
                f"      same command there and this check becomes real."
            )
        elif not Path(absolute).is_file():
            problems.append(
                f"no executable at {absolute}.\n"
                f"      It is referenced by path, not copied, so it must exist and stay put for the\n"
                f"      lifetime of the workflow."
            )
    return absolute, problems, notes


def _check_staged_paths(cuts_path: Path, scripts_dir: Path) -> list[str]:
    """Every path SWIF2 stages by reference must be reachable from a worker node.

    ``-input`` copies FROM the path you give it, at job start, on the node. So
    cuts.json and each wrapper.sh are subject to exactly the same rule as the
    executable: a path under /Users or /tmp is not there.

    This is what makes a generated script non-portable. Generating on a laptop
    bakes in laptop paths; the script is then readable but not submittable. On
    ifarm the same command bakes in /work paths and is fine -- so the fix is
    always "generate it where you will submit it", never "edit the paths".
    """
    bad: list[str] = []
    for label, p in (("cuts.json (--config)", cuts_path.resolve()),
                     ("the wrapper scripts (--scripts-dir)", scripts_dir.resolve())):
        if not any(str(p).startswith(v) for v in _NODE_VISIBLE):
            bad.append(f"{label}: {p}")
    return bad


def _grid_is_placeholder(grid_path: Path) -> bool:
    if not grid_path.is_file():
        return True
    try:
        with grid_path.open() as f:
            doc = json.load(f)
    except (OSError, json.JSONDecodeError):
        return True
    blob = json.dumps(doc).lower()
    return "placeholder" in blob


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="pi0.batch",
        description="Submit the RG-D pi0 Stage A skim to the JLab farm via SWIF2.",
    )
    p.add_argument("farm_config", type=Path, help="a farm JSON, e.g. config/farm/LD2.farm.json")
    p.add_argument("--config", type=Path, default=Path("config"), help="config DIRECTORY (default: config)")
    p.add_argument("--submit", action="store_true", help="actually run the generated swif2 script (default: dry-run)")
    p.add_argument("--stage", choices=("a", "b"), default="a", help="which stage to drive (default: a)")
    p.add_argument("--exe", default="./build/src/stageA_skim/stageA_skim", help="path to stageA_skim ON THE NODE")
    p.add_argument("--stageb-exe", default="./build/src/stageB_bin/stageB_bin", help="path to stageB_bin")
    p.add_argument("--scripts-dir", type=Path, default=Path("batch_scripts"),
                   help="where the wrapper + swif2 scripts are written (default: batch_scripts)")
    p.add_argument("--files-per-job", type=int, default=None, help="override farm/files_per_job")
    p.add_argument("--workflow", default=None, help="override farm/swif2/workflow")
    p.add_argument("--partition", default=None,
                   help="override farm/swif2/partition. 'priority' is the fast-turnaround queue -- "
                        "use it for smoke tests; leave production runs on 'production'.")
    p.add_argument("--time", default=None,
                   help="override farm/swif2/time, e.g. 01:00:00. Converted to SWIF2's seconds form; "
                        "a malformed value is refused here rather than by the scheduler after submit.")
    p.add_argument("--max-jobs", type=int, default=None, help="synthesize at most N jobs (for testing)")
    p.add_argument("--max-events", type=int, default=None,
                   help="cap each skim at N events. FOR SMOKE TESTS: an RG-D train skim is ~200 GB, "
                        "so an unbounded job on one is a long, expensive way to find a typo. Every "
                        "output is then a PREFIX of its input, says so in its own provenance, and is "
                        "refused by Stage B unless you pass --allow-truncated-inputs.")
    p.add_argument("--file-list", type=Path, default=None,
                   help="read the input file list from this file instead of scanning. Useful on a login "
                        "node, or to reproduce an old production exactly.")
    p.add_argument("--slim-dir", type=Path, default=None, help="--stage b: where Stage A's slims landed")
    p.add_argument("--allow-rga-fallback-production", action="store_true",
                   help="proceed even though every output will be stamped gbt.fallback_used")
    args = p.parse_args(argv)

    try:
        cfg = load_farm_config(args.farm_config)
        runs = load_runs(args.config)
    except (FarmConfigError, OSError, json.JSONDecodeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    sw_over = {}
    if args.workflow:
        sw_over["workflow"] = args.workflow
    if args.partition:
        sw_over["partition"] = args.partition
    if args.time:
        format_swif2_time(args.time, where="--time")  # refuse now, not on the farm
        sw_over["time"] = args.time
    if sw_over:
        cfg = replace(cfg, swif2=replace(cfg.swif2, **sw_over))
    fpj = args.files_per_job or cfg.files_per_job
    cuts_path = args.config / "cuts.json"

    if args.stage == "b":
        return _stage_b(args, cfg, cuts_path)

    # ---- discover -------------------------------------------------------
    listing = _read_listing(args.file_list) if args.file_list else None
    try:
        accepted, rejected = discover_inputs(cfg, runs, listing=listing)
    except (FarmConfigError, RunListError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    print(f"farm config    : {cfg.source}")
    print(f"target         : {cfg.target}  ({cfg.polarity})")
    print(f"cuts.json      : {cuts_path}  sha256 {config_sha256(cuts_path)[:16] if cuts_path.is_file() else '(missing)'}")
    print(f"inputs         : {len(cfg.inputs)} spec(s)")
    print(f"files accepted : {len(accepted)}")
    print(f"files rejected : {len(rejected)}")

    if rejected:
        by_reason: dict[str, int] = {}
        for r in rejected:
            by_reason[_generalise(r.reason)] = by_reason.get(_generalise(r.reason), 0) + 1
        print("\nREJECTED, by reason -- none of this is silent:")
        for reason, n in sorted(by_reason.items(), key=lambda kv: -kv[1]):
            print(f"  {n:6d}  {reason}")
        print("  (the run-list filter is the point: inbending and outbending runs share one")
        print("   directory tree, and the slim schema does not record polarity, so a scan")
        print("   without this filter mixes torus polarities and nothing downstream notices.)")

    if not accepted:
        print("\nerror: no input files survived. Nothing to submit.", file=sys.stderr)
        return 1

    chunks = chunk_inputs(accepted, fpj)
    if args.max_jobs:
        chunks = chunks[: args.max_jobs]
        print(f"\n--max-jobs {args.max_jobs}: synthesizing {len(chunks)} of the full set. THIS IS A TEST "
              f"SUBMISSION, not a production.")

    runs_seen = sorted({f.run for f in accepted})
    print(f"\nruns           : {len(runs_seen)}  ({runs_seen[0]}-{runs_seen[-1]})")
    print(f"files per job  : {fpj}")
    print(f"jobs           : {len(chunks)}")
    print(f"workflow       : {cfg.swif2.workflow}")
    print(f"resources      : {cfg.swif2.cores} core(s), {cfg.swif2.ram_gb}G ram, "
          f"{cfg.swif2.disk_gb}G disk, {format_swif2_time(cfg.swif2.time)}")
    if args.max_events:
        print(f"max events     : {args.max_events} PER FILE -- every output is a TRUNCATED PREFIX, "
              f"not a sample.\n                 Stamped into each output's provenance; Stage B refuses "
              f"such a file by default.")
    print(f"output         : {cfg.output_dir}/<run>/<job>/slim_<i>.root")
    print(f"logs           : {cfg.log_dir}/<job>.log")

    # ---- pre-flight -----------------------------------------------------
    exe_abs, exe_problems, exe_notes = _resolve_exe(args.exe)
    for n in exe_notes:
        print(f"\nnote: {n}")

    # Staged paths are a blocker only for --submit. A dry run on a laptop is a
    # legitimate way to read the script; submitting one built from laptop paths
    # is not.
    staged_bad = _check_staged_paths(cuts_path, args.scripts_dir)
    if staged_bad:
        where = "\n".join(f"        {b}" for b in staged_bad)
        msg = (
            "these paths are staged BY REFERENCE and are not on a filesystem the worker nodes\n"
            "      mount, so the jobs cannot read them:\n"
            f"{where}\n"
            "      SWIF2 copies -input FROM these paths, on the node, at job start. A script\n"
            "      generated here has this machine's paths baked in and is readable but NOT\n"
            "      submittable. Generate it where you will submit it -- run the same command on\n"
            "      ifarm from a checkout on /work. Do not edit the paths in the script."
        )
        if args.submit:
            exe_problems = exe_problems + [msg]
        else:
            print(f"\nnote: {msg}")

    blockers = _preflight(cuts_path, len(chunks))
    fallback_only = not exe_problems  # --allow-rga-fallback-production forgives the GBT gap, not a bad exe
    blockers = exe_problems + blockers
    if blockers and not (args.allow_rga_fallback_production and fallback_only):
        # Flush first: the summary above is on stdout and the refusal below is on
        # stderr, and an unflushed stdout puts the refusal before the numbers that
        # explain it.
        sys.stdout.flush()
        print("\n" + "=" * 76, file=sys.stderr)
        print("REFUSING TO SUBMIT", file=sys.stderr)
        print("=" * 76, file=sys.stderr)
        for b in blockers:
            print(f"  * {b}", file=sys.stderr)
        print("\nNothing was submitted. Fix the above. (--allow-rga-fallback-production forgives\n"
              "the photon-model gap; it does not forgive an executable the nodes cannot reach.)",
              file=sys.stderr)
        return 2

    with cuts_path.open() as f:
        if json.load(f).get("photon", {}).get("allow_rga_fallback", False):
            print("\n" + "!" * 76)
            print("photon.allow_rga_fallback IS TRUE. Every photon in this production will be scored")
            print("by an RG-A INBENDING PASS-1 model -- data it was not trained on. Every output will")
            print("carry gbt.fallback_used = TRUE, and the Python stage will refuse to publish from")
            print("it without --allow-unpublishable. This is a study, not a measurement.")
            print("!" * 76)

    for g in ("grid_A_q2_xb.json", "grid_B_z_pt2.json"):
        if _grid_is_placeholder(args.config / "binning" / g):
            print(f"\nnote: config/binning/{g} is still a PLACEHOLDER. That does not affect Stage A "
                  f"(which does not bin), but make_grid must run on these slims before Stage B.")

    # ---- synthesize -----------------------------------------------------
    scripts_dir = args.scripts_dir
    scripts_dir.mkdir(parents=True, exist_ok=True)
    for ch in chunks:
        job = f"{cfg.swif2.workflow}_{ch.name_suffix}"
        (scripts_dir / f"{job}.wrapper.sh").write_text(
            stagea_wrapper(cfg, job_inputs(ch), exe_abs, max_events=args.max_events))

    script = swif2_script(cfg, chunks, scripts_dir=scripts_dir, exe=exe_abs, cuts_path=cuts_path)
    script_path = scripts_dir / f"{cfg.swif2.workflow}.swif2.sh"
    script_path.write_text(script)
    script_path.chmod(0o755)
    print(f"\nwrote {script_path}")
    print(f"      {len(chunks)} wrapper(s) in {scripts_dir}/")

    if not args.submit:
        print("\nDRY RUN -- nothing submitted. Inspect the script above, then re-run with --submit.")
        print("\nsample add-job (job 0):")
        for line in script.splitlines():
            if line.startswith("swif2 add-job"):
                idx = script.splitlines().index(line)
                print("\n".join("  " + l for l in script.splitlines()[idx : idx + 12]))
                break
        return 0

    print(f"\nsubmitting {len(chunks)} jobs ...")
    r = subprocess.run(["bash", str(script_path)])
    if r.returncode != 0:
        print(f"error: swif2 script exited {r.returncode}", file=sys.stderr)
        return 1
    # These are swif2's real command names, checked against `swif2 -help` on
    # ifarm. There is no `list-jobs` -- clas-framework's docs name one, and it
    # errors with "is not a valid swif command".
    print(f"\nsubmitted. Monitor with:\n"
          f"  swif2 status      -workflow {cfg.swif2.workflow}\n"
          f"  swif2 diagnose    -workflow {cfg.swif2.workflow}   # why jobs are not progressing\n"
          f"  swif2 show-job    -workflow {cfg.swif2.workflow} -name <job>\n"
          f"  swif2 retry-jobs  -workflow {cfg.swif2.workflow}   # resubmit problem jobs\n"
          f"  swif2 cancel      -workflow {cfg.swif2.workflow} -delete")
    return 0


def _stage_b(args, cfg, cuts_path: Path) -> int:
    """Print (or run) the single Stage B command that chains every slim.

    Not submitted as a SWIF2 workflow: it is one job, and it needs every Stage A
    output, so there is nothing to schedule around.
    """
    slim_dir = args.slim_dir or Path(cfg.output_dir)
    if not slim_dir.is_dir():
        print(f"error: no slim directory at {slim_dir}. Pass --slim-dir, or run --stage a first.",
              file=sys.stderr)
        return 1
    slims = sorted(str(p) for p in Path(slim_dir).rglob("slim_*.root"))
    if not slims:
        print(f"error: no slim_*.root under {slim_dir}. Stage A has not produced anything yet.",
              file=sys.stderr)
        return 1

    print(f"target      : {cfg.target}")
    print(f"slims found : {len(slims)}  under {slim_dir}")
    print(f"\nStage B runs ONCE over all of them, so the donor pool is drawn from the whole target")
    print(f"rather than one file at a time. stageB_bin validates that every input agrees on target,")
    print(f"config hash and photon model, and REFUSES a mixed set -- chaining an LD2 slim into an Sn")
    print(f"run would put Sn photons in LD2's mixed background and quietly corrupt R_A.")

    args_file = Path(args.scripts_dir) / f"stageB_{cfg.target}.inputs.txt"
    args_file.parent.mkdir(parents=True, exist_ok=True)
    args_file.write_text("\n".join(slims) + "\n")

    cmd = [args.stageb_exe]
    for s in slims:
        cmd += ["--input", s]
    cmd += ["--output", f"{cfg.output_dir}/binned_{cfg.target}.root",
            "--config", str(cuts_path),
            "--grid-a", str(args.config / "binning" / "grid_A_q2_xb.json"),
            "--grid-b", str(args.config / "binning" / "grid_B_z_pt2.json")]

    print(f"\ninput list written to {args_file} ({len(slims)} paths)")
    print("\ncommand:")
    print(f"  {args.stageb_exe} \\")
    print(f"      $(sed 's/^/--input /' {args_file} | tr '\\n' ' ') \\")
    print(f"      --output {cfg.output_dir}/binned_{cfg.target}.root \\")
    print(f"      --config {cuts_path} \\")
    print(f"      --grid-a {args.config}/binning/grid_A_q2_xb.json \\")
    print(f"      --grid-b {args.config}/binning/grid_B_z_pt2.json")

    if not args.submit:
        print("\nDRY RUN -- not executed. Re-run with --submit to run it here.")
        return 0
    print(f"\nrunning Stage B over {len(slims)} slims. This is single-threaded by design (a")
    print("multi-threaded sum of doubles would make sum_q2 depend on the thread count).")
    return subprocess.run(cmd).returncode


if __name__ == "__main__":
    raise SystemExit(main())
