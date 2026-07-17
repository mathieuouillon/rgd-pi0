"""Tests for pi0.farm -- run lists, input discovery, chunking, SWIF2 synthesis.

The run-list filter is the safety-critical part and gets the most attention: it
is the only thing standing between a recursive directory scan and a production
that silently mixes torus polarities.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from pi0.farm import (
    TARGET_TO_RUNLIST,
    Chunk,
    FarmConfigError,
    InputFile,
    chunk_inputs,
    discover_inputs,
    format_swif2_time,
    load_farm_config,
    load_runs,
    run_of,
    runs_for,
    stagea_wrapper,
    swif2_script,
)

CONFIG_DIR = Path(__file__).resolve().parents[2] / "config"


def _dst(run: int, chunk: str = "00000-00004") -> str:
    return (
        f"/mss/clas12/rg-d/production/pass1/recon/LD2/dst/recon/"
        f"{run:06d}/rec_clas_{run:06d}.evio.{chunk}.hipo"
    )


# ---------------------------------------------------------------------------
# Run lists
# ---------------------------------------------------------------------------


def test_runs_json_matches_the_note():
    """tab:runs in the analysis note. If these drift, one of the two is wrong."""
    runs = load_runs(CONFIG_DIR)
    assert len(runs["outbending"]["LD2"]) == 135
    assert len(runs["outbending"]["CxC"]) == 112
    assert len(runs["outbending"]["CuSn"]) == 268
    assert len(runs["inbending"]["LD2"]) == 24
    assert len(runs["inbending"]["CxC"]) == 8
    assert len(runs["inbending"]["CuSn"]) == 27


def test_outbending_lists_are_pairwise_disjoint():
    """Each run carried a single target -- the note's claim, checked."""
    runs = load_runs(CONFIG_DIR)
    ld2 = set(runs["outbending"]["LD2"])
    cxc = set(runs["outbending"]["CxC"])
    cusn = set(runs["outbending"]["CuSn"])
    assert not (ld2 & cxc)
    assert not (ld2 & cusn)
    assert not (cxc & cusn)
    assert len(ld2 | cxc | cusn) == 515


def test_cu_and_sn_share_the_cusn_run_list():
    """Two foils, one assembly, one set of runs, separated offline by vz."""
    runs = load_runs(CONFIG_DIR)
    assert TARGET_TO_RUNLIST["Cu"] == TARGET_TO_RUNLIST["Sn"] == "CuSn"
    assert runs_for(runs, "outbending", "Cu") == runs_for(runs, "outbending", "Sn")
    assert runs_for(runs, "outbending", "Cu") != runs_for(runs, "outbending", "LD2")


def test_inbending_and_outbending_ld2_do_not_overlap():
    runs = load_runs(CONFIG_DIR)
    assert not (runs_for(runs, "inbending", "LD2") & runs_for(runs, "outbending", "LD2"))


def test_unknown_target_and_polarity_refuse():
    runs = load_runs(CONFIG_DIR)
    with pytest.raises(FarmConfigError, match="unknown target"):
        runs_for(runs, "outbending", "CuSn")  # a run list, not a target
    with pytest.raises(FarmConfigError, match="unknown polarity"):
        runs_for(runs, "sideways", "LD2")


def test_run_of_reads_filename_then_parent_dir():
    assert run_of(_dst(18419)) == 18419
    assert run_of("/mss/x/dst/recon/018419/something_else.hipo") == 18419
    assert run_of("/tmp/no_run_here.hipo") is None


# ---------------------------------------------------------------------------
# The filter -- the whole point
# ---------------------------------------------------------------------------


def _cfg(tmp_path: Path, target: str = "LD2", polarity: str = "outbending", **over) -> Path:
    doc = {
        "farm": {
            "target": target,
            "polarity": polarity,
            "inputs": ["/mss/nonexistent"],
            "output_dir": str(tmp_path / "out"),
            "log_dir": str(tmp_path / "out" / "logs"),
            "files_per_job": 5,
            "swif2": {"workflow": f"wf_{target}", "time": "04:00:00",
                      "modules": {"load": ["root/6.36.04"]}, "tags": {"stream": "pi0"}},
            **over,
        }
    }
    p = tmp_path / f"{target}.farm.json"
    p.write_text(json.dumps(doc))
    return p


def test_inbending_files_are_rejected_from_an_outbending_production(tmp_path):
    """THE test. Inbending and outbending LD2 share one /mss tree; the slim schema
    does not record polarity, so nothing downstream could catch this."""
    runs = load_runs(CONFIG_DIR)
    cfg = load_farm_config(_cfg(tmp_path))
    listing = [_dst(18305), _dst(18306), _dst(18419), _dst(18420)]  # 2 inbending, 2 outbending
    acc, rej = discover_inputs(cfg, runs, listing=listing)
    assert [f.run for f in acc] == [18419, 18420]
    assert len(rej) == 2
    assert all("not in the outbending LD2 run list" in r.reason for r in rej)


def test_a_run_in_neither_list_is_rejected(tmp_path):
    """18432 is LD2_trigger_rgd_v2_2_Q2_2_5 -- a Q2>2.5 trigger config, and in
    NEITHER the inbending nor the outbending list. A Q2-biased trigger sculpts the
    Q2 spectrum, and Q2 is a binning axis."""
    runs = load_runs(CONFIG_DIR)
    assert 18432 not in runs_for(runs, "outbending", "LD2")
    assert 18432 not in runs_for(runs, "inbending", "LD2")
    assert 18432 in runs["LD2_trigger_sets"]["LD2_trigger_rgd_v2_2_Q2_2_5"]
    cfg = load_farm_config(_cfg(tmp_path))
    acc, rej = discover_inputs(cfg, runs, listing=[_dst(18432)])
    assert acc == []
    assert len(rej) == 1


def test_exclude_runs_is_honoured(tmp_path):
    runs = load_runs(CONFIG_DIR)
    cfg = load_farm_config(_cfg(tmp_path, exclude_runs=[18419]))
    acc, rej = discover_inputs(cfg, runs, listing=[_dst(18419), _dst(18420)])
    assert [f.run for f in acc] == [18420]
    assert "exclude_runs" in rej[0].reason


def test_a_file_with_no_run_number_is_rejected_not_guessed(tmp_path):
    runs = load_runs(CONFIG_DIR)
    cfg = load_farm_config(_cfg(tmp_path))
    acc, rej = discover_inputs(cfg, runs, listing=["/tmp/mystery.hipo"])
    assert acc == []
    assert "no run number" in rej[0].reason


def test_duplicates_are_rejected_not_double_counted(tmp_path):
    runs = load_runs(CONFIG_DIR)
    cfg = load_farm_config(_cfg(tmp_path))
    acc, rej = discover_inputs(cfg, runs, listing=[_dst(18419), _dst(18419)])
    assert len(acc) == 1
    assert any("duplicate" in r.reason for r in rej)


def test_discovery_order_is_canonical_not_listing_order(tmp_path):
    """The donor pool downstream is order-dependent (reservoir sampling is a
    function of the sequence of offers), so discovery must not inherit whatever
    order the filesystem or the .dat file happened to have."""
    runs = load_runs(CONFIG_DIR)
    cfg = load_farm_config(_cfg(tmp_path))
    a, _ = discover_inputs(cfg, runs, listing=[_dst(18420), _dst(18419)])
    b, _ = discover_inputs(cfg, runs, listing=[_dst(18419), _dst(18420)])
    assert [f.path for f in a] == [f.path for f in b]
    assert [f.run for f in a] == [18419, 18420]


def test_cu_and_sn_accept_the_same_files(tmp_path):
    """They read the same runs and differ only in the vz window."""
    runs = load_runs(CONFIG_DIR)
    listing = [_dst(18561), _dst(18563)]
    cu, _ = discover_inputs(load_farm_config(_cfg(tmp_path, "Cu")), runs, listing=listing)
    sn, _ = discover_inputs(load_farm_config(_cfg(tmp_path, "Sn")), runs, listing=listing)
    assert [f.path for f in cu] == [f.path for f in sn] == listing


# ---------------------------------------------------------------------------
# Input expansion
# ---------------------------------------------------------------------------


@pytest.fixture
def hipo_tree(tmp_path):
    (tmp_path / "sub").mkdir()
    for rel in ("a.hipo", "b.hipo", "sub/c.hipo"):
        (tmp_path / rel).touch()
    (tmp_path / "not_hipo.txt").touch()
    return tmp_path


def test_directory_is_scanned_recursively_for_hipo_only(hipo_tree):
    from pi0.farm import _expand_one
    found = _expand_one(str(hipo_tree))
    assert [Path(f).name for f in found] == ["a.hipo", "b.hipo", "c.hipo"]


def test_literal_file_is_used_as_is(hipo_tree):
    from pi0.farm import _expand_one
    assert _expand_one(str(hipo_tree / "a.hipo")) == [str(hipo_tree / "a.hipo")]


def test_absolute_glob(hipo_tree):
    from pi0.farm import _expand_one
    assert [Path(f).name for f in _expand_one(f"{hipo_tree}/*.hipo")] == ["a.hipo", "b.hipo"]


def test_recursive_glob(hipo_tree):
    from pi0.farm import _expand_one
    assert [Path(f).name for f in _expand_one(f"{hipo_tree}/**/*.hipo")] == \
        ["a.hipo", "b.hipo", "c.hipo"]


def test_relative_glob_resolves_against_cwd(hipo_tree, monkeypatch):
    """A relative pattern has no directory part to split on. Hand-rolled path
    arithmetic returns nothing here instead of erroring; glob.glob does not."""
    from pi0.farm import _expand_one
    monkeypatch.chdir(hipo_tree)
    assert [Path(f).name for f in _expand_one("*.hipo")] == ["a.hipo", "b.hipo"]


def test_a_glob_matching_nothing_is_reported_not_faked(tmp_path):
    """clas-framework pushes an unmatched glob through verbatim as a fake path, so
    the job gets a nonexistent -input and fails on the node."""
    from pi0.farm import _expand_one
    assert _expand_one(str(tmp_path / "nope" / "*.hipo")) == []


def test_a_missing_input_spec_is_reported_with_a_reason(tmp_path):
    """clas-framework's recursive_directory_iterator THROWS on a missing directory
    and the throw is uncaught, killing the tool. This reports instead."""
    runs = load_runs(CONFIG_DIR)
    p = _cfg(tmp_path)
    doc = json.loads(p.read_text())
    doc["farm"]["inputs"] = [str(tmp_path / "does_not_exist")]
    p.write_text(json.dumps(doc))
    cfg = load_farm_config(p)
    acc, rej = discover_inputs(cfg, runs)
    assert acc == []
    assert any("resolved to no files" in r.reason for r in rej)


# ---------------------------------------------------------------------------
# Chunking
# ---------------------------------------------------------------------------


def test_chunks_never_span_runs():
    files = [InputFile(_dst(18419, f"{i:05d}"), 18419) for i in range(7)] + \
            [InputFile(_dst(18420, f"{i:05d}"), 18420) for i in range(3)]
    chunks = chunk_inputs(files, 5)
    for c in chunks:
        assert len({f.run for f in c.files}) == 1
    assert [len(c.files) for c in chunks] == [5, 2, 3]


def test_chunk_indices_are_contiguous_from_zero():
    files = [InputFile(_dst(18419, f"{i:05d}"), 18419) for i in range(11)]
    chunks = chunk_inputs(files, 4)
    assert [c.index for c in chunks] == [0, 1, 2]
    assert sum(len(c.files) for c in chunks) == 11


def test_files_per_job_must_be_positive():
    with pytest.raises(FarmConfigError):
        chunk_inputs([InputFile(_dst(18419), 18419)], 0)


# ---------------------------------------------------------------------------
# Time conversion
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("given,want", [
    ("04:00:00", "14400s"), ("16:00:00", "57600s"), ("00:30:00", "1800s"),
    ("6:00", "21600s"), ("900s", "900s"),
])
def test_format_swif2_time(given, want):
    assert format_swif2_time(given) == want


@pytest.mark.parametrize("bad", ["4h", "", "later", "-1:00:00", "00:00:00", "01:99:00"])
def test_malformed_time_refuses_here_not_on_the_farm(bad):
    """clas-framework returns malformed input unchanged and lets SWIF2 complain
    after submission. Refusing at parse time is the whole difference."""
    with pytest.raises(FarmConfigError):
        format_swif2_time(bad)


# ---------------------------------------------------------------------------
# Config schema
# ---------------------------------------------------------------------------


def test_shipped_farm_configs_all_load():
    for t in ("LD2", "CxC", "Cu", "Sn"):
        cfg = load_farm_config(CONFIG_DIR / "farm" / f"{t}.farm.json")
        assert cfg.target == t
        assert cfg.polarity == "outbending"
        assert format_swif2_time(cfg.swif2.time).endswith("s")


def test_shipped_farm_configs_carry_no_cut_values():
    """Cut values live in cuts.json. A threshold here is a bug."""
    forbidden = ("gbt_threshold", "vz_min", "vz_max", "q2_min", "w_min", "y_max",
                 "sampling_fraction", "mass_window")
    for t in ("LD2", "CxC", "Cu", "Sn"):
        blob = (CONFIG_DIR / "farm" / f"{t}.farm.json").read_text()
        for key in forbidden:
            assert key not in blob, f"{t}.farm.json carries a cut value: {key}"


def test_bare_modules_array_is_refused(tmp_path):
    """A bare array parses fine and silently loads nothing -- the job then fails on
    the node for want of ROOT, not at submit."""
    doc = json.loads(_cfg(tmp_path).read_text())
    doc["farm"]["swif2"]["modules"] = ["root/6.36.04"]
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(doc))
    with pytest.raises(FarmConfigError, match="load"):
        load_farm_config(p)


def test_missing_required_key_refuses_rather_than_defaults(tmp_path):
    doc = json.loads(_cfg(tmp_path).read_text())
    del doc["farm"]["output_dir"]
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(doc))
    with pytest.raises(FarmConfigError, match="output_dir"):
        load_farm_config(p)


def test_target_must_be_a_skim_target_not_a_runlist(tmp_path):
    doc = json.loads(_cfg(tmp_path).read_text())
    doc["farm"]["target"] = "CuSn"
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(doc))
    with pytest.raises(FarmConfigError, match="not a stageA_skim target"):
        load_farm_config(p)


# ---------------------------------------------------------------------------
# SWIF2 synthesis
# ---------------------------------------------------------------------------


def _script(tmp_path, n_files=2, fpj=2):
    runs = load_runs(CONFIG_DIR)
    cfg = load_farm_config(_cfg(tmp_path))
    files = [InputFile(_dst(18419, f"{i:05d}"), 18419) for i in range(n_files)]
    chunks = chunk_inputs(files, fpj)
    return cfg, chunks, swif2_script(cfg, chunks, scripts_dir=tmp_path,
                                     exe="./stageA", cuts_path=CONFIG_DIR / "cuts.json")


def test_swif2_script_is_valid_bash(tmp_path):
    _, _, s = _script(tmp_path)
    p = tmp_path / "s.sh"
    p.write_text(s)
    assert subprocess.run(["bash", "-n", str(p)]).returncode == 0


def test_swif2_uses_trailing_positional_not_dash_shell(tmp_path):
    """clas-framework's docs show -shell; its code emits a trailing positional and
    its own test asserts -shell is absent. Ported from the code."""
    _, _, s = _script(tmp_path)
    assert "-shell" not in s
    assert "bash wrapper.sh" in s


def test_swif2_emits_create_addjob_run_in_order(tmp_path):
    _, chunks, s = _script(tmp_path, n_files=4, fpj=2)
    assert s.index("swif2 create") < s.index("swif2 add-job") < s.index("swif2 run")
    assert s.count("swif2 add-job") == len(chunks) == 2
    assert s.count("swif2 create") == 1 and s.count("swif2 run") == 1


def test_every_job_mkdirs_its_output_dir(tmp_path):
    """SWIF2's -output does not create the parent; without this the copy is lost."""
    _, chunks, s = _script(tmp_path, n_files=4, fpj=2)
    assert s.count("mkdir -p") == len(chunks)


def test_tags_are_emitted_in_sorted_order(tmp_path):
    """clas-framework stores tags in a std::map, so emission is alphabetical
    regardless of declaration order. Matched so a golden comparison holds."""
    runs = load_runs(CONFIG_DIR)
    p = _cfg(tmp_path)
    doc = json.loads(p.read_text())
    doc["farm"]["swif2"]["tags"] = {"zebra": "z", "alpha": "a"}
    p.write_text(json.dumps(doc))
    cfg = load_farm_config(p)
    chunks = chunk_inputs([InputFile(_dst(18419), 18419)], 1)
    s = swif2_script(cfg, chunks, scripts_dir=tmp_path, exe="./x", cuts_path=CONFIG_DIR / "cuts.json")
    assert s.index("-tag alpha") < s.index("-tag zebra")


def test_time_is_converted_in_the_emitted_script(tmp_path):
    _, _, s = _script(tmp_path)
    assert "-time 14400s" in s
    assert "04:00:00" not in s


def test_every_input_and_output_is_staged(tmp_path):
    _, _, s = _script(tmp_path, n_files=2, fpj=2)
    assert s.count("-input input_0.hipo") == 1
    assert s.count("-input input_1.hipo") == 1
    assert "-input 'cuts.json'" in s
    assert "-input 'wrapper.sh'" in s
    assert s.count("-output slim_0.root") == 1
    assert "-output 'log.txt'" in s


# ---------------------------------------------------------------------------
# The job wrapper
# ---------------------------------------------------------------------------


def test_wrapper_is_valid_bash(tmp_path):
    cfg = load_farm_config(_cfg(tmp_path))
    p = tmp_path / "w.sh"
    p.write_text(stagea_wrapper(cfg, 3, "./stageA"))
    assert subprocess.run(["bash", "-n", str(p)]).returncode == 0


def test_wrapper_runs_the_skim_once_per_staged_file(tmp_path):
    """stageA_skim takes ONE file per invocation, so the wrapper loops."""
    cfg = load_farm_config(_cfg(tmp_path))
    w = stagea_wrapper(cfg, 5, "./stageA")
    assert "seq 0 4" in w


def test_wrapper_does_not_set_e_so_one_bad_file_does_not_abandon_the_chunk(tmp_path):
    cfg = load_farm_config(_cfg(tmp_path))
    w = stagea_wrapper(cfg, 3, "./stageA")
    assert "set -uo pipefail" in w
    assert "set -euo pipefail" not in w
    assert "exit $rc" in w


def test_wrapper_explains_exit_3(tmp_path):
    """A wall of bare 'exit 3' is a mystery; on RG-D it is the expected refusal."""
    cfg = load_farm_config(_cfg(tmp_path))
    w = stagea_wrapper(cfg, 1, "./stageA")
    assert "$s -eq 3" in w
    assert "allow_rga_fallback" in w


def test_wrapper_runs_the_exe_in_a_job_dir_that_is_not_the_checkout(tmp_path):
    """THE regression test. A SWIF2 job's CWD is a scratch dir holding only the
    staged inputs. A relative --exe resolves to nothing there and every job exits
    127 -- which `bash -n`, a dry run, and reading the script all miss."""
    cfg = load_farm_config(_cfg(tmp_path))
    job = tmp_path / "jobdir"
    job.mkdir()
    (job / "input_0.hipo").touch()
    (job / "cuts.json").touch()

    # A relative exe, as the old default was.
    (job / "wrapper.sh").write_text(stagea_wrapper(cfg, 1, "./build/src/stageA_skim/stageA_skim"))
    r = subprocess.run(["bash", "wrapper.sh"], cwd=job, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"})
    assert r.returncode == 127, "a relative exe must fail in a job dir -- that is the bug"

    # An absolute exe that exists: the wrapper finds and runs it.
    fake = tmp_path / "stageA_fake"
    fake.write_text("#!/usr/bin/env bash\nexit 0\n")
    fake.chmod(0o755)
    (job / "wrapper.sh").write_text(stagea_wrapper(cfg, 1, str(fake)))
    r = subprocess.run(["bash", "wrapper.sh"], cwd=job, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"})
    assert r.returncode == 0, f"absolute exe should run: {r.stderr}"


def test_wrapper_propagates_the_skim_exit_code(tmp_path):
    cfg = load_farm_config(_cfg(tmp_path))
    job = tmp_path / "j2"
    job.mkdir()
    (job / "input_0.hipo").touch()
    fake = tmp_path / "stageA_exit3"
    fake.write_text("#!/usr/bin/env bash\nexit 3\n")
    fake.chmod(0o755)
    (job / "wrapper.sh").write_text(stagea_wrapper(cfg, 1, str(fake)))
    r = subprocess.run(["bash", "wrapper.sh"], cwd=job, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"})
    assert r.returncode == 3, "SWIF2 must see the job as failed"
    assert "allow_rga_fallback" in r.stderr, "exit 3 must explain itself"


def test_wrapper_missing_staged_input_does_not_pass_silently(tmp_path):
    cfg = load_farm_config(_cfg(tmp_path))
    job = tmp_path / "j3"
    job.mkdir()  # no input_0.hipo
    fake = tmp_path / "stageA_ok"
    fake.write_text("#!/usr/bin/env bash\nexit 0\n")
    fake.chmod(0o755)
    (job / "wrapper.sh").write_text(stagea_wrapper(cfg, 1, str(fake)))
    r = subprocess.run(["bash", "wrapper.sh"], cwd=job, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"})
    assert r.returncode == 4
    assert "MISSING STAGED INPUT" in r.stderr


def test_wrapper_one_bad_file_does_not_abandon_the_rest_of_the_chunk(tmp_path):
    """No `set -e`: file 0 failing must not skip files 1 and 2."""
    cfg = load_farm_config(_cfg(tmp_path))
    job = tmp_path / "j4"
    job.mkdir()
    for i in range(3):
        (job / f"input_{i}.hipo").touch()
    fake = tmp_path / "stageA_flaky"
    # Fails on input_0, succeeds on the rest; records every file it saw.
    fake.write_text('#!/usr/bin/env bash\necho "$2" >> seen.txt\n'
                    '[[ "$2" == "input_0.hipo" ]] && exit 3\nexit 0\n')
    fake.chmod(0o755)
    (job / "wrapper.sh").write_text(stagea_wrapper(cfg, 3, str(fake)))
    r = subprocess.run(["bash", "wrapper.sh"], cwd=job, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin"})
    seen = (job / "seen.txt").read_text().split()
    assert seen == ["input_0.hipo", "input_1.hipo", "input_2.hipo"], "all three must be attempted"
    assert r.returncode == 3, "and the job must still report failure"


# ---------------------------------------------------------------------------
# Executable resolution
# ---------------------------------------------------------------------------


def test_relative_exe_is_refused_with_the_resolved_path_offered():
    from pi0.batch import _resolve_exe
    absolute, problems, _ = _resolve_exe("./build/src/stageA_skim/stageA_skim")
    assert problems and "RELATIVE" in problems[0]
    assert Path(absolute).is_absolute()


def test_absolute_but_node_invisible_exe_is_refused():
    from pi0.batch import _resolve_exe
    _, problems, _ = _resolve_exe("/tmp/build/stageA_skim")
    assert problems and "not on a filesystem the worker nodes" in problems[0]


def test_node_visible_exe_that_exists_is_accepted(tmp_path, monkeypatch):
    from pi0 import batch
    fake = tmp_path / "stageA_skim"
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)
    monkeypatch.setattr(batch, "_NODE_VISIBLE", (str(tmp_path) + "/",))
    absolute, problems, notes = batch._resolve_exe(str(fake))
    assert problems == [] and notes == []
    assert absolute == str(fake)


def test_node_visible_exe_that_is_missing_is_refused_when_the_fs_is_mounted(tmp_path, monkeypatch):
    from pi0 import batch
    monkeypatch.setattr(batch, "_NODE_VISIBLE", (str(tmp_path) + "/",))
    _, problems, _ = batch._resolve_exe(str(tmp_path / "not_built_yet"))
    assert problems and "no executable at" in problems[0]


def test_staged_paths_off_a_node_visible_fs_are_reported(tmp_path):
    """cuts.json and the wrappers are staged BY REFERENCE: SWIF2 copies from those
    paths on the node. A script generated on a laptop has laptop paths baked in."""
    from pi0.batch import _check_staged_paths
    bad = _check_staged_paths(tmp_path / "cuts.json", tmp_path / "batch_scripts")
    assert len(bad) == 2
    assert any("cuts.json" in b for b in bad)
    assert any("wrapper" in b for b in bad)


def test_staged_paths_on_a_node_visible_fs_are_accepted(monkeypatch, tmp_path):
    from pi0 import batch
    monkeypatch.setattr(batch, "_NODE_VISIBLE", (str(tmp_path) + "/",))
    assert batch._check_staged_paths(tmp_path / "cuts.json", tmp_path / "scripts") == []


def test_unmounted_node_fs_is_a_note_not_a_refusal():
    """Generating a script on a laptop to submit from ifarm: /work is not mounted,
    so existence is unknowable. Say so; do not refuse, and do not pretend to know."""
    from pi0.batch import _resolve_exe
    assert not Path("/work").is_dir(), "this test assumes /work is not mounted here"
    absolute, problems, notes = _resolve_exe("/work/clas12/users/x/rgd-pi0/build/stageA_skim")
    assert problems == []
    assert notes and "not mounted on this machine" in notes[0]
    assert absolute == "/work/clas12/users/x/rgd-pi0/build/stageA_skim"


def test_wrapper_loads_modules_and_exports_env(tmp_path):
    p = _cfg(tmp_path)
    doc = json.loads(p.read_text())
    doc["farm"]["env"] = {"ROOT_HIST": "0"}
    p.write_text(json.dumps(doc))
    cfg = load_farm_config(p)
    w = stagea_wrapper(cfg, 1, "./stageA")
    assert "module load root/6.36.04" in w
    assert "export ROOT_HIST=0" in w
