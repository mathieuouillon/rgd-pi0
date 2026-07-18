"""Stage A provenance aggregation in :mod:`pi0.io`.

Stage B writes one ``provenance_stageA_NNN`` directory per input file, unmerged
by design; there is no single ``provenance_stageA``. ``_read_stage_a`` reads all
of them and the publish guard must fire on a fallback or truncation in ANY
input. These exercise that with a fake uproot directory tree (no ROOT needed) --
the seam a real multi-input Stage B file is the first thing to test.
"""

from __future__ import annotations

import pytest

from pi0.io import ProvenanceError, _read_stage_a


class _Named:
    def __init__(self, name: str, title: str) -> None:
        self.all_members = {"fName": name, "fTitle": title}


class _Dir:
    """A stand-in for an uproot directory of TNamed objects."""

    def __init__(self, entries: dict[str, str]) -> None:
        self._objs = {n: _Named(n, t) for n, t in entries.items()}

    def keys(self, recursive: bool = False) -> list[str]:
        return [f"{k};1" for k in self._objs]

    def __getitem__(self, key: str) -> _Named:
        return self._objs[key.split(";")[0]]


class _File:
    """A stand-in for an uproot ReadOnlyDirectory keyed by directory name."""

    file_path = "fake.root"

    def __init__(self, dirs: dict[str, dict[str, str]]) -> None:
        self._dirs = {n: _Dir(e) for n, e in dirs.items()}

    def keys(self, recursive: bool = False) -> list[str]:
        return [f"{k};1" for k in self._dirs]

    def __getitem__(self, key: str) -> _Dir:
        return self._dirs[key.split(";")[0]]


def test_read_stage_a_flags_fallback_in_any_input() -> None:
    """A fallback in one of several inputs taints the aggregate."""
    f = _File(
        {
            "provenance": {"config.sha256": "abc"},
            "provenance_stageA_000": {"gbt.fallback_used": "false", "events.max_events_requested": "all"},
            "provenance_stageA_001": {
                "gbt.fallback_used": "TRUE -- photons scored by a model trained on OTHER data",
                "gbt.model": "RGA_inbending_pass1",
                "events.max_events_requested": "all",
            },
        }
    )
    pa = _read_stage_a(f)
    assert pa.starts_true("gbt.fallback_used")
    assert pa.get("gbt.model") == "RGA_inbending_pass1"


def test_read_stage_a_flags_truncation_in_any_input() -> None:
    """A truncated input taints the aggregate even if others are complete."""
    f = _File(
        {
            "provenance_stageA_000": {"gbt.fallback_used": "false", "events.max_events_requested": "all"},
            "provenance_stageA_001": {
                "gbt.fallback_used": "false",
                "events.max_events_requested": "2000000 -- TRUNCATED RUN",
            },
        }
    )
    pa = _read_stage_a(f)
    assert "TRUNCATED" in (pa.get("events.max_events_requested") or "").upper()
    assert not pa.starts_true("gbt.fallback_used")


def test_read_stage_a_clean_when_all_inputs_clean() -> None:
    f = _File(
        {
            "provenance_stageA_000": {"gbt.fallback_used": "false", "events.max_events_requested": "all"},
            "provenance_stageA_001": {"gbt.fallback_used": "false", "events.max_events_requested": "all"},
        }
    )
    pa = _read_stage_a(f)
    assert not pa.starts_true("gbt.fallback_used")
    assert "TRUNCATED" not in (pa.get("events.max_events_requested") or "").upper()


def test_read_stage_a_missing_block_raises() -> None:
    """A file with no per-input Stage A provenance is refused, not silently accepted."""
    f = _File({"provenance": {"config.sha256": "abc"}, "provenance_text": {}})
    with pytest.raises(ProvenanceError, match="provenance_stageA_NNN"):
        _read_stage_a(f)
