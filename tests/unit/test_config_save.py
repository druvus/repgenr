"""Config.save durability: the provenance file must never be truncated.

``repgenr.yaml`` carries every stage's provenance and resume fingerprint, so a
crash mid-write must leave the previous version intact (write-then-rename).
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from repgenr.core import config as config_mod
from repgenr.core.config import CONFIG_FILENAME, Config


def test_save_load_round_trip(tmp_path: Path) -> None:
    cfg = Config(repgenr_version="2.0.0")
    cfg.record_stage("metadata", tool="gtdb", params={"release": "232.0"})
    cfg.save(tmp_path)

    loaded = Config.load(tmp_path)
    assert loaded.stages["metadata"].tool == "gtdb"
    assert loaded.stages["metadata"].params == {"release": "232.0"}


def test_crash_during_save_preserves_previous_config(tmp_path: Path, monkeypatch) -> None:
    cfg = Config(repgenr_version="2.0.0")
    cfg.record_stage("metadata", tool="gtdb", params={})
    cfg.save(tmp_path)
    before = (tmp_path / CONFIG_FILENAME).read_text()

    def boom(*args, **kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr(config_mod.yaml, "safe_dump", boom)
    cfg.record_stage("genome", tool="datasets", params={})
    with pytest.raises(RuntimeError):
        cfg.save(tmp_path)

    # The previous file is untouched and still parseable.
    assert (tmp_path / CONFIG_FILENAME).read_text() == before
    assert yaml.safe_load(before)["stages"].keys() == {"metadata"}
    # No stray temp files are left behind.
    leftovers = [p.name for p in tmp_path.iterdir() if p.name != CONFIG_FILENAME]
    assert leftovers == []


def test_fingerprint_and_inputs_round_trip(tmp_path: Path) -> None:
    cfg = Config(repgenr_version="2.0.0")
    cfg.record_stage(
        "dereplicate", tool="skder", params={"secondary_ani": 0.99},
        fingerprint="abc123", inputs={"genomes": "d1", "selection.tsv": "d2"},
    )
    cfg.save(tmp_path)

    loaded = Config.load(tmp_path)
    rec = loaded.stages["dereplicate"]
    assert rec.fingerprint == "abc123"
    assert rec.inputs == {"genomes": "d1", "selection.tsv": "d2"}


def test_legacy_yaml_without_inputs_loads(tmp_path: Path) -> None:
    (tmp_path / CONFIG_FILENAME).write_text(
        "schema_version: 1\n"
        "repgenr_version: 2.0.0\n"
        "stages:\n"
        "  metadata:\n"
        "    tool: gtdb\n"
        "    completed: '2026-01-01T00:00:00'\n",
        encoding="utf-8",
    )
    loaded = Config.load(tmp_path)
    rec = loaded.stages["metadata"]
    assert rec.inputs == {}
    assert rec.fingerprint is None
