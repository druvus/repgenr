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
