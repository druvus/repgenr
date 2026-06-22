"""The `repgenr status` command reports pipeline progress from repgenr.yaml."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.config import Config

_runner = CliRunner()


def test_status_no_workdir(tmp_path: Path) -> None:
    result = _runner.invoke(app, ["status", "-wd", str(tmp_path / "missing")])
    assert result.exit_code == 0
    assert "No RepGenR run found" in result.stdout


def test_status_bacterial_progress(tmp_path: Path) -> None:
    cfg = Config()
    cfg.record_stage("metadata", completed="2026-01-01T00:00:00")
    cfg.record_stage("genome", completed="2026-01-01T00:01:00")
    cfg.save(tmp_path)

    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "Pipeline: bacterial" in result.stdout
    assert "[done]    metadata" in result.stdout
    assert "[done]    genome" in result.stdout
    # dereplicate is the first incomplete stage
    assert "[next] dereplicate" in result.stdout
    assert "Next: repgenr dereplicate" in result.stdout


def test_status_viral_detection_and_extras(tmp_path: Path) -> None:
    cfg = Config()
    cfg.record_stage("vmetadata", completed="2026-01-01T00:00:00")
    cfg.record_stage("vgenome", completed="2026-01-01T00:01:00")
    cfg.record_stage("snptype", tool="simple", completed="2026-01-01T00:02:00")
    cfg.save(tmp_path)

    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "Pipeline: viral" in result.stdout
    assert "optional stages run:" in result.stdout
    assert "snptype" in result.stdout


def test_status_all_complete(tmp_path: Path) -> None:
    cfg = Config()
    for stage in ("metadata", "genome", "dereplicate", "phylo", "tree2tax"):
        cfg.record_stage(stage, completed="2026-01-01T00:00:00")
    cfg.save(tmp_path)

    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert "All stages complete" in result.stdout
