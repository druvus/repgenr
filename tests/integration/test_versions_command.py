"""`repgenr versions` surfaces recorded tool versions from repgenr.yaml."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.config import Config

_runner = CliRunner()


def _workdir_with_versions(tmp_path: Path) -> Path:
    cfg = Config()
    cfg.record_stage("vmetadata", tool_versions={"datasets": "16.0.0"}, completed="t")
    cfg.record_stage("vgenome", tool_versions={"mashtree": "1.4.6"}, completed="t")
    cfg.save(tmp_path)
    return tmp_path


def test_versions_stdout(tmp_path: Path) -> None:
    wd = _workdir_with_versions(tmp_path)
    result = _runner.invoke(app, ["versions", "-wd", str(wd)])
    assert result.exit_code == 0
    assert "datasets: 16.0.0" in result.stdout
    assert "mashtree: 1.4.6" in result.stdout


def test_versions_fragment_file(tmp_path: Path) -> None:
    wd = _workdir_with_versions(tmp_path)
    out = tmp_path / "frag.yml"
    result = _runner.invoke(app, ["versions", "-wd", str(wd), "--versions-out", str(out)])
    assert result.exit_code == 0
    # 4-space-indented, sorted -> slots under a process key in versions.yml
    assert out.read_text() == "    datasets: 16.0.0\n    mashtree: 1.4.6\n"


def test_versions_empty_workdir(tmp_path: Path) -> None:
    Config().save(tmp_path)  # no stages recorded
    out = tmp_path / "frag.yml"
    result = _runner.invoke(app, ["versions", "-wd", str(tmp_path), "--versions-out", str(out)])
    assert result.exit_code == 0
    assert out.read_text() == ""  # nothing recorded -> empty fragment
