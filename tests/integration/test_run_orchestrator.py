"""The `repgenr run` orchestrator chains the canonical stages in order."""

from __future__ import annotations

from typer.testing import CliRunner

from repgenr.cli import cmd_run
from repgenr.cli.main import app

_runner = CliRunner()


def _record(monkeypatch) -> list[str]:
    calls: list[str] = []

    def fake_run(stage, workdir, build, *, create=False):
        build()  # exercise the param builder (catches bad kwargs)
        calls.append(stage)

    monkeypatch.setattr(cmd_run, "_run", fake_run)
    return calls


def test_run_bacterial_chain(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(app, [
        "run", "-wd", str(tmp_path), "-d", "rep", "-l", "genus", "-tg", "francisella",
    ])
    assert result.exit_code == 0, result.stdout
    assert calls == ["metadata", "genome", "dereplicate", "phylo", "tree2tax"]
    assert "Pipeline complete" in result.stdout


def test_run_viral_chain(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(app, [
        "run", "-wd", str(tmp_path), "--viral", "-t", "mastadenovirus", "-tg", "Mastadenovirus",
    ])
    assert result.exit_code == 0, result.stdout
    assert calls == ["vmetadata", "vgenome", "dereplicate", "phylo", "tree2tax"]


def test_run_validates_tool(monkeypatch, tmp_path) -> None:
    _record(monkeypatch)
    result = _runner.invoke(app, ["run", "-wd", str(tmp_path), "--tool", "bogus"])
    assert result.exit_code != 0
