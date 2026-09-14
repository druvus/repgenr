"""`repgenr reads` wiring and the reads chain in `status`."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli import cmd_reads
from repgenr.cli.main import app
from repgenr.core.config import Config

_runner = CliRunner()


def test_reads_command_builds_its_params(monkeypatch, tmp_path: Path) -> None:
    seen = {}

    def fake_run(stage, workdir, build, *, create=False):
        seen["stage"], seen["params"] = stage, build()

    monkeypatch.setattr(cmd_reads, "_run", fake_run)
    listing = tmp_path / "acc.txt"
    listing.write_text("SRR1\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            *("reads", "-wd", str(tmp_path / "wd"), "-ts", "Francisella tularensis"),
            *("--accession", "SRR2", "--accession-file", str(listing)),
            *("--platform", "ont", "--max-runs", "5", "--min-bases", "100", "--all-runs"),
        ],
    )
    assert result.exit_code == 0, result.output
    p = seen["params"]
    assert seen["stage"] == "reads"
    assert p.target_species == "Francisella tularensis" and p.accessions == ["SRR2"]
    assert p.accession_file == str(listing) and p.platform == "ont"
    assert (p.max_runs, p.min_bases, p.one_per_sample) == (5, 100, False)


def test_reads_command_rejects_an_unknown_platform(monkeypatch, tmp_path: Path) -> None:
    result = _runner.invoke(
        app, ["reads", "-wd", str(tmp_path / "wd"), "-tg", "x", "--platform", "solid"]
    )
    assert result.exit_code != 0
    assert "--platform" in result.output + str(result.exception or "")


def test_status_reports_the_reads_chain(tmp_path: Path) -> None:
    wd = tmp_path / "wd"
    wd.mkdir()
    cfg = Config()
    cfg.record_stage("reads", params={"taxid": "263"}, completed="t")
    cfg.save(wd)
    result = _runner.invoke(app, ["status", "-wd", str(wd)])
    assert result.exit_code == 0, result.output
    assert "Pipeline: reads" in result.output
    assert "[next] assemble" in result.output
