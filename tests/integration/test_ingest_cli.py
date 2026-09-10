"""`repgenr ingest` wiring, resume, and the local pipeline chain in `status`."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli import base as cli
from repgenr.cli.main import app
from repgenr.core.config import Config

_runner = CliRunner()
_SEQ = ">s\n" + "ACGT" * 10 + "\n"


def _source(tmp_path: Path, names: list[str]) -> Path:
    src = tmp_path / "src"
    src.mkdir(exist_ok=True)
    for name in names:
        (src / name).write_text(_SEQ)
    return src


def test_ingest_cli_populates_workdir_and_records_stage(tmp_path: Path) -> None:
    src = _source(tmp_path, ["Fam_Gen_sp1_GCA_000001.1.fasta", "Fam_Gen_sp2_GCA_000002.1.fasta"])
    wd = tmp_path / "wd"

    result = _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(src)])

    assert result.exit_code == 0, result.output
    assert (wd / "selection.tsv").exists()
    assert (wd / "manifest.sqlite").exists()
    assert sorted(p.name for p in (wd / "genomes").iterdir()) == [
        "Fam_Gen_sp1_GCA_000001.1.fasta",
        "Fam_Gen_sp2_GCA_000002.1.fasta",
    ]
    rec = Config.load(wd).stages["ingest"]
    assert rec.completed
    assert rec.inputs, "ingest must declare its source directory as a resume input"


def test_ingest_cli_resume_skips_then_reruns_on_source_change(tmp_path: Path, monkeypatch) -> None:
    src = _source(tmp_path, ["a.fasta"])
    wd = tmp_path / "wd"
    monkeypatch.setitem(cli._RUN_STATE, "force", False)

    first = _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(src)])
    assert first.exit_code == 0, first.output
    stamp = Config.load(wd).stages["ingest"].completed

    second = _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(src)])
    assert second.exit_code == 0, second.output
    assert Config.load(wd).stages["ingest"].completed == stamp, "unchanged source must skip"

    (src / "b.fasta").write_text(_SEQ)
    third = _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(src)])
    assert third.exit_code == 0, third.output
    assert sorted(p.name for p in (wd / "genomes").iterdir()) == ["a.fasta", "b.fasta"]


def test_ingest_cli_rejects_missing_source(tmp_path: Path) -> None:
    result = _runner.invoke(
        app, ["ingest", "-wd", str(tmp_path / "wd"), "--genomes-dir", str(tmp_path / "nope")]
    )
    assert result.exit_code != 0


def test_status_reports_local_chain_after_ingest(tmp_path: Path) -> None:
    cfg = Config()
    cfg.record_stage("ingest", completed="2026-01-01T00:00:00")
    cfg.save(tmp_path)

    result = _runner.invoke(app, ["status", "-wd", str(tmp_path)])

    assert result.exit_code == 0
    assert "Pipeline: local" in result.stdout
    assert "[done]    ingest" in result.stdout
    assert "[next] dereplicate" in result.stdout
    assert "metadata" not in result.stdout


def test_doctor_accepts_ingested_workdir(tmp_path: Path) -> None:
    src = _source(tmp_path, ["a.fasta", "b.fasta"])
    wd = tmp_path / "wd"
    assert _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(src)]).exit_code == 0

    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])

    assert result.exit_code == 0, result.output
    assert "0 failure(s)" in result.output
