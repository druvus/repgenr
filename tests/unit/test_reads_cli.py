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
            *("--max-bases", "5000000000"),
        ],
    )
    assert result.exit_code == 0, result.output
    p = seen["params"]
    assert seen["stage"] == "reads"
    assert p.target_species == "Francisella tularensis" and p.accessions == ["SRR2"]
    assert p.accession_file == str(listing) and p.platform == "ont"
    assert (p.max_runs, p.min_bases, p.one_per_sample) == (5, 100, False)
    assert p.max_bases == 5_000_000_000


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


def test_assemble_command_builds_its_params(monkeypatch, tmp_path: Path) -> None:
    seen = {}

    def fake_run(stage, workdir, build, *, create=False):
        seen["stage"], seen["params"] = stage, build()

    monkeypatch.setattr(cmd_reads, "_run", fake_run)
    og = tmp_path / "og.fasta"
    og.write_text(">o\nACGT\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            *("assemble", "-wd", str(tmp_path / "wd"), "--assembler", "flye", "-t", "8"),
            *("--jobs", "1", "--memory-gb", "32", "--min-contig-length", "1000"),
            *("--outgroup", str(og), "--keep-reads", "--tool-arg", "mode=nano-raw"),
        ],
    )
    assert result.exit_code == 0, result.output
    p = seen["params"]
    assert seen["stage"] == "assemble" and p.assembler == "flye"
    assert (p.threads, p.jobs, p.memory_gb, p.min_contig_length) == (8, 1, 32, 1000)
    assert p.outgroup == str(og) and p.keep_reads and p.extra == {"mode": "nano-raw"}


def test_assemble_command_rejects_an_unknown_assembler(tmp_path: Path) -> None:
    result = _runner.invoke(app, ["assemble", "-wd", str(tmp_path / "wd"), "--assembler", "velvet"])
    assert result.exit_code != 0
    assert "--assembler" in result.output + str(result.exception or "")
