"""The sketch step of the genome-writing stages and the ``repgenr sketch`` command,
with a fake sourmash (tests/conftest.py)."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core import sketches
from repgenr.core.config import Config
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import READS_TSV, list_fasta, record_name, write_reads
from repgenr.core.errors import MissingBinaryError
from repgenr.stages.ingest import IngestParams
from repgenr.stages.ingest import run as ingest_run

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "unit"))
from assemble_fakes import (  # noqa: E402
    read_row,
    register_fake_assembler,
    unregister_fake_assembler,
)

_LOG = logging.getLogger("test")
_runner = CliRunner()


def _genomes(root: Path, n: int = 3) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    for i in range(n):
        (root / f"Fam_Gen_sp_GCA_{i:06d}.1.fasta").write_text(
            f">g{i}\n{'ACGT' * (i + 3)}\n", encoding="utf-8"
        )
    return root


def _sketch_names(workdir: Path) -> set[str]:
    directory = workdir / "sketches"
    if not directory.is_dir():
        return set()
    return {p.name.removesuffix(".sig.zip") for p in directory.glob("[!.]*.sig.zip")}


def _log(workdir: Path) -> str:
    return (workdir / "repgenr.log").read_text(encoding="utf-8")


def _genome_names(workdir: Path) -> set[str]:
    return {
        record_name(p) for p in list_fasta(workdir / "genomes") + list_fasta(workdir / "outgroup")
    }


# --- ingest -----------------------------------------------------------------------


def test_ingest_sketches_every_genome_and_the_outgroup(tmp_path, fake_sourmash) -> None:
    src = _genomes(tmp_path / "src")
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    ingest_run(ctx, IngestParams(genomes_dir=str(src), outgroup="GCA_000002.1"))
    assert _sketch_names(ctx.workdir) == _genome_names(ctx.workdir)
    assert len(_sketch_names(ctx.workdir)) == 3
    record = ctx.config.stages["ingest"]
    assert record.params["sketches"]["state"] == "done"
    assert record.params["sketches"]["written"] == 3
    assert record.tool_versions == {"sourmash": "4.9.4"}

    # A narrower re-ingest prunes the sketch of the genome that left.
    (src / "Fam_Gen_sp_GCA_000001.1.fasta").unlink()
    fake_sourmash.calls.clear()
    ingest_run(ctx, IngestParams(genomes_dir=str(src)))
    assert _sketch_names(ctx.workdir) == {"Fam_Gen_sp_GCA_000000.1", "Fam_Gen_sp_GCA_000002.1"}
    assert fake_sourmash.calls == []  # the remaining sketches are current
    assert ctx.config.stages["ingest"].params["sketches"]["removed"] == 1


def test_ingest_no_sketch_and_no_sourmash(tmp_path) -> None:
    src = _genomes(tmp_path / "src")
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    ingest_run(ctx, IngestParams(genomes_dir=str(src)))
    assert not (ctx.workdir / "sketches").exists()
    assert ctx.config.stages["ingest"].params["sketches"]["state"] == "unavailable"
    assert "sourmash not found" in _log(ctx.workdir) and "--sketch" in _log(ctx.workdir)
    ingest_run(ctx, IngestParams(genomes_dir=str(src), sketch=False))
    assert ctx.config.stages["ingest"].params["sketches"]["state"] == "off"


def test_ingest_explicit_sketch_without_sourmash_changes_nothing(tmp_path, monkeypatch) -> None:
    def missing(caps):
        raise MissingBinaryError("Required binary 'sourmash' not found on PATH.")

    monkeypatch.setattr(sketches, "preflight", missing)
    src = _genomes(tmp_path / "src")
    wd = tmp_path / "wd"
    result = _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(src), "--sketch"])
    assert result.exit_code == MissingBinaryError.exit_code, result.output
    assert not (wd / "genomes").exists()
    assert "ingest" not in Config.load(wd).stages


def test_ingest_from_workdir_copies_current_sketches(tmp_path, fake_sourmash) -> None:
    src = _genomes(tmp_path / "src")
    first = WorkdirContext(tmp_path / "first", create=True)
    ingest_run(first, IngestParams(genomes_dir=str(src)))
    first.close()
    fake_sourmash.calls.clear()
    second = WorkdirContext(tmp_path / "second", create=True)
    ingest_run(second, IngestParams(from_workdirs=[str(first.workdir)]))
    assert fake_sourmash.calls == []
    assert second.config.stages["ingest"].params["sketches"]["copied"] == 3
    assert _sketch_names(second.workdir) == _genome_names(second.workdir)
    rec = second.manifest.sketch_records()["GCA_000000.1"]
    assert rec.file == "sketches/Fam_Gen_sp_GCA_000000.1.sig.zip"


# --- genome -----------------------------------------------------------------------


def test_genome_stage_sketches_the_downloads(tmp_path, monkeypatch, fake_sourmash) -> None:
    from test_genome_stage import _OUTGROUP, _SELECTED, _fake_run_cmd

    from repgenr.stages import genome

    monkeypatch.setattr(genome, "preflight", lambda caps: {"datasets": "16.0"})
    monkeypatch.setattr(genome, "_check_disk", lambda *a, **k: None)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    ctx.manifest.upsert_many([*_SELECTED, _OUTGROUP])
    _fake_run_cmd(monkeypatch)
    genome.run(ctx, genome.GenomeParams())
    assert len(_sketch_names(ctx.workdir)) == 3
    assert _sketch_names(ctx.workdir) == _genome_names(ctx.workdir)
    record = ctx.config.stages["genome"]
    assert record.params["sketches"]["written"] == 3
    assert record.tool_versions == {"datasets": "16.0", "sourmash": "4.9.4"}

    fake_sourmash.calls.clear()
    genome.run(ctx, genome.GenomeParams(sketch=False))
    assert fake_sourmash.calls == []
    assert ctx.config.stages["genome"].params["sketches"]["state"] == "off"


def test_genome_stage_without_sourmash_writes_no_sketches(tmp_path, monkeypatch) -> None:
    from test_genome_stage import _OUTGROUP, _SELECTED, _fake_run_cmd

    from repgenr.stages import genome

    monkeypatch.setattr(genome, "preflight", lambda caps: {"datasets": "16.0"})
    monkeypatch.setattr(genome, "_check_disk", lambda *a, **k: None)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    ctx.manifest.upsert_many([*_SELECTED, _OUTGROUP])
    _fake_run_cmd(monkeypatch)
    genome.run(ctx, genome.GenomeParams())
    assert not (ctx.workdir / "sketches").exists()
    assert ctx.config.stages["genome"].params["sketches"]["state"] == "unavailable"


# --- vgenome ----------------------------------------------------------------------


def test_vgenome_ncbi_virus_sketches(tmp_path, fake_sourmash) -> None:
    from test_ncbi_virus import _fake_records

    from repgenr.stages.vgenome import VgenomeParams
    from repgenr.stages.vgenome import run as vgenome_run
    from repgenr.viral.ncbi_virus import write_records

    wd = tmp_path / "wd"
    dl = wd / "virus_download_wd"
    dl.mkdir(parents=True)
    recs = _fake_records()
    (dl / "download.fa").write_text("".join(f">{r.accession} d\nACGTACGT\n" for r in recs))
    write_records(dl / "virus_records.json", recs)
    ctx = WorkdirContext(wd, create=True)
    vgenome_run(ctx, VgenomeParams(target_genus="lentivirus", length_all=True, no_outgroup=True))
    assert len(_sketch_names(wd)) == 2 and _sketch_names(wd) == _genome_names(wd)
    assert ctx.config.stages["vgenome"].params["sketches"]["written"] == 2

    vgenome_run(
        ctx,
        VgenomeParams(target_genus="lentivirus", length_all=True, no_outgroup=True, sketch=False),
    )
    assert ctx.config.stages["vgenome"].params["sketches"]["state"] == "off"


def test_vgenome_bvbrc_sketches(tmp_path, fake_sourmash) -> None:
    from test_bvbrc_select import _FASTA, _write_metadata

    from repgenr.stages.vgenome import VgenomeParams
    from repgenr.viral import bvbrc

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    dl = ctx.workdir / "virus_download_wd"
    dl.mkdir(parents=True)
    fasta = dl / "download.fa"
    fasta.write_text(_FASTA)
    base_tsv, ncbi_tsv = _write_metadata(dl)
    params = VgenomeParams(target_genus="mastadenovirus", no_outgroup=True, length_range="250-350")
    bvbrc.run_select(ctx, params, fasta, base_tsv, ncbi_tsv, _LOG)
    assert _sketch_names(ctx.workdir) == {"acc1", "acc2"}
    assert ctx.config.stages["vgenome"].params["sketches"]["written"] == 2


# --- assemble ---------------------------------------------------------------------


@pytest.fixture
def fake_assembler():
    register_fake_assembler()
    yield
    unregister_fake_assembler()


def test_assemble_append_sketches_only_the_new_genomes(tmp_path, fake_assembler, fake_sourmash):
    from test_assemble_stage import _gtdb_workdir

    from repgenr.stages.assemble import AssembleParams
    from repgenr.stages.assemble import run as assemble_run

    wd = tmp_path / "wd"
    ctx = _gtdb_workdir(wd)
    write_reads(wd / READS_TSV, [read_row(tmp_path, "SRR1")])
    assemble_run(ctx, AssembleParams(assembler="fakeasm", append=True))
    assert fake_sourmash.calls == ["Francisellaceae_Francisella_tularensis_SRR1"]
    assert _sketch_names(wd) == {"Francisellaceae_Francisella_tularensis_SRR1"}
    assert ctx.config.stages["assemble"].params["sketches"]["written"] == 1


def test_assemble_sketches_and_a_total_failure_clears_them(
    tmp_path, fake_assembler, fake_sourmash, monkeypatch
) -> None:
    from repgenr.stages import assemble as stage
    from repgenr.stages.assemble import AssembleParams
    from repgenr.stages.assemble import run as assemble_run

    wd = tmp_path / "wd"
    ctx = WorkdirContext(wd, create=True)
    write_reads(wd / READS_TSV, [read_row(tmp_path, "SRR1"), read_row(tmp_path, "SRR2")])
    assemble_run(ctx, AssembleParams(assembler="fakeasm", threads=2))
    assert len(_sketch_names(wd)) == 2

    # Every run of the next call is excused: the genome set and its sketches go.
    stage._clear_genome_set(ctx, _LOG)
    assert not (wd / "sketches").exists()


# --- the sketch command and status ------------------------------------------------


def _ingested(tmp_path: Path) -> Path:
    src = _genomes(tmp_path / "src")
    wd = tmp_path / "wd"
    result = _runner.invoke(
        app, ["ingest", "-wd", str(wd), "--genomes-dir", str(src), "--no-sketch"]
    )
    assert result.exit_code == 0, result.output
    return wd


def test_sketch_command_writes_missing_and_reports(tmp_path, fake_sourmash) -> None:
    wd = _ingested(tmp_path)
    result = _runner.invoke(app, ["sketch", "-wd", str(wd), "-t", "2"])
    assert result.exit_code == 0, result.output
    assert "Sketches: 0 present, 3 written, 0 stale replaced, 0 removed" in _log(wd)
    record = Config.load(wd).stages["sketch"]
    assert record.tool == "sourmash" and record.params["written"] == 3
    assert record.fingerprint

    status = _runner.invoke(app, ["status", "-wd", str(wd)])
    assert "sketches: 3/3" in status.output
    assert "    sketch [sourmash]" in status.output

    # A repeat skips; a sketch removed by hand reruns the command for that one.
    fake_sourmash.calls.clear()
    assert _runner.invoke(app, ["sketch", "-wd", str(wd)]).exit_code == 0
    assert fake_sourmash.calls == []
    (wd / "sketches" / "Fam_Gen_sp_GCA_000001.1.sig.zip").unlink()
    assert "sketches: 2/3" in _runner.invoke(app, ["status", "-wd", str(wd)]).output
    assert _runner.invoke(app, ["sketch", "-wd", str(wd)]).exit_code == 0
    assert fake_sourmash.calls == ["Fam_Gen_sp_GCA_000001.1"]

    # --force sketches every genome again.
    fake_sourmash.calls.clear()
    assert _runner.invoke(app, ["--force", "sketch", "-wd", str(wd)]).exit_code == 0
    assert sorted(fake_sourmash.calls) == sorted(_genome_names(wd))
    params = Config.load(wd).stages["sketch"].params
    assert (params["replaced"], params["forced"]) == (0, 3)


def test_sketch_command_refusals(tmp_path, monkeypatch) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    result = _runner.invoke(app, ["sketch", "-wd", str(empty)])
    assert result.exit_code == 3, result.output

    wd = _ingested(tmp_path)

    def missing(caps):
        raise MissingBinaryError("Required binary 'sourmash' not found on PATH.")

    from repgenr.stages import sketch as sketch_stage

    monkeypatch.setattr(sketch_stage, "preflight", missing)
    result = _runner.invoke(app, ["sketch", "-wd", str(wd)])
    assert result.exit_code == MissingBinaryError.exit_code
    assert "sketch" not in Config.load(wd).stages


def test_status_without_sketches_shows_no_count(tmp_path) -> None:
    wd = _ingested(tmp_path)
    assert "sketches:" not in _runner.invoke(app, ["status", "-wd", str(wd)]).output


def test_a_skipped_stage_with_sketch_names_the_sketch_command(tmp_path, fake_sourmash) -> None:
    src = _genomes(tmp_path / "src")
    wd = tmp_path / "wd"
    args = ["ingest", "-wd", str(wd), "--genomes-dir", str(src), "--sketch"]
    assert _runner.invoke(app, args).exit_code == 0
    assert _runner.invoke(app, args).exit_code == 0
    log = _log(wd)
    assert "already completed" in log and "A skipped stage writes no sketches" in log
