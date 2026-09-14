"""The assemble stage: fetch reads, assemble, label, and write the genome contract."""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    READS_TSV,
    SELECTION_TSV,
    ReadRow,
    SelectionRow,
    read_assembly_stats,
    read_excused_runs,
    read_selection,
    write_reads,
    write_selection,
)
from repgenr.core.errors import WorkdirError
from repgenr.core.integrity import check_genome_completeness
from repgenr.stages.assemble import AssembleParams, run

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assemble_fakes import FakeAssembler as _FakeAssembler  # noqa: E402
from assemble_fakes import FakeClassifier as _FakeClassifier  # noqa: E402
from assemble_fakes import fake_checkm2 as _fake_checkm2  # noqa: E402
from assemble_fakes import read_row as _row  # noqa: E402
from assemble_fakes import (  # noqa: E402
    register_fake_assembler,
    register_fake_classifier,
    unregister_fake_assembler,
    unregister_fake_classifier,
)

_LOG = logging.getLogger("test")


@pytest.fixture
def fake_assembler():
    register_fake_assembler()
    yield
    unregister_fake_assembler()


def _prepare(workdir: Path, rows: list[ReadRow]) -> WorkdirContext:
    ctx = WorkdirContext(workdir, create=True)
    write_reads(workdir / READS_TSV, rows)
    return ctx


def test_assemble_writes_the_genome_contract(workdir: Path, tmp_path: Path, fake_assembler) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    n = run(ctx, AssembleParams(assembler="fakeasm", threads=4, jobs=2))
    assert n == 2
    names = sorted(p.name for p in ctx.genomes_dir.iterdir())
    assert names == [
        "Francisellaceae_Francisella_tularensis_SRR1.fasta",
        "Francisellaceae_Francisella_tularensis_SRR2.fasta",
    ]
    text = (ctx.genomes_dir / names[0]).read_text(encoding="utf-8")
    assert text.startswith(">SRR1_contig1\n") and "tiny" not in text  # filtered and renamed
    rows = read_selection(workdir / SELECTION_TSV)
    assert {r.accession for r in rows} == {"SRR1", "SRR2"}
    assert {r.filename for r in rows} == set(names)
    stats = {s.run_accession: s for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)}
    assert stats["SRR1"].n_contigs == 1 and stats["SRR1"].total_length == 1200
    assert stats["SRR1"].assembler == "fakeasm" and stats["SRR1"].est_coverage == 1000.0
    assert stats["SRR1"].ncbi_taxonomy == "Francisellaceae;Francisella;tularensis"
    manifest = {g.accession: g for g in ctx.manifest.all_genomes()}
    assert manifest["SRR1"].source == "sra" and manifest["SRR1"].species == "tularensis"
    record = ctx.config.stages["assemble"]
    assert record.params["n_assembled"] == 2 and record.params["n_excused"] == 0
    assert record.tool_versions == {"fakeasm": "1.0"}
    assert check_genome_completeness(ctx.genomes_dir, workdir, logger=_LOG) == []
    assert not (ctx.scratch_dir / "assemble" / "SRR1").exists()  # reads and scratch removed


def test_a_second_run_skips_finished_assemblies(
    workdir: Path, tmp_path: Path, fake_assembler
) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm"))
    assert _FakeAssembler.calls == ["SRR1"]
    marker = workdir / "assemblies" / "SRR1" / "assembly.ok"
    assert json.loads(marker.read_text(encoding="utf-8"))["assembler"] == "fakeasm"
    run(ctx, AssembleParams(assembler="fakeasm"))
    assert _FakeAssembler.calls == ["SRR1"]  # not assembled again
    assert len(read_selection(workdir / SELECTION_TSV)) == 1


def test_failures_are_excused_and_the_rest_proceed(
    workdir: Path, tmp_path: Path, fake_assembler
) -> None:
    rows = [
        _row(tmp_path, "SRR1"),
        _row(tmp_path, "NOMIRROR", fastq_urls=(), fastq_md5=(), fastq_bytes=()),
        _row(tmp_path, "BADSUM", fastq_md5=("0" * 32, "0" * 32)),
        _row(tmp_path, "CRASH"),
        _row(tmp_path, "IONT", platform="ION_TORRENT"),
    ]
    _FakeAssembler.fail_runs = frozenset({"CRASH"})
    ctx = _prepare(workdir, rows)
    n = run(ctx, AssembleParams(assembler="auto"))
    assert n == 1
    excused = {
        e.run_accession: (e.step, e.reason) for e in read_excused_runs(workdir / EXCUSED_RUNS_TSV)
    }
    assert excused["NOMIRROR"] == ("fetch", "no_fastq_mirror")
    assert excused["BADSUM"][0] == "fetch" and "download_failed" in excused["BADSUM"][1]
    assert excused["CRASH"][0] == "assemble" and "assembly_failed" in excused["CRASH"][1]
    assert excused["IONT"] == ("assemble", "unsupported_platform")
    assert [r.accession for r in read_selection(workdir / SELECTION_TSV)] == ["SRR1"]
    assert check_genome_completeness(ctx.genomes_dir, workdir, logger=_LOG) == []
    assert ctx.config.stages["assemble"].params["n_excused"] == 4


def test_outgroup_fasta_is_staged(workdir: Path, tmp_path: Path, fake_assembler) -> None:
    og = tmp_path / "Fam_Gen_sp_GCF_000009.1.fasta"
    og.write_text(">og\nACGT\n", encoding="utf-8")
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", outgroup=str(og)))
    assert [p.name for p in ctx.outgroup_dir.iterdir()] == [og.name]
    assert (workdir / "outgroup_accession.txt").read_text().strip() == "GCF_000009.1"
    rows = read_selection(workdir / SELECTION_TSV)
    assert {r.accession for r in rows if r.is_outgroup} == {"GCF_000009.1"}


def test_missing_reads_table_is_a_workdir_error(workdir: Path, fake_assembler) -> None:
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(WorkdirError, match="reads"):
        run(ctx, AssembleParams(assembler="fakeasm"))


def test_threads_are_split_across_jobs(workdir: Path, tmp_path: Path, fake_assembler) -> None:
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    run(ctx, AssembleParams(assembler="fakeasm", threads=8, jobs=2))
    marker = json.loads(
        (workdir / "assemblies" / "SRR1" / "assembly.ok").read_text(encoding="utf-8")
    )
    assert marker["tool_stats"]["threads"] == 4


def test_excused_runs_are_written_even_when_nothing_assembles(workdir, tmp_path, fake_assembler):
    _FakeAssembler.fail_runs = frozenset({"SRR1"})
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    with pytest.raises(WorkdirError, match="None of the 1 runs"):
        run(ctx, AssembleParams(assembler="fakeasm"))
    assert read_excused_runs(workdir / EXCUSED_RUNS_TSV)[0].step == "assemble"


def test_jobs_default_is_one_when_long_reads_are_pending(workdir, tmp_path, fake_assembler, caplog):
    """Memory bounds concurrent assemblies; a long-read run gets the machine alone."""
    rows = [
        _row(tmp_path, "SRR1"),
        _row(
            tmp_path,
            "ONT1",
            platform="OXFORD_NANOPORE",
            layout="SINGLE",
            instrument_model="GridION",
        ),
    ]
    ctx = _prepare(workdir, rows)
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm", threads=8))
    assert any("with 1 concurrent job" in r.message for r in caplog.records)
    assert ctx.config.stages["assemble"].params["jobs"] is None  # the request, not the resolution


def test_jobs_default_is_two_for_short_reads(workdir, tmp_path, fake_assembler, caplog):
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm", threads=8))
    assert any("with 2 concurrent jobs" in r.message for r in caplog.records)


# --- quality and classification -----------------------------------------------------


@pytest.fixture
def fake_classifier():
    register_fake_classifier()
    yield
    unregister_fake_classifier()


def test_checkm2_quality_gates_and_feeds_the_selection(
    workdir, tmp_path, fake_assembler, monkeypatch
) -> None:
    from repgenr.stages import assemble as stage

    monkeypatch.setattr(stage, "preflight_checkm2", lambda: {"checkm2": "1.1.0"})
    monkeypatch.setattr(
        stage,
        "run_checkm2",
        _fake_checkm2(
            {
                "SRR1.fasta": (98.5, 0.4),
                "SRR2.fasta": (40.0, 15.0),
            }
        ),
    )
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    n = run(ctx, AssembleParams(assembler="fakeasm", checkm2_db=str(tmp_path / "db")))
    assert n == 1
    rows = read_selection(workdir / SELECTION_TSV)
    assert [(r.accession, r.completeness, r.contamination) for r in rows] == [("SRR1", 98.5, 0.4)]
    assert ctx.manifest.quality()["Francisellaceae_Francisella_tularensis_SRR1.fasta"] == (
        98.5,
        0.4,
    )
    excused = {e.run_accession: e for e in read_excused_runs(workdir / EXCUSED_RUNS_TSV)}
    assert excused["SRR2"].step == "qc" and "qc_failed" in excused["SRR2"].reason
    stats = {s.run_accession: s for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)}
    assert stats["SRR1"].completeness == 98.5 and "SRR2" not in stats
    assert ctx.config.stages["assemble"].params["checkm2_db"] == str(tmp_path / "db")


def test_classifier_agreement_names_the_genome_with_gtdb_tokens(
    workdir, tmp_path, fake_assembler, fake_classifier
) -> None:
    _FakeClassifier.lineages = {
        "SRR1.fasta": (
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales;"
            "f__Francisellaceae;g__Francisella;s__Francisella tularensis_A"
        ),
        "SRR2.fasta": (
            "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;"
            "g__Bacillus;s__Bacillus subtilis"
        ),
    }
    ctx = _prepare(workdir, [_row(tmp_path, "SRR1"), _row(tmp_path, "SRR2")])
    run(
        ctx,
        AssembleParams(
            assembler="fakeasm",
            classifier="fakecls",
            gtdb_sketch=str(tmp_path / "gtdb-rs226-reps.k31-sc10k.sig.zip"),
            gtdb_lineages=str(tmp_path / "lineages.csv"),
        ),
    )
    rows = {r.accession: r for r in read_selection(workdir / SELECTION_TSV)}
    # agreement at genus: the GTDB tokens name the file
    assert rows["SRR1"].filename == "Francisellaceae_Francisella_tularensis-A_SRR1.fasta"
    assert rows["SRR1"].species == "tularensis-A"
    assert (ctx.genomes_dir / rows["SRR1"].filename).exists()
    # disagreement: NCBI tokens stay and the genome is flagged
    assert rows["SRR2"].filename == "Francisellaceae_Francisella_tularensis_SRR2.fasta"
    stats = {s.run_accession: s for s in read_assembly_stats(workdir / ASSEMBLY_STATS_TSV)}
    assert stats["SRR1"].label_source == "classifier" and stats["SRR1"].taxonomy_flag == ""
    assert stats["SRR2"].label_source == "metadata"
    assert stats["SRR2"].taxonomy_flag == "classifier_disagrees"
    assert stats["SRR2"].gtdb_taxonomy.endswith("s__Bacillus subtilis")
    record = ctx.config.stages["assemble"]
    assert record.params["n_disagree"] == 1
    assert record.tool_versions["gtdb_sketch"] == "gtdb-rs226-reps.k31-sc10k"


# --- appending into an existing working directory -----------------------------------


def _gtdb_workdir(workdir: Path) -> WorkdirContext:
    """A workdir as metadata + genome leave it: two GTDB genomes and an outgroup."""
    from repgenr.core.manifest import GenomeRecord

    ctx = WorkdirContext(workdir, create=True)
    ctx.genomes_dir.mkdir(parents=True)
    rows = [
        SelectionRow("GCF_1", "Fam", "Gen", "sp", False, "Fam_Gen_sp_GCF_1.fasta", 99.0, 0.1),
        SelectionRow("GCF_2", "Fam", "Gen", "sp", False, "Fam_Gen_sp_GCF_2.fasta", 98.0, 0.2),
        SelectionRow("GCF_9", "Fam", "Out", "og", True, "Fam_Out_og_GCF_9.fasta"),
    ]
    for r in rows[:2]:
        (ctx.genomes_dir / r.filename).write_text(">g\nACGT\n", encoding="utf-8")
    ctx.outgroup_dir.mkdir(parents=True)
    (ctx.outgroup_dir / rows[2].filename).write_text(">o\nACGT\n", encoding="utf-8")
    (workdir / "outgroup_accession.txt").write_text("GCF_9\n", encoding="utf-8")
    write_selection(workdir / SELECTION_TSV, rows)
    ctx.manifest.replace_genomes(
        [
            GenomeRecord(
                r.accession,
                r.filename,
                "gtdb",
                r.family,
                r.genus,
                r.species,
                r.is_outgroup,
                completeness=r.completeness,
                contamination=r.contamination,
            )
            for r in rows
        ]
    )
    return ctx


def test_append_adds_assemblies_to_a_gtdb_workdir(workdir, tmp_path, fake_assembler) -> None:
    ctx = _gtdb_workdir(workdir)
    write_reads(workdir / READS_TSV, [_row(tmp_path, "SRR1")])
    n = run(ctx, AssembleParams(assembler="fakeasm", append=True))
    assert n == 1
    rows = {r.accession: r for r in read_selection(workdir / SELECTION_TSV)}
    assert set(rows) == {"GCF_1", "GCF_2", "GCF_9", "SRR1"}
    assert rows["GCF_1"].completeness == 99.0  # existing rows untouched
    assert rows["GCF_9"].is_outgroup and rows["SRR1"].filename.endswith("_SRR1.fasta")
    names = {p.name for p in ctx.genomes_dir.iterdir()}
    assert names == {"Fam_Gen_sp_GCF_1.fasta", "Fam_Gen_sp_GCF_2.fasta", rows["SRR1"].filename}
    assert (workdir / "outgroup_accession.txt").read_text().strip() == "GCF_9"
    sources = {g.accession: g.source for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert sources == {"GCF_1": "gtdb", "GCF_2": "gtdb", "GCF_9": "gtdb", "SRR1": "sra"}
    assert check_genome_completeness(ctx.genomes_dir, workdir, logger=_LOG) == []
    assert ctx.config.stages["assemble"].params["append"] is True


def test_append_replaces_its_own_earlier_rows(workdir, tmp_path, fake_assembler) -> None:
    ctx = _gtdb_workdir(workdir)
    write_reads(workdir / READS_TSV, [_row(tmp_path, "SRR1")])
    run(ctx, AssembleParams(assembler="fakeasm", append=True))
    write_reads(workdir / READS_TSV, [_row(tmp_path, "SRR1", species="holarctica")])
    for marker in (workdir / "assemblies").rglob("assembly.ok"):
        marker.unlink()  # force a re-assembly under the new label
    run(ctx, AssembleParams(assembler="fakeasm", append=True))
    rows = read_selection(workdir / SELECTION_TSV)
    assert [r.accession for r in rows].count("SRR1") == 1
    assert next(r for r in rows if r.accession == "SRR1").species == "holarctica"
    assert not (ctx.genomes_dir / "Francisellaceae_Francisella_tularensis_SRR1.fasta").exists()


def test_append_needs_an_existing_selection(workdir, tmp_path, fake_assembler) -> None:
    from repgenr.core.errors import UserInputError

    ctx = _prepare(workdir, [_row(tmp_path, "SRR1")])
    with pytest.raises(UserInputError, match="--append"):
        run(ctx, AssembleParams(assembler="fakeasm", append=True))
