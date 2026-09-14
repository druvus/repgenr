"""The assemble stage: fetch reads, assemble, label, and write the genome contract."""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
from pathlib import Path

import pytest

from repgenr.assemblers.base import Assembler, AssemblyResult, registry
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    READS_TSV,
    SELECTION_TSV,
    ReadRow,
    read_assembly_stats,
    read_excused_runs,
    read_selection,
    write_reads,
)
from repgenr.core.errors import ToolExecutionError, WorkdirError
from repgenr.core.integrity import check_genome_completeness
from repgenr.core.plugins import ToolCapabilities
from repgenr.stages.assemble import AssembleParams, run

_LOG = logging.getLogger("test")
_GENOME = "ACGT" * 300  # 1200 bp, above the default contig floor


class _FakeAssembler(Assembler):
    """Writes one contig from the reads it was given; records every call."""

    capabilities = ToolCapabilities(name="fakeasm")
    read_types = frozenset({"ILLUMINA", "OXFORD_NANOPORE"})
    calls: list[str] = []
    fail_runs: frozenset[str] = frozenset()

    def preflight(self) -> dict[str, str]:
        return {"fakeasm": "1.0"}

    def assemble(self, reads, out_dir, params, logger) -> AssemblyResult:  # noqa: ANN001
        type(self).calls.append(reads.run_accession)
        if reads.run_accession in type(self).fail_runs:
            raise ToolExecutionError(["fakeasm"], 1, "boom")
        assert all(f.exists() for f in reads.files)
        out_dir.mkdir(parents=True, exist_ok=True)
        contigs = out_dir / "raw.fa"
        contigs.write_text(f">node1\n{_GENOME}\n>tiny\nACGT\n", encoding="utf-8")
        return AssemblyResult(contigs=contigs, tool_stats={"threads": params.threads})


@pytest.fixture
def fake_assembler():
    registry._load()
    registry.register("fakeasm", _FakeAssembler, replace=True)
    _FakeAssembler.calls = []
    _FakeAssembler.fail_runs = frozenset()
    yield
    registry._classes.pop("fakeasm", None)


def _fastq(path: Path) -> tuple[str, str, int]:
    """Write a tiny gzipped FASTQ; return (local url, md5, bytes)."""
    data = gzip.compress(b"@r1\nACGT\n+\nIIII\n")
    path.write_bytes(data)
    return str(path), hashlib.md5(data).hexdigest(), len(data)


def _row(tmp_path: Path, run: str, platform: str = "ILLUMINA", layout: str = "PAIRED", **over):
    files = [
        _fastq(tmp_path / f"{run}_{i}.fastq.gz") for i in ((1, 2) if layout == "PAIRED" else (1,))
    ]
    base = dict(
        run_accession=run,
        biosample=f"SAM{run}",
        bioproject="PRJ1",
        organism="Francisella tularensis",
        taxid="263",
        platform=platform,
        instrument_model="MiSeq",
        layout=layout,
        bases=1_200_000,
        read_count=1000,
        family="Francisellaceae",
        genus="Francisella",
        species="tularensis",
        fastq_urls=tuple(f[0] for f in files),
        fastq_md5=tuple(f[1] for f in files),
        fastq_bytes=tuple(f[2] for f in files),
    )
    base.update(over)
    return ReadRow(**base)


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


class _FakeClassifier:
    """Registered classifier returning a canned GTDB lineage per genome."""

    from repgenr.core.plugins import ToolCapabilities as _TC

    capabilities = _TC(name="fakecls")
    lineages: dict[str, str] = {}

    def preflight(self) -> dict[str, str]:
        return {"fakecls": "1.0"}

    def classify(self, genomes, out_dir, params, logger):  # noqa: ANN001
        from repgenr.classifiers.base import Classification, db_version

        return {
            g.name: Classification(
                taxonomy=type(self).lineages[g.name],
                rank="species",
                score=0.9,
                db_version=db_version(params.db),
            )
            for g in genomes
            if g.name in type(self).lineages
        }


@pytest.fixture
def fake_classifier():
    from repgenr.classifiers.base import Classifier
    from repgenr.classifiers.base import registry as cls_registry

    class Fake(_FakeClassifier, Classifier):
        pass

    cls_registry._load()
    cls_registry.register("fakecls", Fake, replace=True)
    _FakeClassifier.lineages = {}
    yield
    cls_registry._classes.pop("fakecls", None)


def _fake_checkm2(quality: dict[str, tuple[float, float]]):
    def run_checkm2(genomes, out_dir, *, db, threads, logger):
        return {g.name: quality[g.name] for g in genomes if g.name in quality}

    return run_checkm2


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
                "Francisellaceae_Francisella_tularensis_SRR1.fasta": (98.5, 0.4),
                "Francisellaceae_Francisella_tularensis_SRR2.fasta": (40.0, 15.0),
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
        "Francisellaceae_Francisella_tularensis_SRR1.fasta": (
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales;"
            "f__Francisellaceae;g__Francisella;s__Francisella tularensis_A"
        ),
        "Francisellaceae_Francisella_tularensis_SRR2.fasta": (
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
