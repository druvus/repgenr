"""The stateless reads-chain steps behind the Nextflow modules: assemble-run
(one run), genome-qc (a batch of assemblies) and reads-gather (the genome
contract from the per-run results)."""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assemble_fakes import (  # noqa: E402
    FakeAssembler,
    FakeClassifier,
    fake_checkm2,
    read_row,
    register_fake_assembler,
    register_fake_classifier,
    unregister_fake_assembler,
    unregister_fake_classifier,
)

from repgenr.core.contracts import (  # noqa: E402
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    SELECTION_TSV,
    read_assembly_stats,
    read_excused_runs,
    read_selection,
    write_reads,
)
from repgenr.core.errors import UserInputError, WorkdirError  # noqa: E402
from repgenr.stages.assemble_steps import (  # noqa: E402
    AssembleRunParams,
    GenomeQcParams,
    ReadsGatherParams,
    assemble_run,
    genome_qc,
    read_classification,
    read_quality,
    reads_gather,
)

_LOG = logging.getLogger("test")


@pytest.fixture
def fakes():
    register_fake_assembler()
    register_fake_classifier()
    yield
    unregister_fake_assembler()
    unregister_fake_classifier()


def _reads_tsv(tmp_path: Path, rows) -> Path:
    path = tmp_path / "reads.tsv"
    write_reads(path, rows)
    return path


# --- assemble-run ---------------------------------------------------------------


def test_assemble_run_writes_contigs_and_a_marker(tmp_path: Path, fakes) -> None:
    reads = _reads_tsv(tmp_path, [read_row(tmp_path, "SRR1"), read_row(tmp_path, "SRR2")])
    out = tmp_path / "SRR1"
    versions = tmp_path / "v.yml"
    assembled = assemble_run(
        AssembleRunParams(
            reads_tsv=reads,
            run="SRR1",
            out_dir=out,
            assembler="fakeasm",
            threads=3,
            versions_out=versions,
        ),
        _LOG,
    )
    assert assembled is True
    assert FakeAssembler.calls == ["SRR1"]  # only the named run
    text = (out / "contigs.fasta").read_text(encoding="utf-8")
    assert text.startswith(">SRR1_contig1\n") and "tiny" not in text
    marker = json.loads((out / "assembly.ok").read_text(encoding="utf-8"))
    assert marker["assembler"] == "fakeasm" and marker["stats"]["total_length"] == 1200
    assert marker["tool_stats"]["threads"] == 3
    assert "fakeasm: 1.0" in versions.read_text(encoding="utf-8")
    assert not (out / "scratch").exists() and not list(tmp_path.glob("**/asm"))


def test_assemble_run_records_an_excuse_instead_of_failing(tmp_path: Path, fakes) -> None:
    rows = [
        read_row(tmp_path, "NOMIRROR", fastq_urls=(), fastq_md5=(), fastq_bytes=()),
        read_row(tmp_path, "CRASH"),
        read_row(tmp_path, "IONT", platform="ION_TORRENT"),
    ]
    FakeAssembler.fail_runs = frozenset({"CRASH"})
    reads = _reads_tsv(tmp_path, rows)
    for run, step, reason in [
        ("NOMIRROR", "fetch", "no_fastq_mirror"),
        ("CRASH", "assemble", "assembly_failed"),
        ("IONT", "assemble", "unsupported_platform"),
    ]:
        out = tmp_path / run
        assembled = assemble_run(
            AssembleRunParams(reads_tsv=reads, run=run, out_dir=out, assembler="auto"), _LOG
        )
        assert assembled is False
        assert not (out / "assembly.ok").exists()
        excused = read_excused_runs(out / EXCUSED_RUNS_TSV)
        assert len(excused) == 1 and excused[0].step == step and reason in excused[0].reason


def test_assemble_run_rejects_an_unknown_run(tmp_path: Path, fakes) -> None:
    reads = _reads_tsv(tmp_path, [read_row(tmp_path, "SRR1")])
    with pytest.raises(UserInputError, match="SRR9"):
        assemble_run(AssembleRunParams(reads_tsv=reads, run="SRR9", out_dir=tmp_path / "o"), _LOG)


# --- genome-qc -------------------------------------------------------------------


def _assembled(tmp_path: Path, runs: list[str], fakes_reads: Path) -> Path:
    assemblies = tmp_path / "assemblies"
    for run in runs:
        assemble_run(
            AssembleRunParams(
                reads_tsv=fakes_reads, run=run, out_dir=assemblies / run, assembler="fakeasm"
            ),
            _LOG,
        )
    return assemblies


def test_genome_qc_scores_and_classifies_every_assembly(tmp_path, fakes, monkeypatch) -> None:
    from repgenr.stages import assemble as stage

    reads = _reads_tsv(tmp_path, [read_row(tmp_path, "SRR1"), read_row(tmp_path, "SRR2")])
    assemblies = _assembled(tmp_path, ["SRR1", "SRR2"], reads)
    monkeypatch.setattr(stage, "preflight_checkm2", lambda: {"checkm2": "1.1.0"})
    monkeypatch.setattr(
        stage, "run_checkm2", fake_checkm2({"SRR1.fasta": (98.5, 0.4), "SRR2.fasta": (40.0, 15.0)})
    )
    FakeClassifier.lineages = {
        "SRR1.fasta": "d__Bacteria;f__Francisellaceae;g__Francisella;s__Francisella tularensis",
    }
    sketch = tmp_path / "gtdb-rs226.k31.sig.zip"
    sketch.write_bytes(b"x")
    out = tmp_path / "qc"
    n = genome_qc(
        GenomeQcParams(
            assemblies_dir=assemblies,
            out_dir=out,
            checkm2_db="/db/checkm2",
            classifier="fakecls",
            gtdb_sketch=str(sketch),
            versions_out=tmp_path / "v.yml",
        ),
        _LOG,
    )
    assert n == 2
    assert read_quality(out / "quality.tsv") == {"SRR1": (98.5, 0.4), "SRR2": (40.0, 15.0)}
    classified = read_classification(out / "classification.tsv")
    assert list(classified) == ["SRR1"]
    assert classified["SRR1"].taxonomy.endswith("s__Francisella tularensis")
    assert classified["SRR1"].db_version == "gtdb-rs226.k31"
    versions = (tmp_path / "v.yml").read_text(encoding="utf-8")
    assert "checkm2: 1.1.0" in versions and "fakecls: 1.0" in versions


def test_genome_qc_needs_a_database(tmp_path, fakes, monkeypatch) -> None:
    monkeypatch.delenv("CHECKM2DB", raising=False)
    monkeypatch.delenv("REPGENR_GTDB_SKETCH", raising=False)
    reads = _reads_tsv(tmp_path, [read_row(tmp_path, "SRR1")])
    assemblies = _assembled(tmp_path, ["SRR1"], reads)
    with pytest.raises(UserInputError, match="checkm2-db"):
        genome_qc(GenomeQcParams(assemblies_dir=assemblies, out_dir=tmp_path / "qc"), _LOG)


# --- reads-gather -------------------------------------------------------------------


def test_reads_gather_writes_the_genome_contract(tmp_path: Path, fakes) -> None:
    rows = [
        read_row(tmp_path, "SRR1"),
        read_row(tmp_path, "SRR2"),
        read_row(tmp_path, "IONT", platform="ION_TORRENT"),
    ]
    reads = _reads_tsv(tmp_path, rows)
    assemblies = tmp_path / "assemblies"
    for run in ["SRR1", "SRR2", "IONT"]:
        assemble_run(
            AssembleRunParams(reads_tsv=reads, run=run, out_dir=assemblies / run, assembler="auto"),
            _LOG,
        )
    out = tmp_path / "out"
    n = reads_gather(
        ReadsGatherParams(reads_tsv=reads, assemblies_dir=assemblies, out_dir=out), _LOG
    )
    assert n == 2
    names = sorted(p.name for p in (out / "genomes").iterdir())
    assert names == [
        "Francisellaceae_Francisella_tularensis_SRR1.fasta",
        "Francisellaceae_Francisella_tularensis_SRR2.fasta",
    ]
    assert (out / "genomes" / names[0]).read_text(encoding="utf-8").startswith(">SRR1_contig1\n")
    selection = read_selection(out / SELECTION_TSV)
    assert [r.accession for r in selection] == ["SRR1", "SRR2"]
    assert selection[0].completeness is None
    stats = {s.run_accession: s for s in read_assembly_stats(out / ASSEMBLY_STATS_TSV)}
    assert stats["SRR1"].assembler == "fakeasm" and stats["SRR1"].n50 == 1200
    assert stats["SRR1"].label_source == "metadata" and stats["SRR1"].est_coverage == 1000.0
    excused = read_excused_runs(out / EXCUSED_RUNS_TSV)
    assert [(e.run_accession, e.reason) for e in excused] == [("IONT", "unsupported_platform")]
    assert (out / "outgroup_accession.txt").read_text(encoding="utf-8") == ""


def test_reads_gather_applies_the_quality_gate_and_gtdb_names(tmp_path, fakes, monkeypatch):
    from repgenr.stages import assemble as stage

    rows = [read_row(tmp_path, "SRR1"), read_row(tmp_path, "SRR2"), read_row(tmp_path, "SRR3")]
    reads = _reads_tsv(tmp_path, rows)
    assemblies = _assembled(tmp_path, ["SRR1", "SRR2", "SRR3"], reads)
    monkeypatch.setattr(stage, "preflight_checkm2", lambda: {"checkm2": "1.1.0"})
    monkeypatch.setattr(
        stage,
        "run_checkm2",
        fake_checkm2(
            {"SRR1.fasta": (98.5, 0.4), "SRR2.fasta": (40.0, 15.0), "SRR3.fasta": (97.0, 1.0)}
        ),
    )
    FakeClassifier.lineages = {
        "SRR1.fasta": "d__Bacteria;f__Francisellaceae;g__Francisella;s__Francisella novicida",
        "SRR3.fasta": "d__Bacteria;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
    }
    sketch = tmp_path / "gtdb-rs226.k31.sig.zip"
    sketch.write_bytes(b"x")
    qc = tmp_path / "qc"
    genome_qc(
        GenomeQcParams(
            assemblies_dir=assemblies,
            out_dir=qc,
            checkm2_db="/db",
            classifier="fakecls",
            gtdb_sketch=str(sketch),
        ),
        _LOG,
    )
    out = tmp_path / "out"
    n = reads_gather(
        ReadsGatherParams(reads_tsv=reads, assemblies_dir=assemblies, out_dir=out, qc_dir=qc),
        _LOG,
    )
    assert n == 2
    names = sorted(p.name for p in (out / "genomes").iterdir())
    # SRR1 takes the GTDB species; SRR3 keeps the submitted name and is flagged.
    assert names == [
        "Francisellaceae_Francisella_novicida_SRR1.fasta",
        "Francisellaceae_Francisella_tularensis_SRR3.fasta",
    ]
    selection = {r.accession: r for r in read_selection(out / SELECTION_TSV)}
    assert selection["SRR1"].completeness == 98.5 and selection["SRR1"].species == "novicida"
    stats = {s.run_accession: s for s in read_assembly_stats(out / ASSEMBLY_STATS_TSV)}
    assert stats["SRR1"].label_source == "classifier" and stats["SRR1"].taxonomy_flag == ""
    assert stats["SRR3"].taxonomy_flag == "classifier_disagrees"
    assert stats["SRR3"].gtdb_taxonomy.endswith("s__Escherichia coli")
    excused = {e.run_accession: e for e in read_excused_runs(out / EXCUSED_RUNS_TSV)}
    assert excused["SRR2"].step == "qc" and "qc_failed" in excused["SRR2"].reason


def test_reads_gather_fails_when_nothing_is_accepted(tmp_path: Path, fakes) -> None:
    reads = _reads_tsv(tmp_path, [read_row(tmp_path, "IONT", platform="ION_TORRENT")])
    assemblies = tmp_path / "assemblies"
    assemble_run(AssembleRunParams(reads_tsv=reads, run="IONT", out_dir=assemblies / "IONT"), _LOG)
    out = tmp_path / "out"
    with pytest.raises(WorkdirError, match="None of the 1 runs"):
        reads_gather(
            ReadsGatherParams(reads_tsv=reads, assemblies_dir=assemblies, out_dir=out), _LOG
        )
    assert read_excused_runs(out / EXCUSED_RUNS_TSV)[0].run_accession == "IONT"
