"""Contracts for the reads-to-assembly entry path: reads.tsv, assembly_stats.tsv,
excused_runs.tsv, and the completeness guard honouring excused runs."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    MISSING_ACCESSIONS_TXT,
    READS_TSV,
    SELECTION_TSV,
    AssemblyStatsRow,
    ExcusedRun,
    ReadRow,
    SelectionRow,
    read_assembly_stats,
    read_excused_runs,
    read_reads,
    write_assembly_stats,
    write_excused_runs,
    write_reads,
    write_selection,
)
from repgenr.core.errors import WorkdirError
from repgenr.core.integrity import check_genome_completeness, excused_accessions

_LOG = logging.getLogger("test")


def _row(run: str = "SRR1", **over) -> ReadRow:
    base = dict(
        run_accession=run,
        biosample="SAMN1",
        bioproject="PRJNA1",
        organism="Francisella tularensis",
        taxid="263",
        family="Francisellaceae",
        genus="Francisella",
        species="tularensis",
        platform="ILLUMINA",
        instrument_model="Illumina MiSeq",
        layout="PAIRED",
        bases=1_000_000,
        read_count=5000,
        fastq_urls=("ftp://a/x_1.fastq.gz", "ftp://a/x_2.fastq.gz"),
        fastq_md5=("a" * 32, "b" * 32),
        fastq_bytes=(100, 200),
    )
    base.update(over)
    return ReadRow(**base)


def test_reads_table_round_trips(tmp_path: Path) -> None:
    rows = [_row(), _row("SRR2", fastq_urls=(), fastq_md5=(), fastq_bytes=())]
    write_reads(tmp_path / READS_TSV, rows)
    back = read_reads(tmp_path / READS_TSV)
    assert back == rows
    assert back[1].fastq_urls == ()  # a run without an ENA mirror keeps empty lists


def test_assembly_stats_round_trips(tmp_path: Path) -> None:
    row = AssemblyStatsRow(
        run_accession="SRR1",
        filename="Fam_Gen_sp_SRR1.fasta",
        assembler="skesa",
        n_contigs=12,
        total_length=1_900_000,
        n50=250_000,
        largest_contig=400_000,
        est_coverage=52.6,
        completeness=99.1,
        contamination=0.4,
        ncbi_taxonomy="Francisellaceae;Francisella;tularensis",
        gtdb_taxonomy="",
        label_source="metadata",
        taxonomy_flag="",
    )
    write_assembly_stats(tmp_path / ASSEMBLY_STATS_TSV, [row])
    assert read_assembly_stats(tmp_path / ASSEMBLY_STATS_TSV) == [row]


def test_excused_runs_round_trip(tmp_path: Path) -> None:
    rows = [ExcusedRun("SRR9", "fetch", "no_fastq_mirror"), ExcusedRun("SRR8", "qc", "qc_failed")]
    write_excused_runs(tmp_path / EXCUSED_RUNS_TSV, rows)
    assert read_excused_runs(tmp_path / EXCUSED_RUNS_TSV) == rows


def test_excused_accessions_unions_both_sources(tmp_path: Path) -> None:
    (tmp_path / MISSING_ACCESSIONS_TXT).write_text("GCF_1\n\nGCF_2\n", encoding="utf-8")
    write_excused_runs(tmp_path / EXCUSED_RUNS_TSV, [ExcusedRun("SRR9", "fetch", "x")])
    assert excused_accessions(tmp_path) == {"GCF_1", "GCF_2", "SRR9"}
    assert excused_accessions(tmp_path / "nowhere") == set()


def test_completeness_guard_excuses_failed_runs(tmp_path: Path) -> None:
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "Fam_Gen_sp_SRR1.fasta").write_text(">a\nACGT\n", encoding="utf-8")
    write_selection(
        tmp_path / SELECTION_TSV,
        [
            SelectionRow("SRR1", "Fam", "Gen", "sp", False, "Fam_Gen_sp_SRR1.fasta"),
            SelectionRow("SRR2", "Fam", "Gen", "sp", False, "Fam_Gen_sp_SRR2.fasta"),
        ],
    )
    with pytest.raises(WorkdirError):
        check_genome_completeness(genomes, tmp_path, logger=_LOG)
    write_excused_runs(tmp_path / EXCUSED_RUNS_TSV, [ExcusedRun("SRR2", "assemble", "failed")])
    assert check_genome_completeness(genomes, tmp_path, logger=_LOG) == []
