"""``repgenr reads``: select sequencing runs from ENA/SRA for the reads chain."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import (
    DEFAULT_THREADS,
    HELP_TARGET_FAMILY,
    HELP_TARGET_GENUS,
    HELP_TARGET_SPECIES,
    HELP_THREADS,
    _assembler_help,
    _parse_key_values,
    _run,
    app,
)


@app.command()
def reads(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    target_family: str | None = typer.Option(
        None, "-tf", "--target-family", help=HELP_TARGET_FAMILY
    ),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus", help=HELP_TARGET_GENUS),
    target_species: str | None = typer.Option(
        None, "-ts", "--target-species", help=HELP_TARGET_SPECIES
    ),
    accession: list[str] = typer.Option(
        [],
        "--accession",
        help="A run (SRR/ERR/DRR), sample (SAMN.., SRS..) or study (PRJNA.., SRP..) "
        "accession to include (repeatable).",
    ),
    accession_file: Path | None = typer.Option(
        None, "--accession-file", help="File of accessions, one per line (# comments allowed)."
    ),
    platform: str = typer.Option(
        "any", "--platform", help="Keep runs of one platform: any, illumina, ont or pacbio."
    ),
    max_runs: int | None = typer.Option(
        None, "--max-runs", min=1, help="Keep at most N runs, the largest by bases."
    ),
    min_bases: int = typer.Option(
        0, "--min-bases", min=0, help="Drop runs with fewer sequenced bases than this."
    ),
    one_per_sample: bool = typer.Option(
        True,
        "--one-per-sample/--all-runs",
        help="Keep the best run of each sample (long reads before short, then bases), "
        "or every run.",
    ),
) -> None:
    """Select sequencing runs from ENA/SRA by taxon or accession (writes reads.tsv)."""
    from .param_builders import reads_params

    def build():
        return reads_params(
            target_family=target_family,
            target_genus=target_genus,
            target_species=target_species,
            accessions=list(accession),
            accession_file=None if accession_file is None else str(accession_file),
            platform=platform,
            max_runs=max_runs,
            min_bases=min_bases,
            one_per_sample=one_per_sample,
        )

    _run("reads", workdir, build, create=True)


@app.command()
def assemble(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    assembler: str = typer.Option("auto", "--assembler", help=_assembler_help()),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS),
    jobs: int | None = typer.Option(
        None,
        "--jobs",
        min=1,
        help="Runs assembled concurrently, threads split across them (default 2, or 1 when a "
        "long-read run is pending).",
    ),
    memory_gb: int = typer.Option(
        16, "--memory-gb", min=1, help="Memory hint per assembly, in GB, for tools that cap RAM."
    ),
    min_contig_length: int = typer.Option(
        500, "--min-contig-length", min=0, help="Drop contigs shorter than this many bases."
    ),
    outgroup: Path | None = typer.Option(
        None, "--outgroup", help="A FASTA file to set aside as the outgroup for rooting."
    ),
    keep_reads: bool = typer.Option(
        False, "--keep-reads", help="Keep the downloaded FASTQ files after assembling."
    ),
    keep_files: bool = typer.Option(
        False, "--keep-files", help="Keep each run's assembler scratch directory."
    ),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Assembler tuning as key=value (repeatable), e.g. mode=nano-raw."
    ),
) -> None:
    """Fetch and assemble the selected runs; write genomes/ and selection.tsv."""
    from .param_builders import assemble_params

    def build():
        return assemble_params(
            assembler=assembler,
            threads=threads,
            jobs=jobs,
            memory_gb=memory_gb,
            min_contig_length=min_contig_length,
            outgroup=None if outgroup is None else str(outgroup),
            keep_reads=keep_reads,
            keep_files=keep_files,
            extra=_parse_key_values(tool_arg, "--tool-arg"),
        )

    _run("assemble", workdir, build)
