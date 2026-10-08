"""``repgenr reads``: select sequencing runs from ENA/SRA for the reads chain."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import (
    DEFAULT_THREADS,
    HELP_THREADS,
    HELP_WORKDIR,
    HELP_WORKDIR_CREATED,
    PANEL_ENTRY,
    _assembler_help,
    _classifier_help,
    _parse_key_values,
    _polisher_help,
    _run,
    app,
    resolve_threads,
)

_HELP_TAXON = (
    "Restrict the selection to this {}. Only the most specific of -tf/-tg/-ts "
    "given is used; they are not combined."
)


@app.command(rich_help_panel=PANEL_ENTRY)
def reads(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR_CREATED),
    target_family: str | None = typer.Option(
        None, "-tf", "--target-family", help=_HELP_TAXON.format("family")
    ),
    target_genus: str | None = typer.Option(
        None, "-tg", "--target-genus", help=_HELP_TAXON.format("genus")
    ),
    target_species: str | None = typer.Option(
        None, "-ts", "--target-species", help=_HELP_TAXON.format("species")
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
    max_bases: int | None = typer.Option(
        None,
        "--max-bases",
        min=1,
        help="Drop runs with more sequenced bases than this (a guard against whole-host "
        "libraries, which would assemble into a host-dominated genome).",
    ),
    drop_selection: list[str] = typer.Option(
        ["MDA"],
        "--drop-selection",
        help="Drop runs whose ENA library selection is this value (repeatable; default "
        "MDA, whole-genome amplification). Pass 'none' to keep every selection.",
    ),
    one_per_sample: bool = typer.Option(
        True,
        "--one-per-sample/--all-runs",
        help="Keep the best run of each sample (a long-read run with enough bases, else "
        "the largest run), or every run.",
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
            max_bases=max_bases,
            drop_selection=[s for s in drop_selection if s.lower() != "none"],
            one_per_sample=one_per_sample,
        )

    _run("reads", workdir, build, create=True)


@app.command(rich_help_panel=PANEL_ENTRY)
def assemble(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    assembler: str = typer.Option("auto", "--assembler", help=_assembler_help()),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    jobs: int | None = typer.Option(
        None,
        "--jobs",
        min=1,
        help="Runs assembled concurrently, threads split across them (default 2, or 1 when a "
        "long-read run is pending).",
    ),
    memory_gb: int = typer.Option(
        16,
        "--memory-gb",
        min=1,
        help="Memory hint per assembly, in GB, for tools that cap RAM; also the budget that "
        "caps concurrent classifier gathers (about 0.6 GB each).",
    ),
    min_contig_length: int = typer.Option(
        500, "--min-contig-length", min=0, help="Drop contigs shorter than this many bases."
    ),
    polisher: str = typer.Option("auto", "--polisher", help=_polisher_help()),
    polish_rounds: int = typer.Option(
        1, "--polish-rounds", min=1, help="Polishing rounds (racon; medaka runs one)."
    ),
    outgroup: Path | None = typer.Option(
        None, "--outgroup", help="A FASTA file to set aside as the outgroup for rooting."
    ),
    append: bool = typer.Option(
        False,
        "--append",
        help="Add the assemblies to a working directory that already holds a selection "
        "(metadata and genome, or ingest) instead of replacing it.",
    ),
    keep_reads: bool = typer.Option(
        False, "--keep-reads", help="Keep the downloaded FASTQ files after assembling."
    ),
    keep_files: bool = typer.Option(
        False, "--keep-files", help="Keep each run's assembler scratch directory."
    ),
    checkm2_db: Path | None = typer.Option(
        None,
        "--checkm2-db",
        help="CheckM2 DIAMOND database; enables quality scoring (or set CHECKM2DB).",
    ),
    min_completeness: float = typer.Option(
        50.0, "--min-completeness", min=0.0, max=100.0, help="CheckM2 completeness floor."
    ),
    max_contamination: float = typer.Option(
        10.0, "--max-contamination", min=0.0, max=100.0, help="CheckM2 contamination ceiling."
    ),
    classifier: str = typer.Option("auto", "--classifier", help=_classifier_help()),
    gtdb_sketch: Path | None = typer.Option(
        None,
        "--gtdb-sketch",
        help="GTDB sourmash sketch database (.sig.zip); enables classification "
        "(or set REPGENR_GTDB_SKETCH).",
    ),
    gtdb_lineages: Path | None = typer.Option(
        None,
        "--gtdb-lineages",
        help="The lineages CSV published with the sketch (or set REPGENR_GTDB_LINEAGES).",
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
            polisher=polisher,
            polish_rounds=polish_rounds,
            outgroup=None if outgroup is None else str(outgroup),
            append=append,
            keep_reads=keep_reads,
            keep_files=keep_files,
            checkm2_db=None if checkm2_db is None else str(checkm2_db),
            min_completeness=min_completeness,
            max_contamination=max_contamination,
            classifier=classifier,
            gtdb_sketch=None if gtdb_sketch is None else str(gtdb_sketch),
            gtdb_lineages=None if gtdb_lineages is None else str(gtdb_lineages),
            extra=_parse_key_values(tool_arg, "--tool-arg"),
        )

    _run("assemble", workdir, build)
