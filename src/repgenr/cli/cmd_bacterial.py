"""Bacterial lineage commands: metadata, genome, dereplicate."""

from __future__ import annotations

from pathlib import Path

import typer

from ..core.errors import UserInputError
from .base import DEFAULT_THREADS, _require_choice, _require_unit_interval, _run, app


@app.command()
def metadata(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    dataset: str = typer.Option(..., "-d", "--dataset", help="all or rep."),
    level: str = typer.Option(..., "-l", "--level", help="family, genus or species."),
    source: str = typer.Option(
        "tsv", "--source", help="tsv (download full table) or api (GTDB API, target only)."
    ),
    release: str | None = typer.Option(None, "-r", "--release", help="GTDB release (tsv source)."),
    gtdb_version: str | None = typer.Option(
        None, "--gtdb-version", help="bac120/ar53 (tsv source)."
    ),
    target_family: str | None = typer.Option(None, "-tf", "--target-family"),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus"),
    target_species: str | None = typer.Option(None, "-ts", "--target-species"),
    outgroup_accession: str | None = typer.Option(None, "--outgroup-accession"),
    metadata_path: str | None = typer.Option(None, "--metadata-path"),
    nodownload: bool = typer.Option(False, "--nodownload"),
    limit: int | None = typer.Option(None, "--limit"),
) -> None:
    """Select a taxon's genomes from GTDB (full table or the GTDB API)."""
    from ..stages.metadata import MetadataParams

    def build() -> MetadataParams:
        return MetadataParams(
            dataset=dataset, level=level, source=source,
            release=release, version=gtdb_version,
            target_family=target_family, target_genus=target_genus,
            target_species=target_species, outgroup_accession=outgroup_accession,
            metadata_path=metadata_path, nodownload=nodownload, limit=limit,
        )

    _run("metadata", workdir, build, create=True)


@app.command()
def genome(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    accession_list_only: bool = typer.Option(False, "--accession-list-only"),
    keep_files: bool = typer.Option(False, "--keep-files"),
) -> None:
    """Download and organize genomes selected by the metadata stage."""
    from ..stages.genome import GenomeParams

    def build() -> GenomeParams:
        return GenomeParams(accession_list_only=accession_list_only, keep_files=keep_files)

    _run("genome", workdir, build)


@app.command()
def dereplicate(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    tool: str = typer.Option("skder", "--tool", help="auto/skder/drep/galah/sourmash."),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "-af", "--aligned-fraction"),
    threads: int = typer.Option(DEFAULT_THREADS, "-t", "--threads"),
    process_size: int | None = typer.Option(
        None, "-s", "--process-size",
        help="Chunk size; when set and exceeded, two-stage chunking runs for any tool.",
    ),
    num_processes: int = typer.Option(
        1, "-p", "--num-processes", help="Parallel stage-1 chunk workers (threads split across)."
    ),
    pre_primary_ani: float | None = typer.Option(
        None, "--pre-primary-ani",
        help="Stage-1 (intra-chunk) primary ANI; defaults to --primary-ani.",
    ),
    pre_secondary_ani: float | None = typer.Option(
        None, "--pre-secondary-ani",
        help="Stage-1 (intra-chunk) secondary ANI; defaults to --secondary-ani.",
    ),
    reduce: str = typer.Option(
        "none", "--reduce",
        help="Taxonomy-aware reduction after ANI: none, species, or genus "
        "(one representative per taxon).",
    ),
    target_reps: int = typer.Option(
        0, "--target-reps",
        help="Target representative count: search --secondary-ani to land near it "
        "(0 = off; re-runs dereplication per search step).",
    ),
    virus: bool = typer.Option(False, "--virus", help="Pass virus-tuned parameters to the tool."),
) -> None:
    """Cluster genomes by ANI and select representatives."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.dereplicate import DereplicateParams

    def build() -> DereplicateParams:
        _require_choice(tool, {"auto", *_derep_registry.names()}, "--tool")
        _require_choice(reduce, {"none", "species", "genus"}, "--reduce")
        if target_reps < 0:
            raise UserInputError(f"--target-reps must be >= 0, got {target_reps}.")
        _require_unit_interval(primary_ani, "--primary-ani")
        _require_unit_interval(secondary_ani, "--secondary-ani")
        _require_unit_interval(aligned_fraction, "--aligned-fraction")
        _require_unit_interval(pre_primary_ani, "--pre-primary-ani")
        _require_unit_interval(pre_secondary_ani, "--pre-secondary-ani")
        return DereplicateParams(
            tool=tool,
            primary_ani=primary_ani,
            secondary_ani=secondary_ani,
            aligned_fraction=aligned_fraction,
            threads=threads,
            process_size=process_size,
            num_processes=num_processes,
            pre_primary_ani=pre_primary_ani,
            pre_secondary_ani=pre_secondary_ani,
            reduce=reduce,
            target_reps=target_reps,
            extra={"virus": virus} if virus else {},
        )

    _run("dereplicate", workdir, build)
