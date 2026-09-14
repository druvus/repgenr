"""``repgenr reads``: select sequencing runs from ENA/SRA for the reads chain."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import HELP_TARGET_FAMILY, HELP_TARGET_GENUS, HELP_TARGET_SPECIES, _run, app


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
