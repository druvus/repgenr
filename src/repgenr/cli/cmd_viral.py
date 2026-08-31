"""Viral lineage commands: vmetadata, vgenome."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import typer

from .base import _run, app


def _validate_released_after(value: str | None) -> str | None:
    """Reject a malformed --released-after date before any network work runs."""
    if value is None:
        return value
    try:
        datetime.strptime(value, "%m/%d/%Y")
    except ValueError as exc:
        raise typer.BadParameter(
            f"expected MM/DD/YYYY (e.g. 01/31/2024), got '{value}'."
        ) from exc
    return value


@app.command()
def vmetadata(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    target: str | None = typer.Option(None, "-t", "--target", help="Virus taxon/group/family."),
    source: str = typer.Option(
        "ncbi_virus", "--source", help="ncbi_virus (NCBI Virus via datasets) or bvbrc."
    ),
    filter: str = typer.Option("complete genome", "-f", "--filter", help="BV-BRC header tag."),
    host: str | None = typer.Option(None, "--host", help="ncbi_virus: restrict to a host species."),
    complete_only: bool = typer.Option(
        False, "--complete-only", help="ncbi_virus: only COMPLETE sequences."
    ),
    released_after: str | None = typer.Option(
        None, "--released-after", callback=_validate_released_after,
        help="ncbi_virus: MM/DD/YYYY.",
    ),
    list_targets: bool = typer.Option(False, "-l", "--list", help="List BV-BRC targets and exit."),
) -> None:
    """Retrieve viral metadata from NCBI Virus (default) or BV-BRC."""
    from .param_builders import vmetadata_params

    def build():
        return vmetadata_params(
            target=target, filter=filter, list_targets=list_targets,
            source=source, host=host, complete_only=complete_only,
            released_after=released_after,
        )

    _run("vmetadata", workdir, build, create=True)


@app.command()
def vgenome(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus"),
    target_species: str | None = typer.Option(None, "-ts", "--target-species"),
    target_serotype: str | None = typer.Option(None, "-tse", "--target-serotype"),
    target_custom: str | None = typer.Option(None, "-tc", "--target-custom", help="key:value."),
    length_all: bool = typer.Option(False, "--length-all"),
    length_deviation: int = typer.Option(10, "--length-deviation", min=0),
    length_method: str = typer.Option(
        "median_of_medians", "--length-method",
        help="Center of the length window: median_of_medians (default) gives "
        "one vote per species, so an over-sequenced outbreak species cannot "
        "shift the window; mean averages every record and loses that defense.",
    ),
    length_range: str | None = typer.Option(None, "--length-range", help="e.g. 25000-35000."),
    discard: str | None = typer.Option(None, "--discard", help="Comma-separated header tags."),
    no_outgroup: bool = typer.Option(False, "--no-outgroup"),
    group_segments: bool = typer.Option(
        False, "--group-segments",
        help="ncbi_virus: combine an isolate's segments into one genome (segmented viruses).",
    ),
    min_outgroup_genomes: int = typer.Option(5, "--outgroup-candidates-taxid-min-genomes"),
    outgroup_treebuilder: str = typer.Option(
        "mashtree", "--outgroup-treebuilder",
        help="Tree builder used for the outgroup distance matrix.",
    ),
    glance: bool = typer.Option(False, "--glance", help="Print selection and stop."),
    print_fasta_headers: bool = typer.Option(False, "--print-fasta-headers"),
    ignore_duplicates: bool = typer.Option(False, "--ignore-duplicates"),
    keep_files: bool = typer.Option(False, "--keep-files"),
) -> None:
    """Select and organize viral genomes (virus equivalent of genome)."""
    from .param_builders import vgenome_params

    def build():
        return vgenome_params(
            target_genus=target_genus, target_species=target_species,
            target_serotype=target_serotype, target_custom=target_custom,
            length_all=length_all, length_deviation=length_deviation,
            length_method=length_method, length_range=length_range, discard=discard,
            no_outgroup=no_outgroup, group_segments=group_segments,
            outgroup_candidates_taxid_min_genomes=min_outgroup_genomes,
            outgroup_treebuilder=outgroup_treebuilder,
            glance=glance, print_fasta_headers=print_fasta_headers,
            ignore_duplicates=ignore_duplicates, keep_files=keep_files,
        )

    _run("vgenome", workdir, build)
