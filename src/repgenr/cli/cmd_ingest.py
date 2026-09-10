"""``repgenr ingest``: start a working directory from local genome FASTAs."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import _run, app


@app.command()
def ingest(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    genomes_dir: Path = typer.Option(
        ..., "--genomes-dir", help="Directory of genome FASTA files to stage under genomes/."
    ),
    selection: Path | None = typer.Option(
        None,
        "--selection",
        help="selection.tsv (accession, taxonomy, filename, outgroup flag, quality) "
        "naming the genomes to take; default: every FASTA under --genomes-dir, "
        "taxonomy parsed from canonical Family_genus_species_ACCESSION names.",
    ),
    outgroup: str | None = typer.Option(
        None,
        "--outgroup",
        help="Genome to set aside as the outgroup: a filename, stem or accession under "
        "--genomes-dir, or a path to a FASTA file elsewhere.",
    ),
    copy: bool = typer.Option(
        False, "--copy", help="Copy the files into genomes/ instead of symlinking them."
    ),
) -> None:
    """Populate a working directory from local genomes (no download)."""
    from .param_builders import ingest_params

    def build():
        return ingest_params(
            genomes_dir=str(genomes_dir),
            selection=None if selection is None else str(selection),
            outgroup=outgroup,
            copy=copy,
        )

    _run("ingest", workdir, build, create=True)
