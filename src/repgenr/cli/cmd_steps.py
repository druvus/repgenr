"""Stateless data-channel steps: genome-fetch, dereplicate-chunk, dereplicate-merge.

These run as discrete Nextflow process steps (no shared workdir); they read
explicit inputs (a selection.tsv or a file-of-filenames) and write a result
directory, rather than going through the workdir-bound :func:`_run` harness.
"""

from __future__ import annotations

from pathlib import Path

import typer

from ..core.errors import RepGenRError, UserInputError
from ..core.logging import configure_logging
from .base import _RUN_STATE, _read_path_fofn, _require_choice, _require_unit_interval, app


@app.command(name="genome-fetch")
def genome_fetch_cmd(
    selection: Path = typer.Option(
        ..., "--selection", help="selection.tsv from the metadata stage."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir for downloaded genomes."),
    keep_files: bool = typer.Option(False, "--keep-files", help="Keep download intermediates."),
) -> None:
    """Download genomes listed in a selection.tsv (stateless data-channel step)."""
    from ..stages.genome_steps import GenomeFetchParams, genome_fetch

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    try:
        genome_fetch(
            GenomeFetchParams(selection_tsv=selection, out_dir=out_dir, keep_files=keep_files),
            logger,
        )
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


@app.command(name="dereplicate-chunk")
def dereplicate_chunk_cmd(
    genomes_fofn: Path = typer.Option(
        ..., "--genomes-fofn", help="File of genome FASTA paths, one per line."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output directory for the chunk result."),
    tool: str = typer.Option("skder", "--tool", help="skder/drep/galah/sourmash."),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "-af", "--aligned-fraction"),
    threads: int = typer.Option(16, "-t", "--threads"),
    virus: bool = typer.Option(False, "--virus", help="Pass virus-tuned parameters to the tool."),
) -> None:
    """Dereplicate one chunk of genomes (scatter step; writes a chunk result dir)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.derep_steps import ChunkParams, dereplicate_chunk

    # A concrete tool (not 'auto') so every scattered chunk and the merge agree.
    _require_choice(tool, set(_derep_registry.names()), "--tool")
    _require_unit_interval(primary_ani, "--primary-ani")
    _require_unit_interval(secondary_ani, "--secondary-ani")
    _require_unit_interval(aligned_fraction, "--aligned-fraction")
    genomes = _read_path_fofn(genomes_fofn)

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    try:
        dereplicate_chunk(
            ChunkParams(
                tool=tool, genomes=genomes, out_dir=out_dir,
                primary_ani=primary_ani, secondary_ani=secondary_ani,
                aligned_fraction=aligned_fraction, threads=threads,
                extra={"virus": virus} if virus else {},
            ),
            logger,
        )
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


@app.command(name="dereplicate-merge")
def dereplicate_merge_cmd(
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir for the merged result."),
    chunk_dir: list[Path] = typer.Option(
        [], "--chunk-dir", help="A chunk result directory (repeatable)."
    ),
    chunk_fofn: Path | None = typer.Option(
        None, "--chunk-fofn", help="File listing chunk result directories, one per line."
    ),
    tool: str = typer.Option("skder", "--tool", help="skder/drep/galah/sourmash."),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "-af", "--aligned-fraction"),
    threads: int = typer.Option(16, "-t", "--threads"),
    virus: bool = typer.Option(False, "--virus", help="Pass virus-tuned parameters to the tool."),
) -> None:
    """Dereplicate the union of chunk representatives (gather step)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.derep_steps import MergeParams, dereplicate_merge

    _require_choice(tool, set(_derep_registry.names()), "--tool")
    _require_unit_interval(primary_ani, "--primary-ani")
    _require_unit_interval(secondary_ani, "--secondary-ani")
    _require_unit_interval(aligned_fraction, "--aligned-fraction")
    chunk_dirs = list(chunk_dir)
    if chunk_fofn is not None:
        chunk_dirs += _read_path_fofn(chunk_fofn)
    if not chunk_dirs:
        raise UserInputError("Provide at least one --chunk-dir or a --chunk-fofn.")

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    try:
        dereplicate_merge(
            MergeParams(
                tool=tool, chunk_dirs=chunk_dirs, out_dir=out_dir,
                primary_ani=primary_ani, secondary_ani=secondary_ani,
                aligned_fraction=aligned_fraction, threads=threads,
                extra={"virus": virus} if virus else {},
            ),
            logger,
        )
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc
