"""Shared mashtree scaffolding for viral outgroup selection.

Both viral back-ends (NCBI Virus records in :mod:`.selection`, legacy BV-BRC
tables in :mod:`.bvbrc`) pick an outgroup the same way: stage length-compatible
candidate sequences prefixed ``S_`` (selected) / ``O_`` (outgroup candidate)
into a scratch directory, run mashtree for a distance matrix, and take the most
distant ``O_`` entry. The scaffolding lives here; the back-ends keep their own
candidate models and failure policies.

Two tuning constants deliberately differ between the back-ends and are NOT
unified here, because changing either silently changes which genomes existing
users get:

* ``RECORDS_LENGTH_TOLERANCE`` (0.15): the records back-end accepts candidates
  within 15 percent of the selection length-range midpoint.
* ``BVBRC_LENGTH_TOLERANCE`` (0.10): the BV-BRC back-end accepts candidates
  within 10 percent of the candidate taxid's median sequence length.

The back-ends also differ on failure policy (records: warn and proceed without
an outgroup; BV-BRC: raise) -- that stays in each back-end as well.
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Sequence
from pathlib import Path

from ..core.containers import run_tool
from ..core.plugins import preflight
from ..treebuilders.mashtree import MashtreeBuilder

RECORDS_LENGTH_TOLERANCE = 0.15
BVBRC_LENGTH_TOLERANCE = 0.10


def preflight_mashtree() -> dict[str, str]:
    """Resolve mashtree availability and version.

    Container-aware: checks the image when a backend is active, the host binary
    otherwise, and mashtree later runs wherever this check passed.
    """
    return preflight(MashtreeBuilder.capabilities)


def prepare_workdir(ctx) -> tuple[Path, Path]:
    """Reset the outgroup scratch dir; return (outgroup_wd, genomes_dir)."""
    outgroup_wd = ctx.workdir / "virus_outgroup_wd"
    if outgroup_wd.exists():
        shutil.rmtree(outgroup_wd)
    genomes_dir = outgroup_wd / "genomes"
    genomes_dir.mkdir(parents=True)
    return outgroup_wd, genomes_dir


def within_length_tolerance(length: float, center: float, tolerance: float) -> bool:
    """True when ``length`` is within ``center`` plus/minus ``tolerance`` (fractional)."""
    return center * (1 - tolerance) <= length <= center * (1 + tolerance)


def run_mashtree_matrix(
    genome_files: Sequence[Path],
    genomesize: int,
    outgroup_wd: Path,
    logger: logging.Logger,
) -> Path:
    """Run mashtree over the staged candidates; return the distance-matrix path.

    The caller decides how to treat a missing matrix file (warn vs raise).
    """
    matrix = outgroup_wd / "distance_matrix.tsv"
    run_tool(MashtreeBuilder.capabilities,
        ["mashtree", "--genomesize", str(genomesize), "--mindepth", "0",
         "--outmatrix", matrix, *genome_files],
        logger=logger, log_prefix="mashtree", stdout_path=outgroup_wd / "mashtree.dnd",
    )
    return matrix


def cleanup_workdir(outgroup_wd: Path, keep_files: bool) -> None:
    """Remove the outgroup scratch dir unless the user asked to keep it."""
    if not keep_files:
        shutil.rmtree(outgroup_wd, ignore_errors=True)
