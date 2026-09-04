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

from ..core.errors import UserInputError
from ..core.plugins import preflight

RECORDS_LENGTH_TOLERANCE = 0.15
BVBRC_LENGTH_TOLERANCE = 0.10


def resolve_outgroup_builder(tool: str):
    """Resolve the tree builder used for the outgroup distance matrix.

    Any registered tree builder implementing ``distance_matrix`` works
    (``--outgroup-treebuilder``); mashtree is the default.
    """
    from ..treebuilders.base import TreeBuilder, registry

    builder = registry.create(tool)
    if type(builder).distance_matrix is TreeBuilder.distance_matrix:
        supporters = sorted(
            name
            for name in registry.names()
            if not registry.is_broken(name)
            and registry.get(name).distance_matrix is not TreeBuilder.distance_matrix
        )
        raise UserInputError(
            f"Tree builder '{tool}' cannot produce the outgroup distance "
            f"matrix. Tools with distance-matrix support: "
            f"{', '.join(supporters) or 'none'}."
        )
    return builder


def preflight_outgroup_builder(builder) -> dict[str, str]:
    """Resolve the builder's availability and version.

    Container-aware: checks the image when a backend is active, the host binary
    otherwise, and the tool later runs wherever this check passed.
    """
    return preflight(builder.capabilities)


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


def run_distance_matrix(
    builder,
    genome_files: Sequence[Path],
    genomesize: int,
    outgroup_wd: Path,
    logger: logging.Logger,
) -> Path:
    """Run the builder's distance matrix over the staged candidates.

    The caller decides how to treat a missing matrix file (warn vs raise).
    """
    from ..treebuilders.base import TreeParams

    return builder.distance_matrix(
        genome_files,
        outgroup_wd,
        TreeParams(extra={"genomesize": genomesize}),
        logger,
    )


def cleanup_workdir(outgroup_wd: Path, keep_files: bool) -> None:
    """Remove the outgroup scratch dir unless the user asked to keep it."""
    if not keep_files:
        shutil.rmtree(outgroup_wd, ignore_errors=True)
