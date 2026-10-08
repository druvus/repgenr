"""sketch stage: write the genome sketches of an existing working directory.

The genome-writing stages sketch their genomes when sourmash can run
(:mod:`repgenr.core.sketches`). This stage adds the sketches to a working
directory written without them, or brings them up to date: it sketches the
genomes whose sketch is missing or stale and removes the sketches of genomes
no longer in the set. ``--force`` (a global option) sketches every genome
again.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime

from ..core.context import WorkdirContext
from ..core.errors import WorkdirError
from ..core.plugins import preflight
from ..core.sketches import SKETCH_PARAMS, remove_stale, sketch_genomes, sketch_targets
from ..core.sourmash import SOURMASH_TOOL


@dataclass
class SketchParams:
    threads: int = 16


def run(ctx: WorkdirContext, params: SketchParams) -> int:
    logger = ctx.logger
    targets = sketch_targets(ctx)
    if not targets:
        raise WorkdirError(
            f"No genome files recorded in the manifest of {ctx.workdir}; run a stage that "
            "writes a genome set first (genome, vgenome, ingest or assemble)."
        )
    versions = preflight(SOURMASH_TOOL)
    removed = remove_stale(ctx, logger)
    summary = sketch_genomes(ctx, targets, params.threads, logger, force=ctx.force)
    summary.removed = removed
    logger.info("%s", summary.line())
    ctx.config.record_stage(
        "sketch",
        tool="sourmash",
        params={"sketch_params": SKETCH_PARAMS, "genomes": len(targets), **summary.as_dict()},
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return summary.written + summary.replaced + summary.forced + summary.copied
