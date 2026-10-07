"""derep_unpack stage: explode clusters into one directory per representative.

Reads the derep ``clusters.tsv`` contract and copies each cluster's genomes
(optionally excluding the representative) into ``derep/unpacked/<representative>/``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import CLUSTERS_TSV, read_clusters
from ..core.errors import WorkdirError
from ..core.process import link_or_copy, remove_tree


@dataclass
class DerepUnpackParams:
    no_representant: bool = False


def run(ctx: WorkdirContext, params: DerepUnpackParams) -> Path:
    logger = ctx.logger
    clusters_file = ctx.derep_dir / CLUSTERS_TSV
    if not clusters_file.exists():
        raise WorkdirError(f"Missing {clusters_file}. Run the dereplicate stage first.")
    clusters = read_clusters(clusters_file)

    unpack_dir = ctx.derep_dir / "unpacked"
    if unpack_dir.exists():
        remove_tree(unpack_dir)
    unpack_dir.mkdir(parents=True)

    empty = 0
    missing: list[str] = []
    for rep, members in clusters.items():
        targets = list(members)
        if not params.no_representant:
            targets = [rep, *members]
        if not targets:
            empty += 1
            continue
        cluster_dir = unpack_dir / Path(rep).stem
        cluster_dir.mkdir()
        for genome in targets:
            source = ctx.genomes_dir / genome
            if source.exists():
                link_or_copy(source, cluster_dir / genome)
            else:
                missing.append(genome)
    _warn_missing(missing, ctx.genomes_dir, logger)
    if empty:
        logger.info("%d clusters had only a representative and were skipped", empty)
    ctx.config.record_stage(
        "derep_unpack",
        params={**asdict(params), "clusters": len(clusters)},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Unpacked %d clusters into %s", len(clusters), unpack_dir)
    return unpack_dir


_MAX_MISSING_LINES = 10


def _warn_missing(missing: list[str], genomes_dir: Path, logger) -> None:
    """Name the cluster members that are not under ``genomes/``.

    One line per genome, or one line listing them all when there are more
    than ten, so a large gap does not flood the console.
    """
    if not missing:
        return
    if len(missing) > _MAX_MISSING_LINES:
        logger.warning(
            "%d cluster members are not in %s and were left out: %s",
            len(missing),
            genomes_dir,
            ", ".join(missing),
        )
        return
    for genome in missing:
        logger.warning("Cluster member %s is not in %s; left out", genome, genomes_dir)
