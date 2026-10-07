"""derep_unpack stage: explode clusters into one directory per representative.

Reads the derep ``clusters.tsv`` contract and hard-links (or, where the file
system cannot, copies) each cluster's genomes, optionally excluding the
representative, into ``derep/unpacked/<representative>/``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import CLUSTERS_TSV, read_clusters
from ..core.errors import WorkdirError
from ..core.process import link_or_copy, staged_dir


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
    dir_names = _cluster_dir_names(clusters)
    empty = 0
    written = 0
    copy_noted = False
    missing: list[str] = []
    # Built beside the old tree and swapped in when complete, so a failure
    # keeps the previous view instead of a partial one.
    with staged_dir(unpack_dir) as staging:
        for rep, members in clusters.items():
            targets = list(members)
            if not params.no_representant:
                targets = [rep, *members]
            if not targets:
                empty += 1
                continue
            cluster_dir = staging / dir_names[rep]
            cluster_dir.mkdir()
            for genome in targets:
                source = ctx.genomes_dir / genome
                if not source.exists():
                    missing.append(genome)
                    continue
                if not link_or_copy(source, cluster_dir / genome) and not copy_noted:
                    logger.info(
                        "Hard links are not possible from %s to %s; copying files "
                        "instead, which takes longer and uses disk space",
                        ctx.genomes_dir,
                        unpack_dir,
                    )
                    copy_noted = True
                written += 1
    _warn_missing(missing, ctx.genomes_dir, logger)
    if empty:
        logger.info("%d clusters had only a representative and were skipped", empty)
    ctx.config.record_stage(
        "derep_unpack",
        params={**asdict(params), "clusters": len(clusters)},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info(
        "Unpacked %d genome files into %d cluster directories under %s",
        written,
        len(clusters) - empty,
        unpack_dir,
    )
    return unpack_dir


def _cluster_dir_names(clusters: dict[str, list[str]]) -> dict[str, str]:
    """Directory name per representative: the file stem, or the full file name
    when two representatives share a stem (``x.fasta`` and ``x.fna``)."""
    stems = Counter(Path(rep).stem for rep in clusters)
    return {rep: Path(rep).stem if stems[Path(rep).stem] == 1 else rep for rep in clusters}


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
