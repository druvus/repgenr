"""derep_unpack stage: explode clusters into one directory per representative.

Reads the derep ``clusters.tsv`` contract and copies each cluster's genomes
(optionally excluding the representative) into ``derep/unpacked/<representative>/``.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import CLUSTERS_TSV, read_clusters
from ..core.errors import WorkdirError


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
        shutil.rmtree(unpack_dir)
    unpack_dir.mkdir(parents=True)

    empty = 0
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
                shutil.copy2(source, cluster_dir / genome)
    if empty:
        logger.info("%d clusters had only a representative and were skipped", empty)
    logger.info("Unpacked %d clusters into %s", len(clusters), unpack_dir)
    return unpack_dir
