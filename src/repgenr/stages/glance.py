"""glance stage: quick ANI overview with dRep compare.

Ports ``glance.py``: run ``dRep compare`` (mash primary clustering only) on all
genomes, copy out the clustering dendrogram, and plot a boxplot + histogram of
the all-vs-all MASH ANI similarities from ``Mdb.csv``.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import list_fasta
from ..core.errors import UserInputError, WorkdirError


@dataclass
class GlanceParams:
    threads: int = 24
    tool: str = "drep"
    plot_max: float = 1.0
    plot_min: float = 0.0
    keep_files: bool = False


def run(ctx: WorkdirContext, params: GlanceParams) -> Path:
    logger = ctx.logger
    from ..dereplicators.base import Dereplicator, registry

    adapter = registry.create(params.tool)
    if type(adapter).compare is Dereplicator.compare:
        supporters = sorted(
            name
            for name in registry.names()
            if not registry.is_broken(name)
            and registry.get(name).compare is not Dereplicator.compare
        )
        raise UserInputError(
            f"Dereplicator '{params.tool}' does not support glance comparisons. "
            f"Tools with compare support: {', '.join(supporters) or 'none'}."
        )
    adapter.preflight()
    genomes = list_fasta(ctx.genomes_dir)
    if not genomes:
        raise WorkdirError(f"No genomes under {ctx.genomes_dir}")

    glance_wd = ctx.workdir / "glance_wd"
    if glance_wd.exists():
        shutil.rmtree(glance_wd)

    result = adapter.compare(genomes, glance_wd, params.threads, logger)

    out_pdf = ctx.workdir / "glance_clustering_dendrogram.pdf"
    if result.dendrogram is not None:
        shutil.copy2(result.dendrogram, out_pdf)

    if result.similarity_csv is not None:
        _plot(result.similarity_csv, ctx.workdir, params, logger)
    else:
        logger.warning("No similarity table produced; skipping plots")

    if not params.keep_files and glance_wd.exists():
        shutil.rmtree(glance_wd)
    logger.info("Glance outputs written to %s", ctx.workdir)
    return out_pdf


def _plot(mdb: Path, workdir: Path, params: GlanceParams, logger) -> None:
    import csv

    from matplotlib import pyplot as plt

    values = []
    with open(mdb, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo)
        for row in reader:
            if row.get("genome1") == row.get("genome2"):
                continue
            try:
                sim = float(row["similarity"])
            except (KeyError, ValueError):
                continue
            if params.plot_min <= sim <= params.plot_max:
                values.append(sim)

    if not values:
        logger.warning("No similarity values in range; skipping plots")
        return

    title = f"MASH ANI all-vs-all ({len(values)} values)"
    fig, ax = plt.subplots()
    ax.boxplot(values)
    ax.set_ylabel("MASH ANI")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(workdir / "glance_MASH_ANI_similarity_boxplot.png")

    fig, ax = plt.subplots()
    ax.hist(values, bins=100)
    ax.set_ylabel("MASH ANI")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(workdir / "glance_MASH_ANI_similarity_histogram.png")
