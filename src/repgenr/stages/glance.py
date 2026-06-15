"""glance stage: quick ANI overview with dRep compare.

Ports ``glance.py``: run ``dRep compare`` (mash primary clustering only) on all
genomes, copy out the clustering dendrogram, and plot a boxplot + histogram of
the all-vs-all MASH ANI similarities from ``Mdb.csv``.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from ..core.binaries import BinarySpec, check_binaries
from ..core.context import WorkdirContext
from ..core.errors import WorkdirError
from ..core.process import run as run_cmd

_DREP = BinarySpec("dRep", version_args=("--version",))
_FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fas")


@dataclass
class GlanceParams:
    threads: int = 24
    plot_max: float = 1.0
    plot_min: float = 0.0
    keep_files: bool = False


def run(ctx: WorkdirContext, params: GlanceParams) -> Path:
    logger = ctx.logger
    check_binaries((_DREP,))
    genomes = [p for p in ctx.genomes_dir.iterdir() if p.suffix in _FASTA_SUFFIXES]
    if not genomes:
        raise WorkdirError(f"No genomes under {ctx.genomes_dir}")

    glance_wd = ctx.workdir / "glance_wd"
    if glance_wd.exists():
        shutil.rmtree(glance_wd)

    run_cmd(
        ["dRep", "compare", "--SkipSecondary", "-g", *genomes,
         "--processors", str(params.threads), glance_wd],
        logger=logger, log_prefix="drep",
    )

    dendrogram = glance_wd / "figures" / "Primary_clustering_dendrogram.pdf"
    out_pdf = ctx.workdir / "glance_clustering_dendrogram.pdf"
    if dendrogram.exists():
        shutil.copy2(dendrogram, out_pdf)

    mdb = glance_wd / "data_tables" / "Mdb.csv"
    if mdb.exists():
        _plot(mdb, ctx.workdir, params, logger)
    else:
        logger.warning("Mdb.csv not found; skipping plots")

    if not params.keep_files and glance_wd.exists():
        shutil.rmtree(glance_wd)
    logger.info("Glance outputs written to %s", ctx.workdir)
    return out_pdf


def _plot(mdb: Path, workdir: Path, params: GlanceParams, logger) -> None:
    import csv

    from matplotlib import pyplot as plt

    values = []
    with open(mdb, newline="") as fo:
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
