"""glance stage: quick ANI overview with dRep compare.

Ports ``glance.py``: run ``dRep compare`` (mash primary clustering only) on all
genomes, copy out the clustering dendrogram, and plot a boxplot + histogram of
the all-vs-all MASH ANI similarities from ``Mdb.csv``.
"""

from __future__ import annotations

import shutil
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import list_fasta
from ..core.errors import UserInputError, WorkdirError
from ..core.process import remove_tree

GLANCE_OUTPUTS = (
    "glance_clustering_dendrogram.pdf",
    "glance_MASH_ANI_similarity_boxplot.png",
    "glance_MASH_ANI_similarity_histogram.png",
)


@dataclass
class GlanceParams:
    threads: int = 16  # same default as the CLI's -t
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
    genomes = list_fasta(ctx.genomes_dir)
    if not genomes:
        raise WorkdirError(f"No genomes under {ctx.genomes_dir}")
    versions = adapter.preflight()

    glance_wd = ctx.workdir / "glance_wd"
    if glance_wd.exists():
        remove_tree(glance_wd)

    result = adapter.compare(genomes, glance_wd, params.threads, logger)

    # The comparison succeeded: drop the previous outputs before writing new
    # ones, so a plot that is not drawn this time (no similarity in range)
    # does not survive from an earlier run.
    for name in GLANCE_OUTPUTS:
        (ctx.workdir / name).unlink(missing_ok=True)
    out_pdf = ctx.workdir / "glance_clustering_dendrogram.pdf"
    if result.dendrogram is not None:
        shutil.copy2(result.dendrogram, out_pdf)

    if result.similarity_csv is not None:
        _plot(result.similarity_csv, ctx.workdir, params, logger)
    else:
        logger.warning("No similarity table produced; skipping plots")

    if not params.keep_files and glance_wd.exists():
        remove_tree(glance_wd)
    ctx.config.record_stage(
        "glance",
        tool=params.tool,
        params=asdict(params),
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Glance outputs written to %s", ctx.workdir)
    return out_pdf


def _pair_similarities(mdb: Path, low: float, high: float) -> list[float]:
    """Similarity of each unordered genome pair in ``Mdb.csv`` within [low, high].

    dRep lists every pair in both orders and each genome against itself; a
    pair is counted once (its first row) and self-comparisons are skipped.
    """
    import csv

    index: dict[str, int] = {}
    seen: set[int] = set()
    values: list[float] = []
    with open(mdb, encoding="utf-8", newline="") as fo:
        for row in csv.DictReader(fo):
            g1, g2 = row.get("genome1"), row.get("genome2")
            if g1 is None or g2 is None or g1 == g2:
                continue
            i = index.setdefault(g1, len(index))
            j = index.setdefault(g2, len(index))
            key = (min(i, j) << 32) | max(i, j)
            if key in seen:
                continue
            seen.add(key)
            try:
                sim = float(row["similarity"])
            except (KeyError, ValueError):
                continue
            if low <= sim <= high:
                values.append(sim)
    return values


def _plot(mdb: Path, workdir: Path, params: GlanceParams, logger) -> None:
    from matplotlib import pyplot as plt

    values = _pair_similarities(mdb, params.plot_min, params.plot_max)
    if not values:
        logger.warning("No similarity values in range; skipping plots")
        return

    title = f"MASH ANI, all-vs-all ({len(values)} genome pairs)"
    fig, ax = plt.subplots()
    ax.boxplot(values)
    ax.set_xticklabels([""])
    ax.set_ylabel("MASH ANI")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(workdir / "glance_MASH_ANI_similarity_boxplot.png")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(values, bins=100)
    ax.set_xlabel("MASH ANI")
    ax.set_ylabel("Genome pairs")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(workdir / "glance_MASH_ANI_similarity_histogram.png")
    plt.close(fig)
