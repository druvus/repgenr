"""dRep dereplication adapter.

Ports the logic of the old ``derep_worker.py``: optionally decompress gzipped
genomes (prodigal inside dRep needs plain FASTA), run ``dRep dereplicate``, then
parse ``Cdb.csv`` (cluster membership) and ``genomeInformation.csv`` (QC) to
build a normalized :class:`DerepResult`.

dRep does not scale natively to very large sets; the dereplicate stage wraps it
with chunking when needed (``supports_native_scaling=False``).
"""

from __future__ import annotations

import csv
import logging
import shutil
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from .base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
)

_FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fas")


class DrepDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="drep",
        conda=("bioconda::drep",),
        required_binaries=(BinarySpec("dRep", version_args=("--version",)),),
        default_params={"S_algorithm": "fastANI"},
        recommended_max_genomes=2000,
        supports_native_scaling=False,
        threads_param="--processors",
    )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        genomes_dir = out_dir / "genomes"
        genomes_dir.mkdir(exist_ok=True)

        # Stage plain (decompressed) copies so dRep/prodigal can read them.
        staged: list[Path] = []
        for src in genomes:
            staged.append(_stage_genome(src, genomes_dir))

        drep_wd = out_dir / "drep_workdir"
        cmd: list[str | Path] = [
            "dRep", "dereplicate", drep_wd,
            "-g", *staged,
            "--processors", str(params.threads),
            "-sa", str(params.secondary_ani),
            "-pa", str(params.primary_ani),
            "--S_algorithm", params.extra.get("S_algorithm", "fastANI"),
            "--length", str(params.extra.get("length", 0)),
        ]
        if params.extra.get("virus"):
            cmd += [
                "--S_algorithm", "ANImf",
                "--cov_thresh", "0.5",
                "--N50_weight", "0",
                "--size_weight", "1",
                "--ignoreGenomeQuality",
                "--clusterAlg", "single",
            ]
        run_tool(self.capabilities, cmd, logger=logger, log_prefix="drep")

        if not drep_wd.exists():
            raise WorkdirError(
                "dRep working directory was not created; confirm dRep is installed and runs."
            )
        return _parse_drep_output(drep_wd, logger)


def _stage_genome(src: Path, dest_dir: Path) -> Path:
    """Copy ``src`` into ``dest_dir``, decompressing a .gz on the way."""
    if src.suffix == ".gz":
        import gzip

        target = dest_dir / src.with_suffix("").name
        with gzip.open(src, "rb") as fi, open(target, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return target
    target = dest_dir / src.name
    if not target.exists():
        shutil.copy2(src, target)
    return target


def _parse_drep_output(drep_wd: Path, logger: logging.Logger) -> DerepResult:
    derep_genomes = drep_wd / "dereplicated_genomes"
    representatives = sorted(
        p for p in derep_genomes.iterdir() if p.suffix in _FASTA_SUFFIXES
    )
    rep_names = {p.name for p in representatives}

    # Cluster membership from Cdb.csv (genome, secondary_cluster).
    clusters_by_id: dict[str, set[str]] = {}
    genome_cluster: dict[str, str] = {}
    cdb = drep_wd / "data_tables" / "Cdb.csv"
    with open(cdb, newline="") as fo:
        reader = csv.DictReader(fo)
        for row in reader:
            genome = row["genome"]
            clust = row.get("secondary_cluster") or row.get("cluster") or ""
            genome_cluster[genome] = clust
            clusters_by_id.setdefault(clust, set()).add(genome)

    status: dict[str, str] = {}
    clusters: dict[str, list[str]] = {}
    for rep in representatives:
        status[rep.name] = STATUS_REPRESENTATIVE
        clust = genome_cluster.get(rep.name)
        members = clusters_by_id.get(clust, set()) if clust is not None else set()
        contained = sorted(m for m in members if m != rep.name)
        clusters[rep.name] = contained
        for m in contained:
            status[m] = STATUS_CONTAINED

    # Genomes filtered out by dRep QC have centrality 0 in genomeInformation.csv.
    info = _parse_genome_information(drep_wd, logger)
    for genome, centrality in info.items():
        if centrality == 0.0 and genome not in rep_names:
            status.setdefault(genome, STATUS_FAIL_QC)

    return DerepResult(
        representatives=representatives,
        clusters=clusters,
        genome_status=status,
        genome_information=[{"genome": g, "centrality": c} for g, c in info.items()],
    )


def _parse_genome_information(drep_wd: Path, logger: logging.Logger) -> dict[str, float]:
    """Return genome -> centrality, normalizing the viral 5-column variant.

    dRep run with ``--ignoreGenomeQuality`` writes 5 columns
    (genome_path,length,N50,genome,centrality) instead of the usual 7; we read
    the centrality field in either layout.
    """
    path = drep_wd / "data_tables" / "genomeInformation.csv"
    if not path.exists():
        logger.warning("genomeInformation.csv not found; skipping QC status")
        return {}
    out: dict[str, float] = {}
    with open(path, newline="") as fo:
        reader = csv.reader(fo)
        header = next(reader, None)
        if header is None:
            return out
        ncols = len(header)
        for row in reader:
            if ncols == 7:
                genome, centrality = row[0], row[6]
            elif ncols == 5:
                genome, centrality = row[3], row[4]
            else:
                continue
            try:
                out[genome] = float(centrality)
            except ValueError:
                out[genome] = 0.0
    return out
