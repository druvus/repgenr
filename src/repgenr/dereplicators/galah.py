"""galah dereplication adapter.

galah is a fast dRep-style clusterer. CLI used::

    galah cluster --genome-fasta-files <files> --ani <ANI%> \
          --output-cluster-definition <clusters.tsv> \
          --output-representative-list <reps.txt> --threads <n>

The cluster-definition file is ``representative<TAB>member`` rows.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping, Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import strip_fasta_suffix
from ..core.plugins import ToolCapabilities
from ..core.process import write_fofn
from .base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
    write_genome_info,
)


class GalahDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="galah",
        container="quay.io/biocontainers/galah:0.4.2--hc1c3326_2",
        conda=("bioconda::galah",),
        # 0.4: --precluster-ani and --min-aligned-fraction as passed above.
        required_binaries=(BinarySpec("galah", version_args=("--version",), min_version="0.4"),),
        recommended_max_genomes=None,
        supports_native_scaling=True,
    )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        clusters_file = out_dir / "galah_clusters.tsv"
        sani = params.secondary_ani
        ani_pct = sani * 100 if sani <= 1.0 else sani
        pani = params.primary_ani
        precluster_pct = pani * 100 if pani <= 1.0 else pani
        af = params.aligned_fraction
        af_pct = af * 100 if af <= 1.0 else af

        # Pass the genome list via a file (--genome-fasta-list), never on argv:
        # 1000s-10000s of paths would exceed ARG_MAX.
        genome_info = _write_galah_genome_info(genomes, params.quality, out_dir, logger)
        order = list(genomes) if genome_info else _by_size(genomes)
        fofn = write_fofn(order, out_dir / "genomes.fofn")
        genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
        cmd: list[str | Path] = [
            "galah",
            "cluster",
            "--genome-fasta-list",
            fofn,
            "--ani",
            f"{ani_pct:g}",
            # galah pre-clusters at a looser ANI before the exact pass and
            # drops pairs below a minimum aligned fraction; both are percentages.
            "--precluster-ani",
            f"{precluster_pct:g}",
            "--min-aligned-fraction",
            f"{af_pct:g}",
            "--threads",
            str(params.threads),
            "--output-cluster-definition",
            clusters_file,
        ]
        if genome_info is not None:
            cmd += ["--genome-info", genome_info]
            genome_dirs.append(str(genome_info))
        run_tool(
            self.capabilities, cmd, logger=logger, log_prefix="galah", extra_mounts=genome_dirs
        )

        clusters: dict[str, list[str]] = {}
        status: dict[str, str] = {}
        rep_paths: dict[str, Path] = {}
        with open(clusters_file, encoding="utf-8") as fo:
            for line in fo:
                line = line.strip()
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 2:
                    continue
                rep_path, member_path = Path(fields[0]), Path(fields[1])
                rep_name, member_name = rep_path.name, member_path.name
                rep_paths[rep_name] = rep_path
                clusters.setdefault(rep_name, [])
                status.setdefault(rep_name, STATUS_REPRESENTATIVE)
                if member_name != rep_name:
                    clusters[rep_name].append(member_name)
                    status[member_name] = STATUS_CONTAINED

        # galah names each representative by the path it was given, so the
        # input genome is the representative file; the contract layer links it.
        # (Copying every representative into out_dir doubled their disk use.)
        by_name = {Path(g).name: Path(g) for g in genomes}
        representatives = [by_name.get(name, src) for name, src in rep_paths.items()]

        return DerepResult(
            representatives=sorted(representatives),
            clusters=clusters,
            genome_status=status,
        )


def _write_galah_genome_info(
    genomes: Sequence[Path],
    quality: Mapping[str, tuple[float, float]],
    out_dir: Path,
    logger: logging.Logger,
) -> Path | None:
    """Write galah's ``--genome-info`` table when every genome has quality.

    galah then ranks genomes by its quality formula (completeness minus 5 x
    contamination, less small penalties for contig count and ambiguous bases)
    for both the representative and the greedy membership, instead of by
    input order. galah stops when a genome has no row, so a partial table is
    never written. Rows name the genome without its FASTA suffix, ``.gz``
    included (``x.fna.gz`` -> ``x``), which is how galah matches them and the
    same rule as :func:`~repgenr.core.contracts.strip_fasta_suffix` (verified
    for every suffix in ``FASTA_SUFFIXES`` with galah 0.4.2 and 0.5.2; a row
    named ``x.fna`` for ``x.fna.gz`` is not found and galah stops).
    """
    if not genomes or not all(Path(g).name in quality for g in genomes):
        return None
    path = write_genome_info(
        out_dir / "genome_info.csv",
        ((strip_fasta_suffix(Path(g).name), *quality[Path(g).name]) for g in genomes),
    )
    logger.info(
        "galah uses the manifest completeness and contamination of %d genomes (--genome-info)",
        len(genomes),
    )
    return path


def _by_size(genomes: Sequence[Path]) -> list[Path]:
    """Genomes by descending file size, then name.

    Without genome quality galah prefers genomes listed earlier as cluster
    representatives, so a name-sorted list can make a fragment that sorts
    first the representative. The size of a gzipped file is its compressed
    size, a rough proxy only.
    """
    return sorted(genomes, key=lambda g: (-_file_size(g), Path(g).name))


def _file_size(path: Path) -> int:
    try:
        return os.path.getsize(path)
    except OSError:
        return 0
