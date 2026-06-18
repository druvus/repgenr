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
import shutil
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.plugins import ToolCapabilities
from ..core.process import write_fofn
from .base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
)

_FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fas")


class GalahDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="galah",
        conda=("bioconda::galah",),
        required_binaries=(BinarySpec("galah", version_args=("--version",)),),
        recommended_max_genomes=None,
        supports_native_scaling=True,
        threads_param="--threads",
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

        # Pass the genome list via a file (--genome-fasta-list), never on argv:
        # 1000s-10000s of paths would exceed ARG_MAX.
        fofn = write_fofn(genomes, out_dir / "genomes.fofn")
        genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
        cmd: list[str | Path] = [
            "galah", "cluster",
            "--genome-fasta-list", fofn,
            "--ani", f"{ani_pct:g}",
            "--threads", str(params.threads),
            "--output-cluster-definition", clusters_file,
        ]
        run_tool(
            self.capabilities, cmd, logger=logger, log_prefix="galah", extra_mounts=genome_dirs
        )

        clusters: dict[str, list[str]] = {}
        status: dict[str, str] = {}
        rep_paths: dict[str, Path] = {}
        with open(clusters_file) as fo:
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

        # Stage representative FASTAs into out_dir for the contract layer.
        rep_dir = out_dir / "representatives"
        rep_dir.mkdir(exist_ok=True)
        representatives: list[Path] = []
        for name, src in rep_paths.items():
            dest = rep_dir / name
            if src.exists() and not dest.exists():
                shutil.copy2(src, dest)
            representatives.append(dest)

        return DerepResult(
            representatives=sorted(representatives),
            clusters=clusters,
            genome_status=status,
        )
