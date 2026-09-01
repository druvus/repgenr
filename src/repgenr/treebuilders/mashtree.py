"""mashtree tree builder (alignment-free; consumes genomes directly)."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.plugins import ToolCapabilities
from ..core.process import warn_argv_bytes
from .base import InputKind, TreeBuilder, TreeParams, as_genome_list


class MashtreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="mashtree",
        conda=("bioconda::mashtree",),
        required_binaries=(
            BinarySpec("mashtree", version_args=("--version",), min_version="1.2"),
        ),
        recommended_max_genomes=10000,
        accepted_extras=frozenset({"genomesize"}),
    )
    input_kind = InputKind.GENOMES

    def distance_matrix(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        matrix = out_dir / "distance_matrix.tsv"
        cmd: list[str | Path] = ["mashtree"]
        genomesize = params.extra.get("genomesize")
        if genomesize is not None:
            cmd += ["--genomesize", str(int(genomesize))]
        cmd += ["--mindepth", "0", "--outmatrix", matrix, *genomes]
        warn_argv_bytes("mashtree", cmd, logger)
        run_tool(self.capabilities, cmd,
            logger=logger, log_prefix="mashtree",
            stdout_path=out_dir / "mashtree.dnd",
        )
        return matrix

    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        genomes = as_genome_list(msa_or_genomes)
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        cmd: list[str | Path] = ["mashtree", "--numcpus", str(params.threads), *genomes]
        warn_argv_bytes("mashtree", cmd, logger)
        run_tool(self.capabilities,
            cmd,
            logger=logger,
            log_prefix="mashtree",
            stdout_path=tree,
        )
        return tree
