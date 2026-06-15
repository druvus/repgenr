"""mashtree tree builder (alignment-free; consumes genomes directly)."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.plugins import ToolCapabilities
from ..core.process import run
from .base import InputKind, TreeBuilder, TreeParams, as_genome_list


class MashtreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="mashtree",
        required_binaries=(BinarySpec("mashtree", version_args=("--version",)),),
        recommended_max_genomes=10000,
        threads_param="--numcpus",
    )
    input_kind = InputKind.GENOMES

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
        run(
            ["mashtree", "--numcpus", str(params.threads), *genomes],
            logger=logger,
            log_prefix="mashtree",
            stdout_path=tree,
        )
        return tree
