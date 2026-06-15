"""FastTree tree builder (approximate ML from an MSA; fast at scale)."""

from __future__ import annotations

import logging
import shutil
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.plugins import ToolCapabilities
from ..core.process import run
from .base import InputKind, TreeBuilder, TreeParams, as_msa_path


class FasttreeBuilder(TreeBuilder):
    # FastTree / VeryFastTree are both accepted; prefer FastTree on PATH.
    capabilities = ToolCapabilities(
        name="fasttree",
        required_binaries=(BinarySpec("FastTree", version_args=()),),
        recommended_max_genomes=5000,
        threads_param=None,
    )
    input_kind = InputKind.MSA_FASTA

    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        msa = as_msa_path(msa_or_genomes)
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        binary = "FastTree" if shutil.which("FastTree") else "fasttree"
        run(
            [binary, "-nt", "-gtr", msa],
            logger=logger,
            log_prefix="fasttree",
            stdout_path=tree,
        )
        return tree
