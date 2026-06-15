"""IQ-TREE tree builder (maximum likelihood from an MSA)."""

from __future__ import annotations

import logging
import shutil
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import run
from .base import InputKind, TreeBuilder, TreeParams, as_msa_path


class IqtreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="iqtree",
        required_binaries=(BinarySpec("iqtree", version_args=("--version",)),),
        recommended_max_genomes=500,
        threads_param="-T",
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
        work_msa = out_dir / "msa.fasta"
        if work_msa.resolve() != msa.resolve():
            shutil.copy2(msa, work_msa)

        cmd: list[str | Path] = [
            "iqtree", "-T", "auto", "--threads-max", str(params.threads), "-s", work_msa,
        ]
        if params.outgroup:
            cmd += ["-o", params.outgroup]
        if params.bootstrap > 0:
            cmd += ["-B", str(params.bootstrap)]
        run(cmd, logger=logger, log_prefix="iqtree")

        treefile = work_msa.with_suffix(".fasta.treefile")
        if not treefile.exists():
            treefile = Path(str(work_msa) + ".treefile")
        if not treefile.exists():
            raise WorkdirError("IQ-TREE did not produce a .treefile")
        tree = out_dir / "tree.nwk"
        shutil.copy2(treefile, tree)
        return tree
