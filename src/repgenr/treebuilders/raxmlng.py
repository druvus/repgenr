"""RAxML-NG tree builder (maximum likelihood from an MSA)."""

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


class RaxmlNgBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="raxmlng",
        required_binaries=(BinarySpec("raxml-ng", version_args=("--version",)),),
        recommended_max_genomes=1000,
        threads_param="--threads",
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
        prefix = out_dir / "raxml"
        cmd: list[str | Path] = [
            "raxml-ng", "--all",
            "--msa", msa,
            "--model", params.extra.get("model", "GTR+G"),
            # auto{N}: let RAxML-NG pick an efficient thread count up to the
            # budget, avoiding its core-oversubscription guard on small alignments.
            "--threads", f"auto{{{params.threads}}}",
            "--prefix", prefix,
            "--redo",  # overwrite any outputs from a previous run at this prefix
        ]
        if params.bootstrap > 0:
            cmd += ["--bs-trees", str(params.bootstrap)]
        if params.outgroup:
            cmd += ["--outgroup", params.outgroup]
        run(cmd, logger=logger, log_prefix="raxml-ng")

        best = Path(str(prefix) + ".raxml.bestTree")
        if not best.exists():
            raise WorkdirError("RAxML-NG did not produce a bestTree")
        tree = out_dir / "tree.nwk"
        shutil.copy2(best, tree)
        return tree
