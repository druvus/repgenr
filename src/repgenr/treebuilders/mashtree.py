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
        required_binaries=(BinarySpec("mashtree", version_args=("--version",), min_version="1.2"),),
        recommended_max_genomes=10000,
        accepted_extras=frozenset({"genomesize"}),
    )
    input_kind = InputKind.GENOMES

    @staticmethod
    def _command(params: TreeParams, matrix: Path, genomes: Sequence[Path]) -> list[str | Path]:
        """One argument set for both entry points, so the tree and the matrix are
        computed under the same settings. mindepth 0 uses every k-mer of an
        assembly (there is no read depth to filter on); genomesize is an
        optional extra for mash's distance estimate.
        """
        cmd: list[str | Path] = ["mashtree", "--numcpus", str(params.threads), "--mindepth", "0"]
        genomesize = params.extra.get("genomesize")
        if genomesize is not None:
            cmd += ["--genomesize", str(int(genomesize))]
        cmd += ["--outmatrix", matrix, *genomes]
        return cmd

    def _run(
        self, genomes: Sequence[Path], out_dir: Path, params: TreeParams, logger: logging.Logger
    ) -> tuple[Path, Path]:
        """Run mashtree once; the tree goes to stdout and the matrix to a file."""
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        matrix = out_dir / "distance_matrix.tsv"
        cmd = self._command(params, matrix, genomes)
        warn_argv_bytes("mashtree", cmd, logger)
        run_tool(
            self.capabilities,
            cmd,
            logger=logger,
            log_prefix="mashtree",
            stdout_path=tree,
        )
        return tree, matrix

    def distance_matrix(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        _tree, matrix = self._run(genomes, out_dir, params, logger)
        return matrix

    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        tree, _matrix = self._run(as_genome_list(msa_or_genomes), out_dir, params, logger)
        return tree
