"""sourmash tree builder (alignment-free; k-mer distance + neighbor-joining)."""

from __future__ import annotations

import csv
import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import write_fofn
from ..tree.newick import neighbor_joining
from .base import InputKind, TreeBuilder, TreeParams, as_genome_list


class SourmashBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="sourmash",
        conda=("bioconda::sourmash",),
        required_binaries=(BinarySpec("sourmash", version_args=("--version",)),),
        default_params={"ksize": 31, "scaled": 1000},
        recommended_max_genomes=10000,
        threads_param=None,
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
        ksize = int(params.extra.get("ksize", self.capabilities.default_params["ksize"]))
        scaled = int(params.extra.get("scaled", self.capabilities.default_params["scaled"]))

        sig_dir = out_dir / "signatures"
        sig_dir.mkdir(exist_ok=True)
        fofn = write_fofn(genomes, out_dir / "genomes.fofn")
        run_tool(self.capabilities, 
            [
                "sourmash", "sketch", "dna",
                "-p", f"k={ksize},scaled={scaled}",
                "--from-file", fofn,
                "--outdir", sig_dir,
            ],
            logger=logger,
            log_prefix="sourmash",
        )
        # Skip macOS AppleDouble companions ("._*") that appear on exFAT/NTFS volumes.
        sigs = [
            p for p in (sorted(sig_dir.glob("*.sig")) + sorted(sig_dir.glob("*.sig.gz")))
            if not p.name.startswith("._")
        ]
        if not sigs:
            raise WorkdirError("sourmash produced no signatures")

        matrix_csv = out_dir / "compare.csv"
        run_tool(self.capabilities, 
            ["sourmash", "compare", "-k", str(ksize), "--csv", matrix_csv, *sigs],
            logger=logger,
            log_prefix="sourmash",
        )

        labels, similarity = _read_csv(matrix_csv)
        # distance = 1 - similarity
        dist = [[1.0 - similarity[i][j] for j in range(len(labels))] for i in range(len(labels))]
        clean_labels = [_label_to_genome(label, genomes) for label in labels]
        newick = neighbor_joining(clean_labels, dist)

        tree = out_dir / "tree.nwk"
        tree.write_text(newick + "\n")
        return tree


def _read_csv(path: Path) -> tuple[list[str], list[list[float]]]:
    with open(path, newline="") as fo:
        reader = csv.reader(fo)
        labels = next(reader)
        matrix = [[float(x) for x in row] for row in reader]
    return labels, matrix


def _label_to_genome(label: str, genomes: Sequence[Path]) -> str:
    base = Path(label).name
    stems = {g.stem: g.stem for g in genomes}
    if Path(label).stem in stems:
        return Path(label).stem
    return base
