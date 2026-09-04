"""sourmash tree builder (alignment-free; k-mer distance + neighbor-joining)."""

from __future__ import annotations

import csv
import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import atomic_replace
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_int
from ..core.process import write_fofn
from ..tree.newick import neighbor_joining
from .base import InputKind, TreeBuilder, TreeParams, as_genome_list

# Measured (docs/scaling-audit.md): the pure-Python NJ is ~30 s of the 142 s
# total at n=1000 and cubic beyond; refuse sizes that extrapolate to hours.
_NJ_MAX_GENOMES = 5000


class SourmashBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="sourmash",
        conda=("bioconda::sourmash",),
        accepted_extras=frozenset({"ksize", "scaled"}),
        required_binaries=(BinarySpec("sourmash", version_args=("--version",), min_version="4.0"),),
        default_params={"ksize": 31, "scaled": 1000},
        recommended_max_genomes=2000,
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
        if len(genomes) > _NJ_MAX_GENOMES:
            raise WorkdirError(
                f"The sourmash tree builder runs a pure-Python O(n^3) neighbor "
                f"joining; {len(genomes)} genomes extrapolate to many hours "
                f"(measured fit: ~1.6 h at 5000, ~12 h at 10000). Use a tool "
                f"built for this size (e.g. --treebuilder mashtree), or "
                f"dereplicate further first."
            )
        out_dir.mkdir(parents=True, exist_ok=True)
        ksize = parse_extra_int(params.extra, "ksize", self.capabilities.default_params["ksize"])
        scaled = parse_extra_int(params.extra, "scaled", self.capabilities.default_params["scaled"])

        sig_dir = out_dir / "signatures"
        sig_dir.mkdir(exist_ok=True)
        fofn = write_fofn(genomes, out_dir / "genomes.fofn")
        # Genome paths live inside the fofn (not argv); declare their dirs so the
        # container backend binds them (un-resolved abspaths, matching write_fofn).
        genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
        run_tool(
            self.capabilities,
            [
                "sourmash",
                "sketch",
                "dna",
                "-p",
                f"k={ksize},scaled={scaled}",
                "--from-file",
                fofn,
                "--outdir",
                sig_dir,
            ],
            logger=logger,
            log_prefix="sourmash",
            extra_mounts=genome_dirs,
        )
        # Skip macOS AppleDouble companions ("._*") that appear on exFAT/NTFS volumes.
        sigs = [
            p
            for p in (sorted(sig_dir.glob("*.sig")) + sorted(sig_dir.glob("*.sig.gz")))
            if not p.name.startswith("._")
        ]
        if not sigs:
            raise WorkdirError("sourmash produced no signatures")

        matrix_csv = out_dir / "compare.csv"
        compare_fofn = write_fofn(sigs, out_dir / "signatures.fofn")
        run_tool(
            self.capabilities,
            [
                "sourmash",
                "compare",
                "-k",
                str(ksize),
                "--csv",
                matrix_csv,
                "--from-file",
                compare_fofn,
            ],
            logger=logger,
            log_prefix="sourmash",
        )

        labels, similarity = _read_csv(matrix_csv)
        # distance = 1 - similarity
        dist = [[1.0 - similarity[i][j] for j in range(len(labels))] for i in range(len(labels))]
        clean_labels = [_label_to_genome(label, genomes) for label in labels]
        newick = neighbor_joining(clean_labels, dist)

        tree = out_dir / "tree.nwk"
        with atomic_replace(tree) as fo:
            fo.write(newick + "\n")
        return tree


def _read_csv(path: Path) -> tuple[list[str], list[list[float]]]:
    with open(path, encoding="utf-8", newline="") as fo:
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
