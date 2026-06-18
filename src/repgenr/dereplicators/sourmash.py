"""sourmash dereplication adapter (k-mer / MinHash ANI clustering).

Sketches each genome, runs an all-vs-all ``sourmash compare`` to get a
containment/similarity matrix, then greedily picks representatives: walk genomes
ordered by connectivity, and absorb any genome whose similarity to a chosen
representative is above the threshold. Fast and low-memory, useful as a
scalable alternative to alignment-based ANI.
"""

from __future__ import annotations

import csv
import logging
import os
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import write_fofn
from .base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
)


class SourmashDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="sourmash",
        conda=("bioconda::sourmash",),
        required_binaries=(
            BinarySpec("sourmash", version_args=("--version",), min_version="4.0"),
        ),
        default_params={"ksize": 31, "scaled": 1000},
        recommended_max_genomes=None,
        supports_native_scaling=True,
        threads_param=None,
    )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        ksize = int(params.extra.get("ksize", self.capabilities.default_params["ksize"]))
        scaled = int(params.extra.get("scaled", self.capabilities.default_params["scaled"]))

        sig_dir = out_dir / "signatures"
        sig_dir.mkdir(exist_ok=True)
        fofn = write_fofn(genomes, out_dir / "genomes.fofn")

        # The genome paths live inside the fofn, not in argv, so the container
        # backend cannot infer their mounts; declare their directories (using
        # un-resolved abspaths to match write_fofn and the backend's bind logic).
        genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})

        # one signature file per genome, named by genome basename
        run_tool(self.capabilities,
            [
                "sourmash", "sketch", "dna",
                "-p", f"k={ksize},scaled={scaled}",
                "--from-file", fofn,
                "--outdir", sig_dir,
            ],
            logger=logger,
            log_prefix="sourmash",
            extra_mounts=genome_dirs,
        )

        matrix_csv = out_dir / "compare.csv"
        # Skip macOS AppleDouble companions ("._*") that appear on exFAT/NTFS volumes.
        sig_files = [
            p for p in (sorted(sig_dir.glob("*.sig")) + sorted(sig_dir.glob("*.sig.gz")))
            if not p.name.startswith("._")
        ]
        if not sig_files:
            raise WorkdirError(f"sourmash produced no signatures under {sig_dir}")
        # Pass signatures via --from-file, never on argv (ARG_MAX at scale).
        compare_fofn = write_fofn(sig_files, out_dir / "signatures.fofn")
        run_tool(self.capabilities,
            [
                "sourmash", "compare",
                "-k", str(ksize),
                "--csv", matrix_csv,
                "--from-file", compare_fofn,
            ],
            logger=logger,
            log_prefix="sourmash",
        )

        labels, sim = _read_compare_csv(matrix_csv)
        name_by_label = _match_labels_to_genomes(labels, genomes)
        sani = params.secondary_ani
        threshold = sani if sani <= 1.0 else sani / 100

        clusters, status = _greedy_cluster(labels, sim, name_by_label, threshold)

        rep_paths = [p for p in genomes if p.name in clusters]
        return DerepResult(
            representatives=sorted(rep_paths),
            clusters=clusters,
            genome_status=status,
        )


def _read_compare_csv(path: Path) -> tuple[list[str], npt.NDArray[np.float64]]:
    with open(path, newline="") as fo:
        labels = next(csv.reader(fo))
        # Parse the N x N body straight into a contiguous float array. At 1000s of
        # genomes this is far smaller and faster than a Python list-of-lists
        # (a 10k x 10k matrix is ~0.8 GB as float64 vs several GB of Python floats).
        matrix = np.loadtxt(fo, delimiter=",", ndmin=2)
    return labels, matrix


def _match_labels_to_genomes(labels: Sequence[str], genomes: Sequence[Path]) -> dict[str, str]:
    """Map a sourmash column label to a genome basename.

    sourmash labels are signature names (often the file path or basename). Match
    by checking which genome basename the label ends with / contains.
    """
    out: dict[str, str] = {}
    by_name = {g.name: g.name for g in genomes}
    stems = {g.stem: g.name for g in genomes}
    for label in labels:
        base = Path(label).name
        if base in by_name:
            out[label] = base
        elif Path(label).stem in stems:
            out[label] = stems[Path(label).stem]
        else:
            out[label] = base
    return out


def _greedy_cluster(
    labels: list[str],
    sim: npt.NDArray[np.float64] | list[list[float]],
    name_by_label: dict[str, str],
    threshold: float,
) -> tuple[dict[str, list[str]], dict[str, str]]:
    n = len(labels)
    # Boolean adjacency (>= threshold); the per-row vectorized ops below replace
    # the previous pure-Python O(n^2) double loops (much faster at 1000s genomes).
    adj = np.asarray(sim, dtype=float) >= threshold
    np.fill_diagonal(adj, False)
    connectivity = adj.sum(axis=1)
    # prefer well-connected genomes as representatives (stable: -connectivity, idx)
    order = np.lexsort((np.arange(n), -connectivity))

    assigned = np.full(n, -1, dtype=np.int64)  # member idx -> representative idx (-1 = free)
    reps: list[int] = []
    for i in order:
        if assigned[i] != -1:
            continue
        reps.append(int(i))
        assigned[i] = i
        # claim every still-free neighbour of i in one vectorized step
        claim = adj[i] & (assigned == -1)
        assigned[claim] = i

    clusters: dict[str, list[str]] = {}
    status: dict[str, str] = {}
    for rep_idx in reps:
        rep_name = name_by_label[labels[rep_idx]]
        clusters[rep_name] = []
        status[rep_name] = STATUS_REPRESENTATIVE
    for member_idx in range(n):
        rep_idx = int(assigned[member_idx])
        if member_idx == rep_idx:
            continue
        rep_name = name_by_label[labels[rep_idx]]
        member_name = name_by_label[labels[member_idx]]
        clusters[rep_name].append(member_name)
        status[member_name] = STATUS_CONTAINED
    return clusters, status
