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
from collections.abc import Sequence
from pathlib import Path

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
        required_binaries=(BinarySpec("sourmash", version_args=("--version",)),),
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
        )

        matrix_csv = out_dir / "compare.csv"
        # Skip macOS AppleDouble companions ("._*") that appear on exFAT/NTFS volumes.
        sig_files = [
            p for p in (sorted(sig_dir.glob("*.sig")) + sorted(sig_dir.glob("*.sig.gz")))
            if not p.name.startswith("._")
        ]
        if not sig_files:
            raise WorkdirError(f"sourmash produced no signatures under {sig_dir}")
        run_tool(self.capabilities, 
            [
                "sourmash", "compare",
                "-k", str(ksize),
                "--csv", matrix_csv,
                *sig_files,
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


def _read_compare_csv(path: Path) -> tuple[list[str], list[list[float]]]:
    with open(path, newline="") as fo:
        reader = csv.reader(fo)
        labels = next(reader)
        matrix = [[float(x) for x in row] for row in reader]
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
    sim: list[list[float]],
    name_by_label: dict[str, str],
    threshold: float,
) -> tuple[dict[str, list[str]], dict[str, str]]:
    n = len(labels)
    # connectivity = count of neighbours above threshold (prefer well-connected reps)
    connectivity = [
        sum(1 for j in range(n) if j != i and sim[i][j] >= threshold) for i in range(n)
    ]
    order = sorted(range(n), key=lambda i: connectivity[i], reverse=True)

    assigned: dict[int, int] = {}  # member idx -> representative idx
    reps: list[int] = []
    for i in order:
        if i in assigned:
            continue
        reps.append(i)
        assigned[i] = i
        for j in range(n):
            if j not in assigned and sim[i][j] >= threshold:
                assigned[j] = i

    clusters: dict[str, list[str]] = {}
    status: dict[str, str] = {}
    for rep_idx in reps:
        rep_name = name_by_label[labels[rep_idx]]
        clusters[rep_name] = []
        status[rep_name] = STATUS_REPRESENTATIVE
    for member_idx, rep_idx in assigned.items():
        if member_idx == rep_idx:
            continue
        rep_name = name_by_label[labels[rep_idx]]
        member_name = name_by_label[labels[member_idx]]
        clusters[rep_name].append(member_name)
        status[member_name] = STATUS_CONTAINED
    return clusters, status
