"""Canonical inter-stage file contracts.

Each stage publishes a documented, validated artifact set that the next stage
consumes. These are internal engineering contracts (no legacy backward-compat
obligation), so the names and layout are chosen fresh:

    derep/representatives/      representative genome FASTAs
    derep/clusters.tsv          representative<TAB>member (member==representative for self)
    derep/genome_status.tsv     genome<TAB>status(representative|contained|fail_qc)
    align/msa.fasta             multiple sequence alignment (aligner output)
    snp/core_snp.fasta          variant-site alignment (snp typer output)
    tree/tree.nwk               Newick tree
    tree2tax.tsv                child<TAB>parent (FlexTaxD)
    genomes_map.tsv             accession<TAB>leaf

This module owns the writers/readers so producers (adapters) and consumers
(downstream stages) agree on one place.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

CLUSTERS_TSV = "clusters.tsv"
GENOME_STATUS_TSV = "genome_status.tsv"
MSA_FASTA = "msa.fasta"
CORE_SNP_FASTA = "core_snp.fasta"
TREE_NWK = "tree.nwk"
TREE2TAX_TSV = "tree2tax.tsv"
GENOMES_MAP_TSV = "genomes_map.tsv"


def write_clusters(path: Path, clusters: dict[str, list[str]]) -> None:
    """Write representative -> members. Each representative also lists itself."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t")
        writer.writerow(["representative", "member"])
        for rep, members in clusters.items():
            writer.writerow([rep, rep])
            for member in members:
                if member != rep:
                    writer.writerow([rep, member])


def read_clusters(path: Path) -> dict[str, list[str]]:
    """Read representative -> members (members exclude the representative itself)."""
    clusters: dict[str, list[str]] = defaultdict(list)
    with open(path, newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 2:
                continue
            rep, member = row[0], row[1]
            clusters.setdefault(rep, [])
            if member != rep:
                clusters[rep].append(member)
    return dict(clusters)


def write_genome_status(path: Path, status: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t")
        writer.writerow(["genome", "status"])
        for genome, value in sorted(status.items()):
            writer.writerow([genome, value])


def write_tree2tax(path: Path, edges: list[tuple[str, str]]) -> None:
    """Write child -> parent edges (FlexTaxD), de-duplicated, order preserved."""
    seen: set[tuple[str, str]] = set()
    with open(path, "w", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t")
        writer.writerow(["child", "parent"])
        for child, parent in edges:
            if (child, parent) in seen:
                continue
            seen.add((child, parent))
            writer.writerow([child, parent])


def write_genomes_map(path: Path, mapping: list[tuple[str, str]]) -> None:
    """Write accession -> leaf rows."""
    with open(path, "w", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t")
        for accession, leaf in mapping:
            writer.writerow([accession, leaf])
