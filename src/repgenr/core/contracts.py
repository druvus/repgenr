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
from dataclasses import dataclass
from pathlib import Path

CLUSTERS_TSV = "clusters.tsv"
GENOME_STATUS_TSV = "genome_status.tsv"
SELECTION_TSV = "selection.tsv"
MSA_FASTA = "msa.fasta"
CORE_SNP_FASTA = "core_snp.fasta"
TREE_NWK = "tree.nwk"
TREE2TAX_TSV = "tree2tax.tsv"
GENOMES_MAP_TSV = "genomes_map.tsv"


# Recognised genome FASTA extensions, longest-first so suffix stripping is
# unambiguous (``.fasta.gz`` before ``.fasta``). One definition shared by every
# stage and adapter that lists or names genome files.
FASTA_SUFFIXES = (".fasta.gz", ".fasta", ".fa", ".fna", ".fas")


def list_fasta(source: Path) -> list[Path]:
    """Sorted genome FASTA files directly under ``source``.

    Skips dotfiles (e.g. macOS ``._`` AppleDouble) and returns an empty list when
    ``source`` does not exist, so callers need no separate existence check.
    """
    if not source.exists():
        return []
    return sorted(
        p for p in source.iterdir()
        if not p.name.startswith(".") and p.name.endswith(FASTA_SUFFIXES)
    )


def genome_filename(family: str, genus: str, species: str, accession: str) -> str:
    """Canonical genome FASTA filename. One definition so the metadata selection,
    the genome download and every downstream stage agree on the same names.

    ``family``/``genus``/``species`` must be single tokens (no ``_``) so the name
    round-trips through :func:`parse_genome_filename`; the accession may contain
    underscores (e.g. ``GCF_000001.1``, ``NC_001802.1``).
    """
    return f"{family}_{genus}_{species}_{accession}.fasta"


def strip_fasta_suffix(name: str) -> str:
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def parse_genome_filename(name: str) -> tuple[str, str, str, str]:
    """Inverse of :func:`genome_filename`. Returns (family, genus, species,
    accession). The first three ``_``-separated tokens are the taxonomy and
    **everything after** is the accession, so accessions with underscores
    (bacterial ``GCF_x.y``, viral ``NC_x.y``) and without (viral ``MN908947.3``)
    all round-trip. A non-canonical name (< 4 tokens) yields empty taxonomy and
    the whole stem as the accession.
    """
    stem = strip_fasta_suffix(Path(name).name)
    parts = stem.split("_")
    if len(parts) < 4:
        return "", "", "", stem
    return parts[0], parts[1], parts[2], "_".join(parts[3:])


def accession_from_filename(name: str) -> str:
    """Recover the accession from a canonical genome filename or tree leaf."""
    return parse_genome_filename(name)[3]


@dataclass
class SelectionRow:
    """One selected genome: the portable hand-off from metadata to genome.

    Carries the taxonomy (so downstream taxonomy-aware stages need not re-read the
    SQLite manifest) and the canonical output filename used by the downloader.
    """

    accession: str
    family: str
    genus: str
    species: str
    is_outgroup: bool
    filename: str


def write_selection(path: Path, rows: list[SelectionRow]) -> None:
    """Write the metadata selection (accession + taxonomy + filename + outgroup flag)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t")
        writer.writerow(["accession", "family", "genus", "species", "is_outgroup", "filename"])
        for r in rows:
            writer.writerow([
                r.accession, r.family, r.genus, r.species,
                "1" if r.is_outgroup else "0", r.filename,
            ])


def read_selection(path: Path) -> list[SelectionRow]:
    """Read a selection.tsv back into SelectionRow records."""
    rows: list[SelectionRow] = []
    with open(path, newline="") as fo:
        reader = csv.DictReader(fo, delimiter="\t")
        for row in reader:
            rows.append(
                SelectionRow(
                    accession=row["accession"],
                    family=row.get("family", ""),
                    genus=row.get("genus", ""),
                    species=row.get("species", ""),
                    is_outgroup=row.get("is_outgroup", "0") == "1",
                    filename=row["filename"],
                )
            )
    return rows


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
