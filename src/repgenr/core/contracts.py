"""Canonical inter-stage file contracts.

Each stage publishes a documented, validated artifact set that the next stage
consumes. These are internal engineering contracts (no legacy backward-compat
obligation), so the names and layout are chosen fresh:

    derep/representatives/      representative genome FASTAs
    derep/clusters.tsv          representative<TAB>member (member==representative for self)
    derep/genome_status.tsv     genome<TAB>status(representative|contained|fail_qc)
    derep/cluster_summary.tsv   one row per representative: size, species, quality
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
import os
from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any

CLUSTERS_TSV = "clusters.tsv"
GENOME_STATUS_TSV = "genome_status.tsv"
CLUSTER_SUMMARY_TSV = "cluster_summary.tsv"
SELECTION_TSV = "selection.tsv"
MSA_FASTA = "msa.fasta"
CORE_SNP_FASTA = "core_snp.fasta"
TREE_NWK = "tree.nwk"
TREE2TAX_TSV = "tree2tax.tsv"
GENOMES_MAP_TSV = "genomes_map.tsv"
# Accessions the genome stage requested but NCBI returned nothing for;
# the completeness guard excuses them (core.integrity).
MISSING_ACCESSIONS_TXT = "missing_accessions.txt"


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
        p
        for p in source.iterdir()
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
    completeness: float | None = None
    contamination: float | None = None


def _tsv_writer(fo: IO) -> Any:
    """A TSV writer with Unix line endings. ``csv.writer`` defaults to
    ``\\r\\n``, which leaves a stray ``\\r`` on the last column of every row
    for awk/cut consumers of the contract files."""
    return csv.writer(fo, delimiter="\t", lineterminator="\n")


@contextmanager
def atomic_replace(
    path: Path, *, mode: str = "w", encoding: str | None = "utf-8", newline: str | None = None
) -> Iterator[IO]:
    """Open a temp sibling of ``path`` for writing; publish it atomically.

    On clean exit the temp file replaces ``path`` (``os.replace``); on any
    exception the temp is removed and a previous ``path`` is left untouched, so
    a crash mid-write can never truncate a deliverable. Same pattern as
    ``Config.save``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        with open(tmp, mode, encoding=encoding, newline=newline) as fo:
            yield fo
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


@contextmanager
def atomic_path(path: Path) -> Iterator[Path]:
    """Yield a temp sibling path for tools/copies that need a filename.

    On clean exit the temp (if written) replaces ``path``; on exception it is
    removed and the previous ``path`` survives.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        yield tmp
        if tmp.exists():
            os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def write_selection(path: Path, rows: list[SelectionRow]) -> None:
    """Write the metadata selection (accession + taxonomy + filename + outgroup flag)."""
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(
            [
                "accession",
                "family",
                "genus",
                "species",
                "is_outgroup",
                "filename",
                "completeness",
                "contamination",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    r.accession,
                    r.family,
                    r.genus,
                    r.species,
                    "1" if r.is_outgroup else "0",
                    r.filename,
                    "" if r.completeness is None else f"{r.completeness:.2f}",
                    "" if r.contamination is None else f"{r.contamination:.2f}",
                ]
            )


def _opt_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def read_selection(path: Path) -> list[SelectionRow]:
    """Read a selection.tsv back into SelectionRow records."""
    rows: list[SelectionRow] = []
    with open(path, encoding="utf-8", newline="") as fo:
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
                    completeness=_opt_float(row.get("completeness")),
                    contamination=_opt_float(row.get("contamination")),
                )
            )
    return rows


def write_clusters(path: Path, clusters: dict[str, list[str]]) -> None:
    """Write representative -> members. Each representative also lists itself."""
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["representative", "member"])
        for rep, members in clusters.items():
            writer.writerow([rep, rep])
            for member in members:
                if member != rep:
                    writer.writerow([rep, member])


def read_clusters(path: Path) -> dict[str, list[str]]:
    """Read representative -> members (members exclude the representative itself)."""
    clusters: dict[str, list[str]] = defaultdict(list)
    with open(path, encoding="utf-8", newline="") as fo:
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
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["genome", "status"])
        for genome, value in sorted(status.items()):
            writer.writerow([genome, value])


def read_genome_status(path: Path) -> dict[str, str]:
    """Read genome -> status. A missing file yields an empty mapping.

    The counterpart of :func:`write_genome_status`, so a step that consumes a
    contract directory (Nextflow scatter-gather) recovers the statuses of
    genomes that belong to no cluster -- ``fail_qc`` genomes would otherwise be
    invisible to the consumer.
    """
    status: dict[str, str] = {}
    if not path.exists():
        return status
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 2:
                continue
            status[row[0]] = row[1]
    return status


@dataclass
class ClusterSummaryRow:
    """One dereplication cluster seen from its representative.

    ``n_members`` excludes the representative. ``species`` lists the distinct
    species across representative and members, the representative's first.
    Quality columns are ``None`` when the manifest carried no CheckM values;
    the ``member_*`` extremes span scored members only. ``best_member`` is the
    highest-scoring genome in the cluster (representative included) and equals
    ``representative`` when the keeper is already the best; it is empty when no
    genome in the cluster is scored.
    """

    representative: str
    n_members: int
    n_species: int
    species: str
    rep_completeness: float | None = None
    rep_contamination: float | None = None
    member_max_completeness: float | None = None
    member_min_contamination: float | None = None
    best_member: str = ""


_CLUSTER_SUMMARY_COLUMNS = (
    "representative",
    "n_members",
    "n_species",
    "species",
    "rep_completeness",
    "rep_contamination",
    "member_max_completeness",
    "member_min_contamination",
    "best_member",
)


def write_cluster_summary(path: Path, rows: list[ClusterSummaryRow]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(_CLUSTER_SUMMARY_COLUMNS)
        for r in rows:
            writer.writerow(
                [
                    r.representative,
                    r.n_members,
                    r.n_species,
                    r.species,
                    _fmt_opt(r.rep_completeness),
                    _fmt_opt(r.rep_contamination),
                    _fmt_opt(r.member_max_completeness),
                    _fmt_opt(r.member_min_contamination),
                    r.best_member,
                ]
            )


def _fmt_opt(value: float | None) -> str:
    return "" if value is None else repr(value)


def read_cluster_summary(path: Path) -> list[ClusterSummaryRow]:
    rows: list[ClusterSummaryRow] = []
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo, delimiter="\t")
        for rec in reader:
            rows.append(
                ClusterSummaryRow(
                    representative=rec["representative"],
                    n_members=int(rec["n_members"]),
                    n_species=int(rec["n_species"]),
                    species=rec["species"],
                    rep_completeness=_opt_float(rec.get("rep_completeness")),
                    rep_contamination=_opt_float(rec.get("rep_contamination")),
                    member_max_completeness=_opt_float(rec.get("member_max_completeness")),
                    member_min_contamination=_opt_float(rec.get("member_min_contamination")),
                    best_member=rec.get("best_member", ""),
                )
            )
    return rows


def write_tree2tax(path: Path, edges: list[tuple[str, str]]) -> None:
    """Write child -> parent edges (FlexTaxD), de-duplicated, order preserved."""
    seen: set[tuple[str, str]] = set()
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["child", "parent"])
        for child, parent in edges:
            if (child, parent) in seen:
                continue
            seen.add((child, parent))
            writer.writerow([child, parent])


def write_genomes_map(path: Path, mapping: list[tuple[str, str]]) -> None:
    """Write accession -> leaf rows."""
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        for accession, leaf in mapping:
            writer.writerow([accession, leaf])
