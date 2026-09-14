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
# Segment-grouped viral isolates (vgenome --group-segments): the synthetic
# isolate token used as the genome's accession -> its member accessions, so
# tree2tax can list the real accessions under the isolate's leaf.
SEGMENTS_TSV = "segments.tsv"
# Accessions the genome stage requested but NCBI returned nothing for;
# the completeness guard excuses them (core.integrity).
MISSING_ACCESSIONS_TXT = "missing_accessions.txt"
# The reads-to-assembly entry path: the sequencing runs the reads stage
# selected, the per-assembly metrics, and the runs the assemble stage could not
# turn into a genome (excused by the completeness guard, like missing
# accessions).
READS_TSV = "reads.tsv"
ASSEMBLY_STATS_TSV = "assembly_stats.tsv"
EXCUSED_RUNS_TSV = "excused_runs.tsv"


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


def sanitise_taxon_tokens(family: str, genus: str, species: str) -> tuple[str, str, str]:
    """Turn taxonomy names into the single tokens a canonical filename holds.

    The species drops its genus prefix (``Francisella tularensis`` ->
    ``tularensis``); every token loses its spaces and turns underscores into
    hyphens, so :func:`parse_genome_filename` splits the name back on ``_``.
    One rule for every producer (GTDB TSV and API, the reads path), so genomes
    of one species carry one token whatever selected them.
    """

    def clean(value: str) -> str:
        return value.replace(" ", "").replace("_", "-")

    species = species.replace(genus, "") if genus else species
    return clean(family), clean(genus), clean(species)


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


@dataclass(frozen=True)
class ReadRow:
    """One sequencing run selected by the reads stage (the reads.tsv contract).

    ``fastq_urls``/``fastq_md5``/``fastq_bytes`` are parallel per-file tuples,
    empty when the archive holds no FASTQ mirror for the run. A URL may also
    be a local path, which the assemble stage copies instead of downloading.
    """

    run_accession: str
    biosample: str
    bioproject: str
    organism: str
    taxid: str
    platform: str
    instrument_model: str
    layout: str
    bases: int
    read_count: int
    # Sanitised filename tokens resolved from the taxid by the reads stage.
    family: str = ""
    genus: str = ""
    species: str = ""
    fastq_urls: tuple[str, ...] = ()
    fastq_md5: tuple[str, ...] = ()
    fastq_bytes: tuple[int, ...] = ()
    # ENA's library_selection (RANDOM, MDA, PCR, ...); MDA marks whole-genome
    # amplification, which assembles poorly.
    library_selection: str = ""


_READS_COLUMNS = [
    "run_accession",
    "biosample",
    "bioproject",
    "organism",
    "taxid",
    "platform",
    "instrument_model",
    "layout",
    "bases",
    "read_count",
    "family",
    "genus",
    "species",
    "fastq_urls",
    "fastq_md5",
    "fastq_bytes",
    "library_selection",
]


def write_reads(path: Path, rows: list[ReadRow]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(_READS_COLUMNS)
        for r in rows:
            writer.writerow(
                [
                    r.run_accession,
                    r.biosample,
                    r.bioproject,
                    r.organism,
                    r.taxid,
                    r.platform,
                    r.instrument_model,
                    r.layout,
                    r.bases,
                    r.read_count,
                    r.family,
                    r.genus,
                    r.species,
                    ";".join(r.fastq_urls),
                    ";".join(r.fastq_md5),
                    ";".join(str(b) for b in r.fastq_bytes),
                    r.library_selection,
                ]
            )


def read_reads(path: Path) -> list[ReadRow]:
    rows: list[ReadRow] = []
    with open(path, encoding="utf-8", newline="") as fo:
        for rec in csv.DictReader(fo, delimiter="\t"):
            split = lambda s: tuple(x for x in s.split(";") if x)  # noqa: E731
            rows.append(
                ReadRow(
                    run_accession=rec["run_accession"],
                    biosample=rec["biosample"],
                    bioproject=rec["bioproject"],
                    organism=rec["organism"],
                    taxid=rec["taxid"],
                    platform=rec["platform"],
                    instrument_model=rec["instrument_model"],
                    layout=rec["layout"],
                    bases=int(rec["bases"] or 0),
                    read_count=int(rec["read_count"] or 0),
                    family=rec["family"],
                    genus=rec["genus"],
                    species=rec["species"],
                    fastq_urls=split(rec["fastq_urls"]),
                    fastq_md5=split(rec["fastq_md5"]),
                    fastq_bytes=tuple(int(b) for b in split(rec["fastq_bytes"])),
                    # Absent from tables written before the column existed.
                    library_selection=rec.get("library_selection", "") or "",
                )
            )
    return rows


@dataclass(frozen=True)
class AssemblyStatsRow:
    """Per-assembly metrics and labels (assembly_stats.tsv), one row per genome."""

    run_accession: str
    filename: str
    assembler: str
    n_contigs: int
    total_length: int
    n50: int
    largest_contig: int
    est_coverage: float | None = None
    completeness: float | None = None
    contamination: float | None = None
    ncbi_taxonomy: str = ""
    gtdb_taxonomy: str = ""
    label_source: str = "metadata"
    taxonomy_flag: str = ""


_ASSEMBLY_STATS_COLUMNS = [
    "run_accession",
    "filename",
    "assembler",
    "n_contigs",
    "total_length",
    "n50",
    "largest_contig",
    "est_coverage",
    "completeness",
    "contamination",
    "ncbi_taxonomy",
    "gtdb_taxonomy",
    "label_source",
    "taxonomy_flag",
]


def write_assembly_stats(path: Path, rows: list[AssemblyStatsRow]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(_ASSEMBLY_STATS_COLUMNS)
        for r in rows:
            writer.writerow(
                [
                    r.run_accession,
                    r.filename,
                    r.assembler,
                    r.n_contigs,
                    r.total_length,
                    r.n50,
                    r.largest_contig,
                    "" if r.est_coverage is None else f"{r.est_coverage:.2f}",
                    "" if r.completeness is None else f"{r.completeness:.2f}",
                    "" if r.contamination is None else f"{r.contamination:.2f}",
                    r.ncbi_taxonomy,
                    r.gtdb_taxonomy,
                    r.label_source,
                    r.taxonomy_flag,
                ]
            )


def read_assembly_stats(path: Path) -> list[AssemblyStatsRow]:
    rows: list[AssemblyStatsRow] = []
    with open(path, encoding="utf-8", newline="") as fo:
        for rec in csv.DictReader(fo, delimiter="\t"):
            rows.append(
                AssemblyStatsRow(
                    run_accession=rec["run_accession"],
                    filename=rec["filename"],
                    assembler=rec["assembler"],
                    n_contigs=int(rec["n_contigs"]),
                    total_length=int(rec["total_length"]),
                    n50=int(rec["n50"]),
                    largest_contig=int(rec["largest_contig"]),
                    est_coverage=_opt_float(rec["est_coverage"]),
                    completeness=_opt_float(rec["completeness"]),
                    contamination=_opt_float(rec["contamination"]),
                    ncbi_taxonomy=rec["ncbi_taxonomy"],
                    gtdb_taxonomy=rec["gtdb_taxonomy"],
                    label_source=rec["label_source"],
                    taxonomy_flag=rec["taxonomy_flag"],
                )
            )
    return rows


@dataclass(frozen=True)
class ExcusedRun:
    """A selected run that produced no genome: which step gave up and why."""

    run_accession: str
    step: str  # fetch | assemble | qc | classify
    reason: str


def write_excused_runs(path: Path, rows: list[ExcusedRun]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["run_accession", "step", "reason"])
        for r in rows:
            writer.writerow([r.run_accession, r.step, r.reason])


def read_excused_runs(path: Path) -> list[ExcusedRun]:
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)
        return [ExcusedRun(row[0], row[1], row[2]) for row in reader if len(row) >= 3]


def write_segments(path: Path, members: dict[str, list[str]]) -> None:
    """Write isolate token -> member accession rows (one row per member)."""
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["isolate", "accession"])
        for token, accessions in members.items():
            for accession in accessions:
                writer.writerow([token, accession])


def read_segments(path: Path) -> dict[str, list[str]]:
    """Read ``segments.tsv`` back into isolate token -> member accessions."""
    members: dict[str, list[str]] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                members.setdefault(row[0], []).append(row[1])
    return members


def write_genomes_map(path: Path, mapping: list[tuple[str, str]]) -> None:
    """Write accession -> leaf rows."""
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        for accession, leaf in mapping:
            writer.writerow([accession, leaf])
