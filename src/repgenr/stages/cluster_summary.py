"""cluster_summary stage: one row per dereplication representative.

Summarises ``derep/clusters.tsv`` into ``derep/cluster_summary.tsv`` with the
cluster size, the species it spans (from the manifest taxonomy, else parsed
from the canonical filenames) and
the manifest CheckM quality of the keeper against its members. The dereplicate
stage writes the file itself; this stage regenerates it for an existing
working directory without rerunning the dereplicator.
"""

from __future__ import annotations

import sqlite3
from collections import Counter
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    ClusterSummaryRow,
    parse_genome_filename,
    read_clusters,
    sanitise_taxon_tokens,
    write_cluster_summary,
)
from ..core.errors import WorkdirError
from .derep_keeper import quality_score

Quality = Mapping[str, tuple[float, float]]
# filename -> (genus, species) as recorded in the manifest or selection.tsv.
Taxonomy = Mapping[str, tuple[str, str]]

# Species names listed in the ``species`` column before "+N more".
SPECIES_LIST_MAX = 5


@dataclass
class ClusterSummaryParams:
    pass


def summarise_clusters(
    clusters: Mapping[str, list[str]], quality: Quality, taxonomy: Taxonomy | None = None
) -> list[ClusterSummaryRow]:
    """Build the summary rows, largest cluster first, then by representative name.

    ``taxonomy`` gives the species of each genome; a genome it lacks (or
    holds without a species) falls back to its canonical filename.
    """
    taxonomy = taxonomy or {}
    rows = [_summarise(rep, members, quality, taxonomy) for rep, members in clusters.items()]
    rows.sort(key=lambda r: (-r.n_members, r.representative))
    return rows


def _taxon(name: str, taxonomy: Taxonomy) -> tuple[str, str]:
    """(genus, species) tokens of a genome; species is blank when unknown.

    The manifest supplies the taxonomy; a genome without one, or whose
    species is blank, falls back to its canonical filename. A manifest
    species with a blank genus (``ingest --selection`` requires neither)
    takes the genus from the species when it is a binomial, else from the
    filename, so it matches canonical members.
    """
    file_genus, file_species = parse_genome_filename(name)[1:3]
    genus, species = taxonomy.get(name, ("", ""))
    if species:
        if not genus and " " in species.strip():
            # A binomial without a genus column: take the genus from it.
            genus, species = species.strip().split(" ", 1)
        genus = genus or file_genus
        _, genus, species = sanitise_taxon_tokens("", genus, species)
        if species:
            return genus, species
    return file_genus, file_species


def _species_column(rep: str, others: list[str], taxonomy: Taxonomy) -> tuple[int, str]:
    """Distinct species and the capped, comma-separated list of their names.

    The keeper's species comes first, then the others by member count and
    name. A genome without a species (non-canonical filename, no taxonomy)
    adds none. An epithet shared by two genera is written with its genus.
    """
    taxa = [_taxon(n, taxonomy) for n in (rep, *others)]
    # An epithet known without its genus (manifest row with a blank genus and
    # a non-canonical filename) belongs to the one genus the cluster holds it
    # under, when there is exactly one.
    genera = {sp: {g for g, s in taxa if s == sp and g} for _, sp in taxa if sp}
    taxa = [
        (g or next(iter(genera[sp])) if sp and len(genera[sp]) == 1 else g, sp) for g, sp in taxa
    ]
    counts = Counter(t for t in taxa if t[1])
    keeper = taxa[0]
    ordered = sorted(counts, key=lambda t: (t != keeper, -counts[t], t[1], t[0]))
    genera_per_epithet = Counter(sp for _, sp in counts)
    names = [sp if genera_per_epithet[sp] == 1 or not g else f"{g} {sp}" for g, sp in ordered]
    if len(names) > SPECIES_LIST_MAX:
        names = [*names[:SPECIES_LIST_MAX], f"+{len(names) - SPECIES_LIST_MAX} more"]
    return len(counts), ",".join(names)


def _summarise(
    rep: str, members: list[str], quality: Quality, taxonomy: Taxonomy
) -> ClusterSummaryRow:
    others = [m for m in members if m != rep]
    n_species, species = _species_column(rep, others, taxonomy)

    rep_q = quality.get(rep)
    scored = [(m, quality[m]) for m in others if m in quality]
    best_member = ""
    candidates = [(rep, rep_q)] if rep_q is not None else []
    candidates.extend(scored)
    if candidates:
        # ``max`` keeps the first of equal scores, so a tie leaves the keeper.
        best_member = max(candidates, key=lambda c: quality_score(*c[1]))[0]

    return ClusterSummaryRow(
        representative=rep,
        n_members=len(others),
        n_species=n_species,
        species=species,
        rep_completeness=None if rep_q is None else rep_q[0],
        rep_contamination=None if rep_q is None else rep_q[1],
        member_max_completeness=max((q[0] for _, q in scored), default=None),
        member_min_contamination=min((q[1] for _, q in scored), default=None),
        best_member=best_member,
    )


def taxonomy_lookup(ctx: WorkdirContext) -> dict[str, tuple[str, str]]:
    """Map each genome filename to its manifest (genus, species)."""
    try:
        records = ctx.manifest.all_genomes(include_outgroup=True)
    except (sqlite3.OperationalError, OSError):
        # No manifest (data-channel path, tests): the filenames supply species.
        return {}
    return {r.filename: (r.genus or "", r.species or "") for r in records if r.filename}


def run(ctx: WorkdirContext, params: ClusterSummaryParams) -> Path:
    logger = ctx.logger
    clusters_file = ctx.derep_dir / CLUSTERS_TSV
    if not clusters_file.exists():
        raise WorkdirError(f"Missing {clusters_file}. Run the dereplicate stage first.")
    clusters = read_clusters(clusters_file)
    if not clusters:
        logger.warning("%s lists no clusters; the summary holds only its header", clusters_file)
    from .dereplicate import quality_lookup

    quality = quality_lookup(ctx)
    if not quality:
        logger.info("No assembly quality in the manifest; quality columns are left blank")
    rows = summarise_clusters(clusters, quality, taxonomy_lookup(ctx))
    out = ctx.derep_dir / CLUSTER_SUMMARY_TSV
    write_cluster_summary(out, rows)
    ctx.config.record_stage(
        "cluster_summary",
        params={**asdict(params), "clusters": len(rows)},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Summarised %d clusters into %s", len(rows), out)
    return out
