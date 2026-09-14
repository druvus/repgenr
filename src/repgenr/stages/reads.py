"""reads stage: select sequencing runs from ENA/SRA and write ``reads.tsv``.

The first stage of the reads chain (``reads -> assemble -> dereplicate ->
phylo -> tree2tax``). Runs are found by a taxon query (``tax_tree`` over the
resolved taxid) or by explicit run, sample or study accessions, filtered by
platform and size, reduced to the best run per sample, and labelled with the
NCBI family, genus and species of their taxid. Network only; nothing is
downloaded here.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..core import ena
from ..core.context import WorkdirContext
from ..core.contracts import READS_TSV, ReadRow, sanitise_taxon_tokens, write_reads
from ..core.errors import UserInputError
from ..core.ncbi_taxonomy import get_taxon_data_from_entrez

PLATFORMS = ("any", "illumina", "ont", "pacbio")
_ENA_PLATFORM = {"illumina": "ILLUMINA", "ont": "OXFORD_NANOPORE", "pacbio": "PACBIO_SMRT"}
_LONG_READ = {"OXFORD_NANOPORE", "PACBIO_SMRT"}


@dataclass
class ReadsParams:
    target_family: str | None = None
    target_genus: str | None = None
    target_species: str | None = None
    accessions: list[str] = field(default_factory=list)
    accession_file: str | None = None
    platform: str = "any"  # any | illumina | ont | pacbio
    max_runs: int | None = None
    min_bases: int = 0
    # Drop runs above this many bases: a whole-host library (tens of Gb for a
    # small endosymbiont) would assemble into a host-dominated genome.
    max_bases: int | None = None
    # ENA library_selection values to drop (case-insensitive). MDA (multiple
    # displacement amplification) gives chimeric, uneven assemblies.
    drop_selection: list[str] = field(default_factory=lambda: ["MDA"])
    # Keep the best run of each sample: a long-read run with enough bases, else
    # the largest run.
    one_per_sample: bool = True


def run(ctx: WorkdirContext, params: ReadsParams) -> int:
    logger = ctx.logger
    accessions = [*params.accessions, *_read_accession_file(params.accession_file)]
    target = params.target_species or params.target_genus or params.target_family
    if not accessions and not target:
        raise UserInputError(
            "Select runs by taxon (-tf/-tg/-ts) or by accession (--accession, --accession-file)."
        )

    taxid: str | None = None
    records: list[dict] = []
    if target:
        hit = ena.resolve_taxon(target)
        taxid = hit.taxid
        logger.info(
            "Resolved %r to %s (%s, taxid %s)", target, hit.scientific_name, hit.rank, taxid
        )
        records += ena.search_runs(ena.taxon_query(taxid))
    if accessions:
        records += ena.search_runs(ena.accession_query(accessions))
    rows = _dedupe(ena.to_read_rows(records))
    candidates = len(rows)
    logger.info("ENA returned %d whole-genome sequencing runs", candidates)

    rows = _filter(rows, params, logger)
    if params.one_per_sample:
        rows = _best_per_sample(rows)
    rows.sort(key=lambda r: -r.bases)
    if params.max_runs is not None:
        rows = rows[: params.max_runs]
    if not rows:
        raise UserInputError(
            f"No sequencing runs selected ({candidates} candidates before the platform, size "
            "and per-sample filters). Loosen --platform/--min-bases/--drop-selection or "
            "check the taxon."
        )

    rows = _label(rows, logger)
    write_reads(ctx.workdir / READS_TSV, rows)
    ctx.config.record_stage(
        "reads",
        params={
            **asdict(params),
            "taxid": taxid,
            "candidates": candidates,
            "selected": len(rows),
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Selected %d sequencing runs; wrote %s", len(rows), READS_TSV)
    return len(rows)


def _read_accession_file(path: str | None) -> list[str]:
    if path is None:
        return []
    file = Path(path)
    if not file.is_file():
        raise UserInputError(f"--accession-file {path} is not a file.")
    return [
        line.strip()
        for line in file.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]


def _dedupe(rows: list[ReadRow]) -> list[ReadRow]:
    seen: set[str] = set()
    out = []
    for row in rows:
        if row.run_accession not in seen:
            seen.add(row.run_accession)
            out.append(row)
    return out


def _filter(rows: list[ReadRow], params: ReadsParams, logger: logging.Logger) -> list[ReadRow]:
    wanted = _ENA_PLATFORM.get(params.platform)
    kept = [
        r for r in rows if (wanted is None or r.platform == wanted) and r.bases >= params.min_bases
    ]
    if params.max_bases is not None:
        over = [r for r in kept if r.bases > params.max_bases]
        if over:
            logger.info(
                "Dropping %d run(s) above --max-bases %d (largest %d bp, %s): likely "
                "whole-host libraries",
                len(over),
                params.max_bases,
                max(r.bases for r in over),
                max(over, key=lambda r: r.bases).run_accession,
            )
        kept = [r for r in kept if r.bases <= params.max_bases]
    drop = {s.upper() for s in params.drop_selection}
    if drop:
        amplified = [r for r in kept if r.library_selection.upper() in drop]
        if amplified:
            logger.info(
                "Dropping %d run(s) by library selection (%s): %s",
                len(amplified),
                ", ".join(sorted(drop)),
                ", ".join(r.run_accession for r in amplified[:5])
                + (" ..." if len(amplified) > 5 else ""),
            )
        kept = [r for r in kept if r.library_selection.upper() not in drop]
    return kept


# A long-read run is preferred over the sample's short-read runs only when it
# carries enough sequence to assemble on its own: at least this many bases,
# and at least this fraction of the largest short-read run. Below that it is
# usually a scaffolding or test run next to the real data (seen on Wolbachia:
# a 45 kb PacBio run and a 338 Mb ONT run beside 9 Gb and 11 Gb Illumina runs).
LONG_READ_MIN_BASES = 100_000_000
LONG_READ_MIN_FRACTION = 0.1


def _best_per_sample(rows: list[ReadRow]) -> list[ReadRow]:
    by_sample: dict[str, list[ReadRow]] = {}
    for row in rows:
        by_sample.setdefault(row.biosample or row.run_accession, []).append(row)
    return [_best_run(runs) for runs in by_sample.values()]


def _best_run(runs: list[ReadRow]) -> ReadRow:
    short_max = max((r.bases for r in runs if r.platform not in _LONG_READ), default=0)

    def rank(r: ReadRow) -> tuple[bool, int]:
        usable_long = (
            r.platform in _LONG_READ
            and r.bases >= LONG_READ_MIN_BASES
            and r.bases >= LONG_READ_MIN_FRACTION * short_max
        )
        return (usable_long, r.bases)

    return max(runs, key=rank)


def _label(rows: list[ReadRow], logger: logging.Logger) -> list[ReadRow]:
    """Fill the family/genus/species tokens from each run's NCBI taxid."""
    taxids = sorted({r.taxid for r in rows if r.taxid})
    data, missing, _alts = get_taxon_data_from_entrez(taxids, logger)
    if missing:
        logger.warning("No NCBI lineage for %d taxid(s); their tokens are 'unknown'", len(missing))
    labelled = []
    for row in rows:
        entry = data.get(row.taxid) or {}
        taxdata = entry.get("taxdata") or {}
        names = [
            (taxdata.get(level) or {}).get("name") or "" for level in ("family", "genus", "species")
        ]
        family, genus, species = sanitise_taxon_tokens(*names)
        labelled.append(
            ReadRow(
                **{
                    **asdict(row),
                    "family": family or "unknown",
                    "genus": genus or "unknown",
                    "species": species or "unknown",
                }
            )
        )
    return labelled
