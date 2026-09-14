"""ENA Portal client: sequencing-run discovery and taxon resolution.

ENA mirrors SRA, needs no key, and answers one JSON query with the run
metadata and the FASTQ locations and checksums, which is why the reads stage
reads from it. Network access only; the JSON is easy to freeze as a fixture.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from urllib.parse import quote

from . import http
from .contracts import ReadRow
from .errors import UserInputError

PORTAL_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
TAXONOMY_URL = "https://www.ebi.ac.uk/ena/taxonomy/rest"

RUN_FIELDS = (
    "run_accession",
    "sample_accession",
    "study_accession",
    "tax_id",
    "scientific_name",
    "instrument_platform",
    "instrument_model",
    "library_layout",
    "library_strategy",
    "library_source",
    "base_count",
    "read_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "first_public",
)

# Whole-genome sequencing of the organism itself; excludes amplicons,
# transcriptomes and metagenomes.
_WGS_FILTER = 'library_strategy="WGS" AND library_source="GENOMIC"'

# Accession prefixes and the portal field they are filtered on.
_ACCESSION_FIELDS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"^[SED]RR\d+$"), "run_accession"),
    (re.compile(r"^[SED]RS\d+$"), "secondary_sample_accession"),
    (re.compile(r"^SAM[NED]A?\d+$"), "sample_accession"),
    (re.compile(r"^PRJ[NED][A-Z]\d+$"), "study_accession"),
    (re.compile(r"^[SED]RP\d+$"), "secondary_study_accession"),
)


@dataclass(frozen=True)
class TaxonHit:
    taxid: str
    scientific_name: str
    rank: str


def taxon_query(taxid: str) -> str:
    """Every WGS run under a taxon subtree."""
    return f"tax_tree({taxid}) AND {_WGS_FILTER}"


def accession_query(accessions: list[str]) -> str:
    """Runs named directly or through their sample or study accessions."""
    terms = []
    for acc in accessions:
        for pattern, field in _ACCESSION_FIELDS:
            if pattern.match(acc):
                terms.append(f'{field}="{acc}"')
                break
        else:
            raise UserInputError(
                f"Unrecognised accession {acc!r}: expected a run (SRR/ERR/DRR), a "
                "sample (SAMN/SAME/SAMD or SRS/ERS/DRS), or a study (PRJNA/PRJEB/PRJDB "
                "or SRP/ERP/DRP)."
            )
    return " OR ".join(terms)


def resolve_taxon(name: str) -> TaxonHit:
    """Resolve a taxon name (synonyms included) to one taxid, or explain why not."""
    hits = http.get_json(f"{TAXONOMY_URL}/any-name/{quote(name)}")
    if not hits:
        raise UserInputError(f"Nothing in the ENA taxonomy matches {name!r}.")
    if len(hits) > 1:
        listing = "; ".join(
            f"{h.get('scientificName')} ({h.get('rank')}, taxid {h.get('taxId')})" for h in hits
        )
        raise UserInputError(
            f"{name!r} matches {len(hits)} taxa: {listing}. Use the scientific name of one."
        )
    hit = hits[0]
    return TaxonHit(str(hit["taxId"]), str(hit["scientificName"]), str(hit.get("rank", "")))


def search_runs(query: str, *, page: int = 10000) -> list[dict]:
    """All read_run records matching ``query``, paged through the portal."""
    records: list[dict] = []
    offset = 0
    while True:
        batch = http.get_json(
            PORTAL_URL,
            params={
                "result": "read_run",
                "query": query,
                "fields": ",".join(RUN_FIELDS),
                "format": "json",
                "limit": page,
                "offset": offset,
            },
        )
        records.extend(batch)
        if len(batch) < page:
            return records
        offset += page


def _split(value: str | None) -> tuple[str, ...]:
    return tuple(part for part in (value or "").split(";") if part)


def _https(url: str) -> str:
    """ENA lists FTP paths without a scheme; the same host serves them over HTTPS."""
    if url.startswith(("http://", "https://")):
        return url
    return "https://" + url.removeprefix("ftp://")


def to_read_rows(records: list[dict]) -> list[ReadRow]:
    """Normalise portal records into the reads contract (taxonomy tokens unset)."""
    rows = []
    for rec in records:
        rows.append(
            ReadRow(
                run_accession=rec["run_accession"],
                biosample=rec.get("sample_accession", "") or "",
                bioproject=rec.get("study_accession", "") or "",
                organism=rec.get("scientific_name", "") or "",
                taxid=str(rec.get("tax_id", "") or ""),
                platform=rec.get("instrument_platform", "") or "",
                instrument_model=rec.get("instrument_model", "") or "",
                layout=rec.get("library_layout", "") or "",
                bases=int(rec.get("base_count") or 0),
                read_count=int(rec.get("read_count") or 0),
                fastq_urls=tuple(_https(u) for u in _split(rec.get("fastq_ftp"))),
                fastq_md5=_split(rec.get("fastq_md5")),
                fastq_bytes=tuple(int(b) for b in _split(rec.get("fastq_bytes"))),
            )
        )
    return rows
