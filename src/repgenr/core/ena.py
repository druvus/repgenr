"""ENA Portal client: sequencing-run discovery and taxon resolution.

ENA mirrors SRA, needs no key, and answers one query with the run metadata
and the FASTQ locations and checksums, which is why the reads stage reads
from it. Network access only; the records are easy to freeze as a fixture.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from urllib.parse import quote

from . import http
from .contracts import ReadRow
from .errors import UserInputError, WorkdirError

PORTAL_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
COUNT_URL = "https://www.ebi.ac.uk/ena/portal/api/count"
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
    "library_selection",
    "base_count",
    "read_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "first_public",
)

# A search above this many runs is reported before it starts: genera such as
# Salmonella (878821 runs on 2026-10-09) take minutes and gigabytes.
LARGE_SEARCH = 100_000
# Approximate memory of one run record: a fixed part and one per field
# (measured with tracemalloc on Bacillus, 19772 runs, on 2026-10-09: about
# 1500 bytes per run with the 17 RUN_FIELDS and 400 with four fields). The
# normalised rows exist beside the records, so the estimate doubles it. Used
# only to state an expected size in the warning.
_BYTES_PER_RUN = 100
_BYTES_PER_FIELD = 90

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


def is_whole_genome(record: dict) -> bool:
    """Whether a portal record is WGS of genomic source (what _WGS_FILTER selects)."""
    return (
        str(record.get("library_strategy", "")).upper() == "WGS"
        and str(record.get("library_source", "")).upper() == "GENOMIC"
    )


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


def resolve_taxon(name: str, *, division: str | None = None) -> TaxonHit:
    """Resolve a taxon name (synonyms included) to one taxid, or explain why not.

    ``division`` (an ENA taxonomy division such as 'PRO' for prokaryotes)
    chooses among several taxa of the same name, such as the bacterial genus
    Bacillus and the stick-insect genus Bacillus.
    """
    hits: list[dict] = list(http.get_json(f"{TAXONOMY_URL}/any-name/{quote(name)}"))
    if not hits:
        raise UserInputError(f"Nothing in the ENA taxonomy matches {name!r}.")
    if division is not None and len(hits) > 1:
        hits = [h for h in hits if h.get("division") == division] or hits
    if len(hits) > 1:
        listing = "; ".join(
            f"{h.get('scientificName')} ({h.get('rank')}, taxid {h.get('taxId')})" for h in hits
        )
        raise UserInputError(
            f"{name!r} matches {len(hits)} taxa: {listing}. Use the scientific name of one."
        )
    hit = hits[0]
    return TaxonHit(str(hit["taxId"]), str(hit["scientificName"]), str(hit.get("rank", "")))


def count_runs(query: str) -> int:
    """How many read_run records match ``query`` (the portal count endpoint)."""
    answer = http.get_json(
        COUNT_URL, params={"result": "read_run", "query": query, "format": "json"}
    )
    try:
        return int(answer["count"])
    except (KeyError, TypeError, ValueError) as exc:
        raise WorkdirError(f"Unexpected answer from {COUNT_URL}: {answer!r}") from exc


def announce_search(
    query: str, taxon: str, logger: logging.Logger, *, fields: Sequence[str] = RUN_FIELDS
) -> int | None:
    """Count the runs a search will return and log the count before it starts.

    Above :data:`LARGE_SEARCH` runs a warning names the count, the taxon and
    the expected memory; the search proceeds either way. A failed count is
    logged and returns None, since the search itself reports a network fault.
    """
    try:
        n = count_runs(query)
    except WorkdirError as exc:
        logger.warning("Could not count the ENA runs under %s before the search (%s)", taxon, exc)
        return None
    logger.info(
        "ENA counts %d whole-genome sequencing runs under %s; fetching their records", n, taxon
    )
    if n > LARGE_SEARCH:
        size = 2 * n * (_BYTES_PER_RUN + _BYTES_PER_FIELD * len(fields))
        logger.warning(
            "%s has %d ENA runs (more than %d): fetching their records may take "
            "minutes and about %.1f GB of memory. A species or a narrower genus is "
            "a smaller query.",
            taxon,
            n,
            LARGE_SEARCH,
            size / 1e9,
        )
    return n


def search_runs(query: str, *, fields: Sequence[str] = RUN_FIELDS) -> list[dict]:
    """All read_run records matching ``query``, in one request.

    The portal search takes no ``offset`` (it answers 400 'Unsupported param
    offset'), so paging stopped at the first 10000 records; ``limit=0``
    returns every match (19772 runs for tax_tree(1386), Bacillus, in about
    6 s on 2026-10-09). The answer is read as TSV line by line: it carries
    the same string values as the JSON format without repeating the field
    names, and the body is never held whole. Each record is a dict of the
    requested ``fields`` (an empty string where the portal has no value).
    """
    lines = http.iter_lines(
        PORTAL_URL,
        params={
            "result": "read_run",
            "query": query,
            "fields": ",".join(fields),
            "format": "tsv",
            "limit": 0,
        },
    )
    return read_tsv_records(lines, fields=fields, url=PORTAL_URL)


def read_tsv_records(
    lines: Iterable[str], *, fields: Sequence[str] | None = None, url: str = PORTAL_URL
) -> list[dict]:
    """Records from the portal's TSV answer: a header line, then one run per line.

    With ``fields``, a header other than the requested fields (such as an
    HTML page served with status 200) raises :class:`WorkdirError`.

    A line with another number of columns than the header (a body cut short)
    raises :class:`WorkdirError` instead of returning a partial record.
    """
    it = iter(lines)
    header: list[str] | None = None
    for line in it:
        if line:
            header = line.split("\t")
            break
    if header is None:
        return []
    if fields is not None and header != list(fields):
        raise WorkdirError(
            f"Unexpected answer from {url}: the first line is {header[0][:80]!r} where "
            f"the header of the fields {', '.join(fields)} was expected."
        )
    width = len(header)
    records = []
    for line in it:
        if not line:
            continue
        values = line.split("\t")
        if len(values) != width:
            raise WorkdirError(
                f"Malformed record from {url}: {len(values)} columns where the header "
                f"has {width} ({line[:80]!r}); the answer may have been cut short."
            )
        records.append(dict(zip(header, values, strict=True)))
    return records


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
                library_selection=rec.get("library_selection", "") or "",
            )
        )
    return rows
