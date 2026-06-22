"""NCBI Entrez taxonomy enrichment for viral taxon IDs.

Ports the ``get_taxon_data_from_entrez`` logic from the old ``vmetadata.py``:
batch taxonomy IDs to the Entrez efetch endpoint, parse the returned XML by
splitting on tags (no external XML dependency), and build a per-taxid view of
the taxonomic lineage. Network access is required.
"""

from __future__ import annotations

import logging
import os
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from time import sleep

from ..core import http
from ..core.errors import WorkdirError

# Order in which to present taxonomic names. The last "real" levels are followed
# by the custom "undefined_strain" bucket and the sub-lineage levels.
TAXNAMES_ORDERED: list[str] = [
    "superkingdom", "clade", "kingdom", "phylum", "class", "order",
    "family", "genus", "species", "serotype", "no rank",
]
UNDEFINED_STRAIN = "undefined_strain"
TAXNAMES_ORDERED.append(UNDEFINED_STRAIN)
TAXNAMES_ORDERED += ["subphylum", "subfamily"]

_ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def _send_query(taxon_ids: Iterable[str]) -> str:
    # NCBI asks high-volume callers to identify themselves (tool/email) and
    # raises the rate limit from 3 to 10 req/s when an API key is supplied.
    params: dict[str, object] = {
        "db": "taxonomy", "id": list(taxon_ids), "retmode": "xml", "tool": "repgenr",
    }
    email = os.environ.get("NCBI_EMAIL")
    if email:
        params["email"] = email
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    # Shared retry/backoff session with status checking (a throttle/error page is
    # raised, not silently fed to the XML parser as empty taxonomy).
    return http.get_text(_ENTREZ_BASE + "efetch.fcgi", params=params)


def _request_delay() -> float:
    """Seconds to wait between Entrez requests (3 req/s, or 10 with an API key)."""
    return 0.11 if os.environ.get("NCBI_API_KEY") else 0.34


def get_taxon_data_from_entrez(
    tax_ids_list: Iterable[str],
    logger: logging.Logger,
    num_ids_per_query: int = 100,
    attempts: int = 3,
) -> tuple[dict, set, dict]:
    """Return (taxids_data, taxids_missing_data, taxids_alts)."""
    tax_ids_list = list(tax_ids_list)
    tax_ids_to_parse = list(tax_ids_list)

    taxids_data: dict[str, dict] = {}
    taxids_alts: dict[str, list[str]] = {}

    for loop in range(attempts):
        if not tax_ids_to_parse:
            continue

        sublists = [
            tax_ids_to_parse[i : i + num_ids_per_query]
            for i in range(0, len(tax_ids_to_parse), num_ids_per_query)
        ]
        for enum, sublist in enumerate(sublists):
            logger.info("Submitting Entrez sublist %d (%d taxids)", enum, len(sublist))
            xml_text = _send_query(sublist)
            for taxon_el in _iter_taxa(xml_text):
                parsed = _parse_taxon_element(
                    taxon_el, tax_ids_list, tax_ids_to_parse, taxids_alts, logger
                )
                if parsed is not None:
                    taxid, data = parsed
                    taxids_data[taxid] = data
            sleep(_request_delay())  # stay under the Entrez rate limit

        missing = set(tax_ids_list).difference(taxids_data)
        tax_ids_to_parse = list(missing)
        if missing:
            logger.info("Entrez: %d taxids missing after loop %d, retrying", len(missing), loop)
            sleep(1)

    missing = set(tax_ids_list).difference(taxids_data)
    if missing:
        logger.warning("Entrez: failed to obtain data for %d taxids", len(missing))
        for taxid in missing:
            empty = {
                name: {"taxid": None, "name": None, "level": None} for name in TAXNAMES_ORDERED
            }
            taxids_data[taxid] = {"taxid": taxid, "name": None, "taxdata": empty}
    else:
        logger.info("Entrez data parsed successfully")
    return taxids_data, missing, taxids_alts


def _iter_taxa(xml_text: str) -> list[ET.Element]:
    """Parse an efetch taxonomy response and return its top-level <Taxon> elements.

    An unparseable body (an HTML error/throttle page) or an unexpected root
    (e.g. NCBI's ``<eFetchResult><ERROR>``) raises :class:`WorkdirError` with a
    snippet, rather than silently yielding empty taxonomy.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        snippet = xml_text[:200].replace("\n", " ")
        raise WorkdirError(f"Entrez returned a non-XML response ({exc}): {snippet!r}") from exc
    if root.tag != "TaxaSet":
        snippet = xml_text[:200].replace("\n", " ")
        raise WorkdirError(f"Unexpected Entrez response root <{root.tag}>: {snippet!r}")
    return root.findall("Taxon")


def _parse_taxon_element(el, tax_ids_list, tax_ids_to_parse, taxids_alts, logger):
    taxid = el.findtext("TaxId") or ""

    taxid_alts: list[str] = []
    aka = el.find("AkaTaxIds")
    if aka is not None:
        taxid_alts = [t.text for t in aka.findall("TaxId") if t.text]

    aliases = [taxid, *taxid_alts]
    for alias in aliases:
        taxids_alts[alias] = aliases

    if taxid not in tax_ids_to_parse:
        for alt in taxid_alts:
            if alt in tax_ids_to_parse:
                taxid = alt
                break
    if taxid not in tax_ids_list:
        logger.info("Entrez returned an unexpected taxid %s; ignoring", taxid)
        return None

    scientific_name = el.findtext("ScientificName") or ""
    rank = el.findtext("Rank") or ""

    level_data: dict[str, dict] = {rank: {"taxid": taxid, "name": scientific_name, "level": rank}}
    lineage = el.find("LineageEx")
    if lineage is not None:
        for sub in lineage.findall("Taxon"):
            chunk_taxid = sub.findtext("TaxId") or ""
            if chunk_taxid in taxids_alts:
                for alt in taxids_alts[chunk_taxid]:
                    if alt in set(tax_ids_list):
                        chunk_taxid = alt
                        break
            chunk_level = sub.findtext("Rank") or ""
            level_data[chunk_level] = {
                "taxid": chunk_taxid, "name": sub.findtext("ScientificName") or "",
                "level": chunk_level,
            }

    for name in TAXNAMES_ORDERED:
        level_data.setdefault(name, {"taxid": None, "name": None, "level": None})

    if not any(d["taxid"] == taxid for d in level_data.values() if d["taxid"] is not None):
        level_data[UNDEFINED_STRAIN] = {
            "taxid": taxid, "name": scientific_name, "level": UNDEFINED_STRAIN
        }

    return taxid, {"taxid": taxid, "name": scientific_name, "taxdata": level_data}
