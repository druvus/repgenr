"""NCBI Entrez taxonomy enrichment for viral taxon IDs.

Ports the ``get_taxon_data_from_entrez`` logic from the old ``vmetadata.py``:
batch taxonomy IDs to the Entrez efetch endpoint, parse the returned XML by
splitting on tags (no external XML dependency), and build a per-taxid view of
the taxonomic lineage. Network access is required.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from time import sleep

import requests

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


def _send_query(taxon_ids: Iterable[str]) -> list[str]:
    params = {"db": "taxonomy", "id": list(taxon_ids), "retmode": "xml"}
    response = requests.get(_ENTREZ_BASE + "efetch.fcgi", params=params, timeout=120)
    return response.text.split("\n")


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
        entrez_response: list[str] = []
        for enum, sublist in enumerate(sublists):
            logger.info("Submitting Entrez sublist %d (%d taxids)", enum, len(sublist))
            entrez_response += _send_query(sublist)
            sleep(0.5)  # stay under the Entrez rate limit

        for chunk_raw in _group_taxa(entrez_response):
            parsed = _parse_taxon(chunk_raw, tax_ids_list, tax_ids_to_parse, taxids_alts, logger)
            if parsed is not None:
                taxid, data = parsed
                taxids_data[taxid] = data

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


def _group_taxa(lines: list[str]) -> list[str]:
    groups: list[str] = []
    for line in lines:
        if not line:
            continue
        if line[0] != " " and "<Taxon>" in line:
            groups.append("")
        elif line[0] == " " and groups:
            groups[-1] += line
    return groups


def _parse_taxon(chunk_raw, tax_ids_list, tax_ids_to_parse, taxids_alts, logger):
    chunk = chunk_raw.split("<LineageEx>")[0]
    taxid = chunk.split("<TaxId>")[1].split("</TaxId>")[0]

    taxid_alts: list[str] = []
    if "<AkaTaxIds>" in chunk_raw:
        aka = chunk_raw.split("<AkaTaxIds>")[1].split("</AkaTaxIds>")[0]
        for piece in aka.split("</TaxId>"):
            if "<TaxId>" in piece:
                taxid_alts.append(piece.split("<TaxId>")[1])

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

    scientific_name = chunk.split("<ScientificName>")[1].split("</ScientificName>")[0]
    rank = chunk.split("<Rank>")[1].split("</Rank>")[0]

    level_data: dict[str, dict] = {rank: {"taxid": taxid, "name": scientific_name, "level": rank}}
    if "<LineageEx>" in chunk_raw:
        lineage = chunk_raw.split("<LineageEx>")[1].split("</LineageEx>")[0]
        for piece in lineage.split("</Taxon>"):
            if "<Taxon>" not in piece:
                continue
            sub = piece.split("<Taxon>")[1]
            chunk_taxid = sub.split("<TaxId>")[1].split("</TaxId>")[0]
            if chunk_taxid in taxids_alts:
                for alt in taxids_alts[chunk_taxid]:
                    if alt in set(tax_ids_list):
                        chunk_taxid = alt
                        break
            chunk_name = sub.split("<ScientificName>")[1].split("</ScientificName>")[0]
            chunk_level = sub.split("<Rank>")[1].split("</Rank>")[0]
            level_data[chunk_level] = {
                "taxid": chunk_taxid, "name": chunk_name, "level": chunk_level
            }

    for name in TAXNAMES_ORDERED:
        level_data.setdefault(name, {"taxid": None, "name": None, "level": None})

    if not any(d["taxid"] == taxid for d in level_data.values() if d["taxid"] is not None):
        level_data[UNDEFINED_STRAIN] = {
            "taxid": taxid, "name": scientific_name, "level": UNDEFINED_STRAIN
        }

    return taxid, {"taxid": taxid, "name": scientific_name, "taxdata": level_data}
