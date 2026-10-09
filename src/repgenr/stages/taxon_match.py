"""Compare a submitted family, genus and species with a GTDB lineage.

Shared by the post-assembly classifier step of ``assemble`` and the reads
screen (``assemble --screen-reads``), so both apply the same genus rule and
the same ``genus_renamed`` tolerance.
"""

from __future__ import annotations

import re

from ..core.contracts import sanitise_taxon_tokens

AGREE = "agree"
GENUS_RENAMED = "genus_renamed"
DISAGREE = "disagree"

# A GTDB placeholder suffix on a species epithet (coli_A, sanitised to coli-A).
_GTDB_SUFFIX = re.compile(r"-[A-Z]+$")


def same_taxon(submitted: str, gtdb: str) -> bool:
    """Whether two family (or genus) tokens agree, ignoring a GTDB suffix."""
    if not submitted or submitted == "unknown" or not gtdb:
        return False
    return _GTDB_SUFFIX.sub("", submitted) == _GTDB_SUFFIX.sub("", gtdb)


def same_epithet(submitted: str, gtdb: str) -> bool:
    """Whether two species tokens name the same epithet, ignoring a GTDB suffix."""
    if not submitted or submitted == "unknown" or not gtdb:
        return False
    return _GTDB_SUFFIX.sub("", submitted) == _GTDB_SUFFIX.sub("", gtdb)


def gtdb_tokens(lineage: str) -> tuple[str, str, str]:
    """Family, genus and species tokens from a GTDB lineage string."""
    ranks = {}
    for chunk in lineage.split(";"):
        chunk = chunk.strip()
        if len(chunk) > 3 and chunk[1:3] == "__":
            ranks[chunk[0]] = chunk[3:]
    return sanitise_taxon_tokens(ranks.get("f", ""), ranks.get("g", ""), ranks.get("s", ""))


def genus_agreement(submitted: tuple[str, str, str], gtdb: tuple[str, str, str]) -> str:
    """AGREE when the genus tokens are equal; GENUS_RENAMED when the genus
    differs while the family and the species epithet agree (an NCBI name GTDB
    has moved to another genus of the same family); DISAGREE otherwise.

    The family is required because an epithet alone is shared across
    unrelated genera (Klebsiella pneumoniae, Streptococcus pneumoniae).
    """
    if gtdb[1] and gtdb[1] == submitted[1]:
        return AGREE
    if same_taxon(submitted[0], gtdb[0]) and same_epithet(submitted[2], gtdb[2]):
        return GENUS_RENAMED
    return DISAGREE
