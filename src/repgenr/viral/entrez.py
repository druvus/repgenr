"""Compatibility re-export: the Entrez taxonomy helper moved to
:mod:`repgenr.core.ncbi_taxonomy` so entry paths other than the viral one can
label genomes by NCBI taxid. Import from there in new code."""

from __future__ import annotations

from ..core import ncbi_taxonomy as _impl
from ..core.ncbi_taxonomy import (  # noqa: F401
    TAXNAMES_ORDERED,
    UNDEFINED_STRAIN,
    _iter_taxa,
    _parse_taxon_element,
    _request_delay,
    _send_query,
    get_taxon_data_from_entrez,
)

__all__ = [name for name in dir(_impl) if not name.startswith("__")]
