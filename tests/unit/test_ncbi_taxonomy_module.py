"""The Entrez taxonomy helper lives in core and the viral package re-exports it."""

from __future__ import annotations


def test_core_module_owns_the_helper_and_viral_re_exports_it() -> None:
    from repgenr.core import ncbi_taxonomy
    from repgenr.viral import entrez

    assert entrez.get_taxon_data_from_entrez is ncbi_taxonomy.get_taxon_data_from_entrez
    assert entrez.TAXNAMES_ORDERED is ncbi_taxonomy.TAXNAMES_ORDERED
