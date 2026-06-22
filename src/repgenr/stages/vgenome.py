"""viral genome stage.

From the viral metadata written by the vmetadata stage, select sequences by
taxonomy and length, write them to ``genomes/``, and determine an outgroup so the
bacterial-style derep/phylo/tree2tax stages work unchanged. Two back-ends share
this entry point and are dispatched on which metadata is present:

* NCBI Virus records (the default) -> :mod:`repgenr.viral.selection`
* legacy BV-BRC tables             -> :mod:`repgenr.viral.bvbrc`

The selection logic lives in those modules; this stage owns the parameters and
the dispatch.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..core.context import WorkdirContext
from ..core.errors import WorkdirError


@dataclass
class VgenomeParams:
    target_genus: str | None = None
    target_species: str | None = None
    target_serotype: str | None = None
    target_custom: str | None = None
    length_all: bool = False
    length_deviation: int = 10
    length_method: str = "median_of_medians"
    length_range: str | None = None
    discard: str | None = None
    no_outgroup: bool = False
    group_segments: bool = False  # ncbi_virus: combine an isolate's segments into one genome
    outgroup_candidates_taxid_min_genomes: int = 5
    glance: bool = False
    print_fasta_headers: bool = False
    ignore_duplicates: bool = False
    keep_files: bool = False
    extra: dict = field(default_factory=dict)


def run(ctx: WorkdirContext, params: VgenomeParams) -> int:
    logger = ctx.logger
    download_wd = ctx.workdir / "virus_download_wd"
    fasta = download_wd / "download.fa"
    records_json = download_wd / "virus_records.json"
    if not fasta.exists():
        raise WorkdirError("Viral metadata missing. Run the vmetadata stage first.")

    # NCBI Virus records present -> default back-end; otherwise the legacy
    # BV-BRC tables must be present.
    if records_json.exists():
        from ..viral.selection import run_records

        return run_records(ctx, params, download_wd, fasta, records_json, logger)

    base_tsv = download_wd / "metadata_base.tsv"
    ncbi_tsv = download_wd / "metadata_ncbi.tsv"
    if not base_tsv.exists() or not ncbi_tsv.exists():
        raise WorkdirError("Viral metadata missing. Run the vmetadata stage first.")

    from ..viral.bvbrc import run_select

    return run_select(ctx, params, fasta, base_tsv, ncbi_tsv, logger)
