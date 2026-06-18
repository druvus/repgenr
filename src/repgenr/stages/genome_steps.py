"""Discrete genome-download step for data-channel orchestration.

The shared-workdir ``genome`` stage reads selected accessions from the SQLite
manifest and downloads them into the working directory. For the data-channel
pipeline the same download is exposed as a stateless step that reads a portable
``selection.tsv`` (produced by the metadata stage) and writes genome FASTAs into
an output directory -- no manifest, no shared workdir.

Output layout under ``out_dir``::

    genomes/   downloaded selected genomes (named by selection.tsv)
    outgroup/  the outgroup genome, if the selection includes one
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from ..core.contracts import read_selection
from ..core.errors import WorkdirError
from .genome import download_accessions


@dataclass
class GenomeFetchParams:
    selection_tsv: Path
    out_dir: Path
    keep_files: bool = False


def genome_fetch(params: GenomeFetchParams, logger: logging.Logger) -> int:
    """Download the genomes listed in ``selection.tsv`` into ``out_dir``.

    Returns the number of selected (non-outgroup) genomes requested.
    """
    if not params.selection_tsv.exists():
        raise WorkdirError(f"genome-fetch: selection file not found: {params.selection_tsv}")
    rows = read_selection(params.selection_tsv)
    if not rows:
        raise WorkdirError(f"genome-fetch: selection file is empty: {params.selection_tsv}")

    selected = [r for r in rows if not r.is_outgroup]
    outgroup = [r for r in rows if r.is_outgroup]
    scratch = params.out_dir / "scratch"

    if selected:
        filenames = {r.accession: r.filename for r in selected}
        download_accessions(
            list(filenames), filenames, params.out_dir / "genomes", scratch,
            logger, params.keep_files,
        )
    if outgroup:
        og_filenames = {r.accession: r.filename for r in outgroup}
        download_accessions(
            list(og_filenames), og_filenames, params.out_dir / "outgroup", scratch,
            logger, params.keep_files,
        )

    logger.info(
        "genome-fetch: requested %d genomes + %d outgroup into %s",
        len(selected), len(outgroup), params.out_dir,
    )
    return len(selected)
