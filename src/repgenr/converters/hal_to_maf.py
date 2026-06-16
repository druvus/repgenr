"""HAL -> MAF wrapper (Cactus output normalization).

Cactus emits a HAL alignment; ``hal2maf`` projects it to MAF against a reference
genome, which :mod:`repgenr.converters.maf_to_fasta` then turns into an
MSA-FASTA.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.containers import run_tool
from ..core.plugins import ToolCapabilities


def hal_to_maf(
    hal_path: str | Path,
    reference: str,
    out_path: str | Path,
    logger: logging.Logger,
    caps: ToolCapabilities | None = None,
) -> Path:
    """Project a HAL to MAF with ``hal2maf``.

    ``caps`` carries the container image when the caller (Cactus) runs
    containerized; ``hal2maf`` ships in the Cactus image. ``--noAncestors`` drops
    ancestral nodes; the Minigraph-Cactus ``_MINIGRAPH_`` backbone pseudo-genome
    is filtered downstream in :func:`maf_to_fasta` (excluding it here via
    ``--targetGenomes`` is brittle: hal2maf errors if any listed genome was
    dropped from the graph).
    """
    hal_path = Path(hal_path)
    out_path = Path(out_path)
    caps = caps or ToolCapabilities(name="hal2maf")
    run_tool(
        caps,
        ["hal2maf", "--refGenome", reference, "--noAncestors", hal_path, out_path],
        logger=logger,
        log_prefix="hal2maf",
    )
    return out_path
