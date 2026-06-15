"""HAL -> MAF wrapper (Cactus output normalization).

Cactus emits a HAL alignment; ``hal2maf`` projects it to MAF against a reference
genome, which :mod:`repgenr.converters.maf_to_fasta` then turns into an
MSA-FASTA.
"""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.process import run


def hal_to_maf(
    hal_path: str | Path,
    reference: str,
    out_path: str | Path,
    logger: logging.Logger,
) -> Path:
    hal_path = Path(hal_path)
    out_path = Path(out_path)
    run(
        ["hal2maf", "--refGenome", reference, "--noAncestors", hal_path, out_path],
        logger=logger,
        log_prefix="hal2maf",
    )
    return out_path
