"""GFA -> MSA-FASTA (graph pangenome, optional).

Turning a pangenome graph (Cactus/PGGB GFA) into a column-wise MSA requires a
graph-aware tool such as ``odgi``. When ``odgi`` is available we build a graph
and emit a per-path MSA; otherwise we raise a clear error pointing at the
HAL -> MAF -> FASTA route, which is the supported default for Cactus.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from ..core.errors import MissingBinaryError
from ..core.process import run


def gfa_to_fasta(
    gfa_path: str | Path,
    out_path: str | Path,
    logger: logging.Logger,
) -> Path:
    gfa_path = Path(gfa_path)
    out_path = Path(out_path)
    if shutil.which("odgi") is None:
        raise MissingBinaryError(
            "GFA -> MSA conversion needs 'odgi'. For Cactus, prefer the "
            "HAL -> MAF -> FASTA route (default) instead."
        )
    og = out_path.with_suffix(".og")
    run(["odgi", "build", "-g", gfa_path, "-o", og], logger=logger, log_prefix="odgi")
    # odgi paths -f emits path sequences; this is the graph's MSA view.
    run(["odgi", "paths", "-i", og, "-f", out_path], logger=logger, log_prefix="odgi")
    return out_path
