"""Helpers shared by the viral genome-selection back-ends.

Both selection paths (NCBI Virus records in :mod:`repgenr.viral.selection` and
legacy BV-BRC in :mod:`repgenr.viral.bvbrc`) parse the same ``--target-*``
options and pick an outgroup from the same mashtree distance matrix, so those
two pieces live here. ``VgenomeParams`` is only referenced for typing (the path
modules import it lazily), so there is no import cycle with the stage module.
"""

from __future__ import annotations

from pathlib import Path
from statistics import mean, median, stdev
from typing import TYPE_CHECKING

from ..core.binaries import BinarySpec

if TYPE_CHECKING:
    import logging

    from ..stages.vgenome import VgenomeParams

MASHTREE = BinarySpec("mashtree", version_args=("--version",), min_version="1.2")


def parse_targets(params: VgenomeParams) -> dict[str, list[str]]:
    """Collect the requested taxonomy levels into ``{level: [values]}``."""
    out: dict[str, list[str]] = {}
    if params.target_genus:
        out["genus"] = [x.strip().lower() for x in params.target_genus.split(",")]
    if params.target_species:
        out["species"] = [x.strip().lower() for x in params.target_species.split(",")]
    if params.target_serotype:
        out["serotype"] = [x.strip().lower() for x in params.target_serotype.split(",")]
    if params.target_custom:
        out["custom"] = [x.strip().lower() for x in params.target_custom.split(",")]
    return out


def select_outgroup_from_matrix(matrix: Path, logger: logging.Logger) -> str | None:
    """Pick an outgroup label from a mashtree distance matrix.

    The matrix mixes selected genomes (``S_*``) and outgroup candidates
    (``O_*``). Returns the candidate that is consistently most distant from the
    selection (lowest max distance, then lowest min distance) and clears a
    mean-minus-stdev (floored at the median) distance threshold; falls back to
    the most distant candidate, or None when there is nothing to compare.
    """
    header: list[str] = []
    rows: list[list[str]] = []
    with open(matrix) as fo:
        for enum, line in enumerate(fo):
            parts = line.rstrip("\n").split("\t")
            if enum == 0:
                header = parts
            else:
                rows.append(parts)

    # group statistics for an informational warning
    o_vs_s: dict[str, dict[str, float]] = {}
    s_vs_o_all: list[float] = []
    for idx, seq1 in enumerate(header):
        if idx == 0 or not seq1.startswith("S"):
            continue
        for row in rows:
            seq2 = row[0]
            if not seq2.startswith("O"):
                continue
            val = float(row[idx])
            o_vs_s.setdefault(seq2, {})[seq1] = val
            s_vs_o_all.append(val)
    if not o_vs_s or not s_vs_o_all:
        return None

    threshold = mean(s_vs_o_all) - (stdev(s_vs_o_all) if len(s_vs_o_all) > 1 else 0.0)
    threshold = max(threshold, median(s_vs_o_all))

    # candidates sorted by (lowest max distance, then lowest min distance)
    summary = {
        o: {"min": min(vals.values()), "max": max(vals.values())} for o, vals in o_vs_s.items()
    }
    ordered = sorted(summary.items(), key=lambda x: (-x[1]["max"], x[1]["min"]))
    for candidate, dists in ordered:
        if dists["min"] >= threshold:
            return candidate
    return ordered[0][0] if ordered else None
