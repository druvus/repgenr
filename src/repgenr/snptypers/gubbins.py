"""Gubbins recombination masking of a core-SNP alignment (optional)."""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.binaries import BinarySpec, check_binaries
from ..core.errors import WorkdirError
from ..core.process import run

GUBBINS_BINARY = BinarySpec("run_gubbins.py", version_args=("--version",))


def mask_recombination(core_snp_fasta: Path, out_dir: Path, logger: logging.Logger) -> Path:
    """Run Gubbins; return the recombination-filtered polymorphic-sites FASTA."""
    check_binaries((GUBBINS_BINARY,))
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "gubbins"
    run(
        ["run_gubbins.py", "--prefix", prefix, core_snp_fasta],
        logger=logger,
        cwd=out_dir,
        log_prefix="gubbins",
    )
    filtered = Path(str(prefix) + ".filtered_polymorphic_sites.fasta")
    if not filtered.exists():
        raise WorkdirError("Gubbins did not produce a filtered polymorphic sites FASTA")
    return filtered
