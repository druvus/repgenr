"""Gubbins recombination masking of a core-SNP alignment (optional)."""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities, preflight

GUBBINS_BINARY = BinarySpec("run_gubbins.py", version_args=("--version",))
_CAPABILITIES = ToolCapabilities(
    name="gubbins",
    required_binaries=(GUBBINS_BINARY,),
    conda=("bioconda::gubbins",),
)


def mask_recombination(core_snp_fasta: Path, out_dir: Path, logger: logging.Logger) -> Path:
    """Run Gubbins; return the recombination-filtered polymorphic-sites FASTA."""
    preflight(_CAPABILITIES)  # checks host binary, or engine when containerized
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "gubbins"
    run_tool(
        _CAPABILITIES,
        ["run_gubbins.py", "--prefix", prefix, core_snp_fasta],
        logger=logger,
        cwd=out_dir,
        log_prefix="gubbins",
    )
    filtered = Path(str(prefix) + ".filtered_polymorphic_sites.fasta")
    if not filtered.exists():
        raise WorkdirError("Gubbins did not produce a filtered polymorphic sites FASTA")
    return filtered
