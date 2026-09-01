"""Gubbins recombination masking of a whole-genome alignment."""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from .base import Masker, MaskParams


class GubbinsMasker(Masker):
    capabilities = ToolCapabilities(
        name="gubbins",
        required_binaries=(BinarySpec("run_gubbins.py", version_args=("--version",)),),
        conda=("bioconda::gubbins",),
    )

    def mask(
        self,
        full_alignment: Path,
        out_dir: Path,
        params: MaskParams,
        logger: logging.Logger,
    ) -> Path:
        """Run Gubbins on the whole-genome alignment; return the
        recombination-filtered polymorphic-sites FASTA."""
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = out_dir / "gubbins"
        argv: list[str | Path] = [
            "run_gubbins.py", "--threads", str(params.threads),
            "--prefix", prefix, full_alignment,
        ]
        run_tool(
            self.capabilities,
            argv,
            logger=logger,
            cwd=out_dir,
            log_prefix="gubbins",
        )
        filtered = Path(str(prefix) + ".filtered_polymorphic_sites.fasta")
        if not filtered.exists():
            raise WorkdirError("Gubbins did not produce a filtered polymorphic sites FASTA")
        return filtered
