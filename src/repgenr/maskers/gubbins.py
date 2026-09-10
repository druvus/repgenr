"""Gubbins recombination masking of a whole-genome alignment."""

from __future__ import annotations

import logging
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from .base import Masker, MaskParams

_ALLOWED = set("ACGTNacgtn-")


def sanitise_alignment(source: Path, dest: Path, logger: logging.Logger) -> Path:
    """Copy ``source`` with every character outside ACGTN- replaced by N.

    Gubbins accepts only those symbols; SNP typers can leave IUPAC ambiguity
    codes (e.g. a heterozygous call rendered by ``bcftools consensus``) in a
    whole-genome alignment. Returns ``dest`` (or ``source`` when unchanged).
    """
    replaced = 0
    with open(source, encoding="utf-8") as fi, open(dest, "w", encoding="utf-8") as fo:
        for line in fi:
            if line.startswith(">"):
                fo.write(line)
                continue
            body = line.rstrip("\n")
            fixed = "".join(c if c in _ALLOWED else "N" for c in body)
            replaced += sum(1 for a, b in zip(body, fixed, strict=True) if a != b)
            fo.write(fixed + "\n")
    if replaced:
        logger.warning(
            "Replaced %d ambiguous base(s) with N before Gubbins (only ACGTN- are accepted).",
            replaced,
        )
        return dest
    dest.unlink(missing_ok=True)
    return source


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
        cleaned = sanitise_alignment(full_alignment, out_dir / "input_alignment.fasta", logger)
        argv: list[str | Path] = [
            "run_gubbins.py",
            "--threads",
            str(params.threads),
            "--prefix",
            prefix,
            cleaned,
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
