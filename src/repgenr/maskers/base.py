"""Masker interface: post-process a core-SNP alignment (e.g. recombination
removal) before tree building.

A masker takes the SNP typer's core-SNP FASTA and returns a filtered FASTA.
Selected via ``--mask <name>`` on the snptype/phylo commands; ``none`` skips
masking without touching the registry.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[Masker] = Registry("repgenr.maskers")


class Masker(ABC):
    """Base class for core-SNP alignment maskers."""

    capabilities: ToolCapabilities

    def preflight(self) -> dict[str, str]:
        """Confirm required binaries are present; return resolved versions."""
        return preflight(self.capabilities)

    @abstractmethod
    def mask(
        self,
        core_snp_fasta: Path,
        out_dir: Path,
        logger: logging.Logger,
    ) -> Path:
        """Filter ``core_snp_fasta``; return the path to the filtered FASTA."""
        raise NotImplementedError
