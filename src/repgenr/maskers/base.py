"""Masker interface: post-process a whole-genome alignment (e.g. recombination
removal) before tree building.

A masker takes the SNP typer's whole-genome alignment, not the SNP alignment,
and returns a filtered variable-site FASTA that replaces the typer's core-SNP
alignment. Selected via ``--mask <name>`` on the snptype/phylo commands;
``none`` skips masking without touching the registry.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[Masker] = Registry("repgenr.maskers")


@dataclass
class MaskParams:
    """Tuning passed to every ``Masker.mask`` call.

    ``exclude`` names records (by FASTA id) that the masker must leave out of
    its own analysis but keep in the output: the outgroup. Recombination
    scanners expect closely related isolates and abort on a distant one.
    """

    threads: int = 16
    exclude: frozenset[str] = frozenset()


class Masker(ABC):
    """Base class for whole-genome alignment maskers."""

    capabilities: ToolCapabilities

    def preflight(self) -> dict[str, str]:
        """Confirm required binaries are present; return resolved versions."""
        return preflight(self.capabilities)

    @abstractmethod
    def mask(
        self,
        full_alignment: Path,
        out_dir: Path,
        params: MaskParams,
        logger: logging.Logger,
    ) -> Path:
        """Mask recombinant regions of ``full_alignment``; return the masked
        variable-site FASTA that replaces the typer's core-SNP alignment."""
        raise NotImplementedError
