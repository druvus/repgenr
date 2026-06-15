"""SNP typer interface.

A SNP typer takes genomes (and usually a reference) and produces a core-SNP
alignment plus optional VCF and SNP distance matrix. The core-SNP FASTA is both
a standalone typing deliverable and an alternative MSA source for the phylo
stage (interchangeable with an aligner's MSA).
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[SnpTyper] = Registry("repgenr.snptypers")


@dataclass
class SnpParams:
    threads: int = 16
    reference: Path | None = None
    mask: str = "none"  # none | gubbins
    extra: dict = field(default_factory=dict)


@dataclass
class SnpResult:
    core_snp_fasta: Path
    vcf: Path | None = None
    snp_distance_matrix: Path | None = None
    masked: bool = False


class SnpTyper(ABC):
    capabilities: ToolCapabilities
    requires_reference: bool = True

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    @abstractmethod
    def call(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: SnpParams,
        logger: logging.Logger,
    ) -> SnpResult:
        raise NotImplementedError
