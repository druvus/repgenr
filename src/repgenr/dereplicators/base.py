"""Dereplicator interface.

An adapter takes a set of genome FASTAs and returns a :class:`DerepResult`
(representatives + cluster membership + per-genome status). The adapter never
writes contract files itself -- the dereplicate stage normalizes the result into
``derep/representatives/`` + ``clusters.tsv`` + ``genome_status.tsv``. This keeps
every dereplicator interchangeable with zero downstream change.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[Dereplicator] = Registry("repgenr.dereplicators")

# status values used in DerepResult.genome_status
STATUS_REPRESENTATIVE = "representative"
STATUS_CONTAINED = "contained"
STATUS_FAIL_QC = "fail_qc"


@dataclass
class DerepParams:
    """Normalized dereplication parameters shared across tools.

    ``extra`` carries tool-specific overrides keyed by adapter name.
    """

    primary_ani: float = 0.90
    secondary_ani: float = 0.99
    aligned_fraction: float = 0.50
    threads: int = 16
    extra: dict = field(default_factory=dict)


@dataclass
class DerepResult:
    """Normalized output every dereplicator must return."""

    representatives: list[Path]
    clusters: dict[str, list[str]]  # representative filename -> contained filenames
    genome_status: dict[str, str]  # genome filename -> status
    genome_information: list[dict] | None = None  # optional checkM-like QC rows


class Dereplicator(ABC):
    """Base class for dereplication adapters."""

    capabilities: ToolCapabilities

    def preflight(self) -> dict[str, str]:
        """Confirm required binaries are present; return resolved versions."""
        return preflight(self.capabilities)

    @abstractmethod
    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        """Cluster ``genomes`` and return representatives + membership."""
        raise NotImplementedError
