"""Polisher interface: correct a draft assembly with the reads it was built from.

Long-read assemblies carry systematic indel errors that break genes, which a
completeness check then reads as missing or duplicated markers. A polisher
declares the platforms (``read_types``) it corrects; :func:`select_polisher`
picks the first available adapter accepting a run, in a fixed preference
order, and returns None for reads that need no polishing (Illumina, PacBio
HiFi). Selected via ``--polisher <name>`` on the assemble command, or ``auto``.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path

from ..assemblers.base import ReadSet
from ..core.plugins import Registry, ToolCapabilities, preflight, tool_available

registry: Registry[Polisher] = Registry("repgenr.polishers")

# Auto-selection order among the adapters that accept a run.
_PREFERENCE = ("medaka", "racon")


@dataclass
class PolishParams:
    threads: int = 16
    rounds: int = 1
    # Tool tuning from ``--tool-arg``; each adapter declares the keys it reads
    # in ``capabilities.accepted_extras``.
    extra: dict = field(default_factory=dict)


@dataclass
class PolishResult:
    contigs: Path  # the corrected FASTA (filtered and renamed by the stage)
    rounds: int
    tool_stats: dict = field(default_factory=dict)


class Polisher(ABC):
    """Base class for assembly polishers."""

    capabilities: ToolCapabilities
    # ENA instrument platforms whose assemblies the tool corrects.
    read_types: frozenset[str] = frozenset()

    def preflight(self) -> dict[str, str]:
        """Check the tool's binaries; return resolved versions for provenance."""
        return preflight(self.capabilities)

    def accepts(self, reads: ReadSet) -> bool:
        return reads.platform in self.read_types

    @abstractmethod
    def polish(
        self,
        reads: ReadSet,
        draft: Path,
        out_dir: Path,
        params: PolishParams,
        logger: logging.Logger,
    ) -> PolishResult:
        """Correct ``draft`` with ``reads``; write under ``out_dir``."""


def select_polisher(reg: Registry[Polisher], reads: ReadSet) -> str | None:
    """The preferred available adapter for a run, or None when none applies."""
    ranked = sorted(
        reg.names(),
        key=lambda n: (_PREFERENCE.index(n) if n in _PREFERENCE else len(_PREFERENCE), n),
    )
    for name in ranked:
        if reg.is_broken(name):
            continue
        cls = reg.get(name)
        probe = cls.__new__(cls)  # accepts() needs no adapter state
        if probe.accepts(reads) and tool_available(cls.capabilities):
            return name
    return None


def read_dirs(reads: ReadSet, *paths: Path) -> list[str]:
    """Directories the container backend must mount: the reads and any drafts."""
    return sorted({str(Path(f).resolve().parent) for f in (*reads.files, *paths)})
