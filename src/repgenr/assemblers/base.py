"""Assembler interface: turn one sequencing run's reads into a contig FASTA.

An assembler declares the platforms (``read_types``) and library layouts it
accepts; :func:`select_assembler` picks the first available adapter matching a
run, in a fixed preference order. Selected via ``--assembler <name>`` on the
assemble command, or ``auto``.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight, tool_available

registry: Registry[Assembler] = Registry("repgenr.assemblers")

# Auto-selection order among the adapters that accept a run's platform and layout.
_PREFERENCE = ("skesa", "shovill", "flye")


@dataclass(frozen=True)
class ReadSet:
    """One run's reads as the archive describes them: platform, instrument,
    layout (``PAIRED``/``SINGLE``) and the local FASTQ files, in pair order."""

    run_accession: str
    platform: str  # ENA instrument_platform: ILLUMINA, OXFORD_NANOPORE, PACBIO_SMRT, ...
    instrument_model: str
    layout: str
    files: tuple[Path, ...]
    bases: int = 0


@dataclass
class AssembleParams:
    threads: int = 16
    # A memory hint for tools that cap their own RAM (skesa, shovill).
    memory_gb: int = 16
    # Tool tuning from ``--tool-arg``; each adapter declares the keys it reads
    # in ``capabilities.accepted_extras``.
    extra: dict = field(default_factory=dict)


@dataclass
class AssemblyResult:
    contigs: Path  # the tool's raw contig FASTA (filtered and renamed by the stage)
    tool_stats: dict = field(default_factory=dict)


class Assembler(ABC):
    """Base class for read assemblers."""

    capabilities: ToolCapabilities
    # ENA instrument platforms the tool assembles, and the layouts it accepts.
    read_types: frozenset[str] = frozenset()
    layouts: frozenset[str] = frozenset({"PAIRED", "SINGLE"})

    def preflight(self) -> dict[str, str]:
        """Confirm required binaries are present; return resolved versions."""
        return preflight(self.capabilities)

    def accepts(self, reads: ReadSet) -> bool:
        return reads.platform in self.read_types and reads.layout in self.layouts

    @abstractmethod
    def assemble(
        self,
        reads: ReadSet,
        out_dir: Path,
        params: AssembleParams,
        logger: logging.Logger,
    ) -> AssemblyResult:
        """Assemble ``reads`` under ``out_dir``; return the contig FASTA."""
        raise NotImplementedError


def select_assembler(reg: Registry[Assembler], reads: ReadSet) -> str | None:
    """The preferred available adapter for a run, or None when nothing accepts it."""
    ranked = sorted(
        reg.names(),
        key=lambda n: (_PREFERENCE.index(n) if n in _PREFERENCE else len(_PREFERENCE), n),
    )
    for name in ranked:
        if reg.is_broken(name):
            continue
        cls = reg.get(name)
        probe = cls.__new__(cls)  # accepts() needs no adapter state
        if Assembler.accepts(probe, reads) and tool_available(cls.capabilities):
            return name
    return None
