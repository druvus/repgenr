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


_PAIR_SUFFIXES = (("_1", "_2"), ("_R1", "_R2"))


def split_reads(reads: ReadSet) -> tuple[tuple[Path, Path] | None, list[Path]]:
    """The read pair of a run, and the files that are not part of it.

    ENA lists up to three FASTQ files for a paired run: ``<run>_1`` and
    ``<run>_2`` plus an orphan ``<run>`` file of reads whose mate was dropped.
    The pair is found by the ``_1``/``_2`` (or ``_R1``/``_R2``) stem suffix;
    two files with no such suffix are taken as the pair as given. Everything
    else is returned as single-end files, in the order listed.
    """
    files = [Path(f) for f in reads.files]

    def stem(path: Path) -> str:
        name = path.name
        for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"):
            if name.endswith(ext):
                return name[: -len(ext)]
        return path.stem

    for one, two in _PAIR_SUFFIXES:
        firsts = [f for f in files if stem(f).endswith(one)]
        seconds = [f for f in files if stem(f).endswith(two)]
        for r1 in firsts:
            base = stem(r1)[: -len(one)]
            r2 = next((f for f in seconds if stem(f) == base + two), None)
            if r2 is not None:
                singles = [f for f in files if f not in (r1, r2)]
                return (r1, r2), singles
    if len(files) == 2:
        return (files[0], files[1]), []
    return None, files


def _read_dirs(reads: ReadSet) -> list[str]:
    """The directories holding a run's reads, for the container backend's mounts."""
    return sorted({str(Path(f).resolve().parent) for f in reads.files})


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
