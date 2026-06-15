"""TreeBuilder interface (tree inference, decoupled from alignment).

MSA-based builders (iqtree, FastTree, RAxML-NG) consume an MSA-FASTA from an
aligner or a SNP typer's core-SNP alignment. Alignment-free builders (mashtree,
sourmash) consume genome files directly and need no MSA source.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[TreeBuilder] = Registry("repgenr.treebuilders")


class InputKind(Enum):
    MSA_FASTA = "msa_fasta"  # needs an MSA source (aligner or snptyper)
    GENOMES = "genomes"      # alignment-free; consumes genome files directly


def as_msa_path(value: Path | Sequence[Path]) -> Path:
    """Narrow a build() argument to a single MSA-FASTA path."""
    if isinstance(value, Path):
        return value
    if isinstance(value, str):
        return Path(value)
    raise TypeError("MSA tree builders expect a single MSA-FASTA path")


def as_genome_list(value: Path | Sequence[Path]) -> list[Path]:
    """Narrow a build() argument to a list of genome paths."""
    if isinstance(value, (str, Path)):
        return [Path(value)]
    return [Path(p) for p in value]


@dataclass
class TreeParams:
    threads: int = 16
    outgroup: str | None = None
    bootstrap: int = 0
    extra: dict = field(default_factory=dict)


class TreeBuilder(ABC):
    capabilities: ToolCapabilities
    input_kind: InputKind = InputKind.MSA_FASTA

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    @abstractmethod
    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        """Build a tree; return the path to a Newick file."""
        raise NotImplementedError
