"""Classifier interface: assign a GTDB lineage to an assembled genome.

Used by the assemble stage to verify the organism a sequencing run was
submitted under. A classifier takes a reference database the user supplies
(``--gtdb-sketch`` and ``--gtdb-lineages`` for sourmash) and returns one
lineage string per genome in GTDB form (``d__...;p__...;...;s__...``).
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[Classifier] = Registry("repgenr.classifiers")


@dataclass
class ClassifyParams:
    db: Path
    lineages: Path | None = None
    threads: int = 16
    extra: dict = field(default_factory=dict)


@dataclass(frozen=True)
class Classification:
    taxonomy: str  # GTDB lineage string
    rank: str = ""
    score: float | None = None
    db_version: str = ""


def db_version(db: Path | str) -> str:
    """A reference database's identity for provenance: its basename without extensions."""
    name = Path(db).name
    for suffix in (".sig.zip", ".sbt.zip", ".zip", ".sig"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


class Classifier(ABC):
    """Base class for genome classifiers."""

    capabilities: ToolCapabilities

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    @abstractmethod
    def classify(
        self,
        genomes: list[Path],
        out_dir: Path,
        params: ClassifyParams,
        logger: logging.Logger,
    ) -> dict[str, Classification]:
        """Lineage per genome filename; a genome with no match is absent."""
        raise NotImplementedError
