"""Classifier interface: assign a GTDB lineage to an assembled genome.

Used by the assemble stage to verify the organism a sequencing run was
submitted under. A classifier takes a reference database the user supplies
(``--gtdb-sketch`` and ``--gtdb-lineages`` for sourmash) and returns one
lineage string per genome in GTDB form (``d__...;p__...;...;s__...``).
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[Classifier] = Registry("repgenr.classifiers")


@dataclass
class ClassifyParams:
    db: Path
    lineages: Path | None = None
    threads: int = 16
    # A memory budget in GB for the classifier's concurrent work; None leaves
    # concurrency to the thread count alone.
    memory_gb: float | None = None
    extra: dict = field(default_factory=dict)
    # Genome path -> a sourmash signature file of that genome with the
    # parameters of the workdir sketches (k=21, 31, 51; scaled=1000), its
    # signatures named by the genome's record name. Given only to a
    # classifier whose sketch_request() asks for parameters it holds.
    sketches: Mapping[Path, Path] | None = None


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
    # Whether classify() needs ClassifyParams.lineages; checked before any
    # assembly so a missing file is reported at once.
    needs_lineages: bool = False

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    def sketch_request(self, extra: Mapping[str, object]) -> tuple[int, int] | None:
        """The (ksize, scaled) of the query sketches this classifier uses.

        None (the default) means it reads no genome sketches; one that
        returns a pair receives ``ClassifyParams.sketches``.
        """
        return None

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
