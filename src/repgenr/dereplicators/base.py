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
from collections.abc import Collection, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..core.errors import WorkdirError
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
class CompareResult:
    """Output of an all-vs-all comparison (the glance stage's contract).

    ``similarity_csv`` holds pairwise rows with at least the columns
    ``genome1``, ``genome2``, ``similarity``; ``dendrogram`` is an optional
    pre-rendered clustering figure.
    """

    similarity_csv: Path | None = None
    dendrogram: Path | None = None


@dataclass
class DerepResult:
    """Normalized output every dereplicator must return."""

    representatives: list[Path]
    clusters: dict[str, list[str]]  # representative filename -> contained filenames
    genome_status: dict[str, str]  # genome filename -> status
    genome_information: list[dict] | None = None  # optional checkM-like QC rows


_VALID_STATUS = frozenset({STATUS_REPRESENTATIVE, STATUS_CONTAINED, STATUS_FAIL_QC})


def check_result_complete(result: DerepResult, genome_names: Collection[str]) -> None:
    """Refuse a result that silently drops genomes.

    Every input genome must carry a known status; every representative must be
    marked as such; every contained genome must belong to exactly one cluster.
    Adapters that parse tool tables can otherwise return an empty membership
    (e.g. after an upstream column-layout change) and the stage would complete
    with genomes missing from every deliverable.
    """
    names = set(genome_names)
    status = result.genome_status
    missing = sorted(names - status.keys())
    if missing:
        raise WorkdirError(
            f"Dereplication left {len(missing)} of {len(names)} genome(s) without a "
            f"status (e.g. {', '.join(missing[:3])}). The dereplication result is incomplete."
        )
    bad = sorted(f"{g}={s}" for g, s in status.items() if s not in _VALID_STATUS)
    if bad:
        raise WorkdirError(f"Unknown dereplication status value(s): {', '.join(bad[:3])}")

    rep_names = {p.name for p in result.representatives}
    unmarked = sorted(r for r in rep_names if status.get(r) != STATUS_REPRESENTATIVE)
    if unmarked:
        raise WorkdirError(
            f"{len(unmarked)} representative(s) are not marked as such in genome_status "
            f"(e.g. {', '.join(unmarked[:3])})."
        )

    home: dict[str, int] = {}
    for rep, members in result.clusters.items():
        for m in members:
            if m != rep:
                home[m] = home.get(m, 0) + 1
    orphans = sorted(
        g for g, s in status.items() if s == STATUS_CONTAINED and home.get(g, 0) == 0
    )
    if orphans:
        raise WorkdirError(
            f"{len(orphans)} contained genome(s) have no representative "
            f"(e.g. {', '.join(orphans[:3])}). The adapter returned an empty cluster table."
        )
    doubled = sorted(g for g, n in home.items() if n > 1)
    if doubled:
        raise WorkdirError(
            f"{len(doubled)} genome(s) appear in more than one cluster "
            f"(e.g. {', '.join(doubled[:3])})."
        )


class Dereplicator(ABC):
    """Base class for dereplication adapters."""

    capabilities: ToolCapabilities

    def preflight(self) -> dict[str, str]:
        """Confirm required binaries are present; return resolved versions."""
        return preflight(self.capabilities)

    def compare(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        threads: int,
        logger: logging.Logger,
    ) -> CompareResult:
        """Optional capability: all-vs-all comparison for ``repgenr glance``.

        Adapters that can produce a pairwise similarity table override this;
        the default signals the capability is absent.
        """
        raise NotImplementedError(
            f"Dereplicator '{self.capabilities.name}' does not support glance "
            "comparisons (no compare() implementation)."
        )

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
