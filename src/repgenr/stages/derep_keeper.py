"""Quality-aware representative selection applied after any dereplicator.

Dereplicators pick representatives by connectivity, first-seen order or tool
defaults, which favours the most-sequenced genotype in a cluster. This step
re-picks each cluster's representative by assembly quality when the manifest
carries CheckM-style completeness and contamination, and leaves the adapter's
choice in place otherwise.

The score follows dRep's default weights for the terms that can be computed
without running tools: completeness - 5 x contamination + 0.5 x log10(N50),
with the N50 read from the genome FASTA. Equal scores are broken by higher
completeness, then lower contamination, then higher N50, then the genome
name, so the choice does not depend on the dereplicator. A genome without
quality never becomes the representative; a scored genome replaces an
unscored representative only when it is high quality (MIMAG: completeness
above 90, contamination below 5).
"""

from __future__ import annotations

import gzip
import logging
import math
import zlib
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path

from ..core.errors import UserInputError
from ..core.process import is_gzip
from ..dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE, DerepResult

CONTAMINATION_WEIGHT = 5.0
N50_WEIGHT = 0.5
# MIMAG high-quality draft thresholds (Bowers et al. 2017).
HQ_MIN_COMPLETENESS = 90.0
HQ_MAX_CONTAMINATION = 5.0
# Names listed in a log line before "and N more".
_LOG_NAMES_MAX = 5

Quality = Mapping[str, tuple[float, float]]
# Genome filename -> N50 in bases, or None when the file cannot be found.
N50Of = Callable[[str], "int | None"]


def quality_score(completeness: float, contamination: float) -> float:
    """CheckM part of the score: completeness minus five times contamination."""
    return completeness - CONTAMINATION_WEIGHT * contamination


def keeper_score(completeness: float, contamination: float, n50: int | None) -> float:
    """dRep-style score: ``quality_score`` plus 0.5 x log10(N50).

    An unknown N50 adds nothing.
    """
    score = quality_score(completeness, contamination)
    if n50:
        score += N50_WEIGHT * math.log10(n50)
    return score


def is_high_quality(completeness: float, contamination: float) -> bool:
    """MIMAG high quality: completeness above 90 and contamination below 5."""
    return completeness > HQ_MIN_COMPLETENESS and contamination < HQ_MAX_CONTAMINATION


def genome_n50(path: Path) -> int:
    """N50 of the sequences in a (possibly gzipped) FASTA file, in one pass."""
    lengths: list[int] = []
    current = -1
    try:
        fo = (
            gzip.open(path, "rt", encoding="utf-8", errors="replace")
            if is_gzip(path)
            else open(path, encoding="utf-8", errors="replace")
        )
        with fo:
            for line in fo:
                if line.startswith(">"):
                    if current >= 0:
                        lengths.append(current)
                    current = 0
                elif current >= 0:
                    current += len(line.strip())
    except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise UserInputError(
            f"Genome {path} is a truncated or corrupt gzip file ({exc}); "
            "replace it with a complete copy."
        ) from exc
    if current >= 0:
        lengths.append(current)
    lengths.sort(reverse=True)
    half = sum(lengths) / 2
    running = 0
    for length in lengths:
        running += length
        if running >= half:
            return length
    return 0


class N50Lookup:
    """Genome filename -> N50, read from the first directory holding the file.

    Each genome is read at most once per lookup, and only when asked for, so
    a run reads only the genomes whose scores are compared. ``known`` seeds
    values computed elsewhere (a chunk's ``genome_n50.tsv``). A genome found
    in no directory has no N50 (``None``), which adds nothing to its score.
    """

    def __init__(self, dirs: Iterable[Path], *, known: Mapping[str, int] | None = None) -> None:
        self._dirs = list(dirs)
        self._cache: dict[str, int | None] = dict(known or {})

    def __call__(self, name: str) -> int | None:
        if name not in self._cache:
            self._cache[name] = self._read(name)
        return self._cache[name]

    def _read(self, name: str) -> int | None:
        for d in self._dirs:
            path = d / name
            if path.is_file():
                return genome_n50(path)
        return None

    def computed(self) -> dict[str, int]:
        """The N50 values known so far (genomes that were found)."""
        return {name: n50 for name, n50 in self._cache.items() if n50 is not None}


def _no_n50(_name: str) -> int | None:
    return None


def rank_key(name: str, quality: Quality, n50: N50Of) -> tuple[float, float, float, int, str]:
    """Sort key, best first, for a scored genome: score, completeness,
    contamination, N50, name."""
    completeness, contamination = quality[name]
    value = n50(name)
    return (
        -keeper_score(completeness, contamination, value),
        -completeness,
        contamination,
        -(value or 0),
        name,
    )


def choose_keeper(
    incumbent: str,
    others: Iterable[str],
    quality: Quality,
    n50: N50Of | None = None,
) -> str:
    """The genome that should represent ``incumbent`` and ``others``.

    Only scored genomes compete. When the incumbent (the tool's pick) is
    scored, the best scored genome wins, the incumbent included, ranked by
    :func:`rank_key`. When it is not, a scored genome replaces it only if it
    is high quality (:func:`is_high_quality`); otherwise the incumbent stays.
    """
    n50_of = n50 or _no_n50
    pool = list(dict.fromkeys([incumbent, *others]))
    scored = [n for n in pool if n in quality]
    if incumbent not in quality:
        scored = [n for n in scored if is_high_quality(*quality[n])]
    if not scored:
        return incumbent
    return min(scored, key=lambda n: rank_key(n, quality, n50_of))


def best_score(name: str, quality: Quality, n50: N50Of | None = None) -> float | None:
    """:func:`keeper_score` of a genome, or None when it has no quality."""
    q = quality.get(name)
    if q is None:
        return None
    return keeper_score(q[0], q[1], (n50 or _no_n50)(name))


def log_names(names: list[str]) -> str:
    """At most five names, then "and N more"."""
    shown = ", ".join(names[:_LOG_NAMES_MAX])
    extra = len(names) - _LOG_NAMES_MAX
    return shown if extra <= 0 else f"{shown} and {extra} more"


def log_unpromoted(left: list[str], logger: logging.Logger, what: str = "cluster(s)") -> None:
    """Report groups whose unscored representative was kept although scored
    genomes were available, none of them high quality."""
    if left:
        logger.info(
            "Keeper: %d %s keep the tool's representative, which has no quality, "
            "because no scored genome in them is high quality (completeness > %g, "
            "contamination < %g): %s",
            len(left),
            what,
            HQ_MIN_COMPLETENESS,
            HQ_MAX_CONTAMINATION,
            log_names(sorted(left)),
        )


def rescore_representatives(
    result: DerepResult,
    quality: Quality,
    logger: logging.Logger,
    n50: N50Of | None = None,
) -> tuple[DerepResult, int]:
    """Return a copy of ``result`` whose representatives are chosen by
    :func:`choose_keeper`, and the number of clusters whose representative
    changed. ``n50`` gives each genome's N50; without it the score has no
    N50 term.
    """
    if not quality:
        return result, 0

    rep_paths = {p.name: p for p in result.representatives}
    new_reps: list[Path] = []
    new_clusters: dict[str, list[str]] = {}
    status = dict(result.genome_status)
    swaps = 0
    left_to_tool: list[str] = []

    for rep_name, members in result.clusters.items():
        best_name = choose_keeper(rep_name, members, quality, n50)
        if best_name == rep_name:
            if rep_name not in quality and any(m in quality for m in members):
                left_to_tool.append(rep_name)
            new_reps.append(rep_paths[rep_name])
            new_clusters[rep_name] = list(members)
            continue
        swaps += 1
        new_reps.append(rep_paths[rep_name].with_name(best_name))
        new_clusters[best_name] = sorted(n for n in [rep_name, *members] if n != best_name)
        status[best_name] = STATUS_REPRESENTATIVE
        status[rep_name] = STATUS_CONTAINED
        current = best_score(rep_name, quality, n50)
        logger.info(
            "Keeper: %s replaces %s (score %.2f vs %s)",
            best_name,
            rep_name,
            best_score(best_name, quality, n50),
            "n/a" if current is None else f"{current:.2f}",
        )

    log_unpromoted(left_to_tool, logger)
    if swaps:
        logger.info(
            "Quality-aware keeper changed %d of %d representatives", swaps, len(new_clusters)
        )
    return DerepResult(
        representatives=sorted(new_reps),
        clusters=new_clusters,
        genome_status=status,
        genome_information=result.genome_information,
    ), swaps
