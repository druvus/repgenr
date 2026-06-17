"""Micro-benchmark for the dereplication pure-Python hot paths.

Profiles the two O(n^2) Python paths that matter when scaling to 1000s-10000s
genomes: sourmash greedy clustering (dense N x N matrix) and the two-stage
chunk-composition. Run before/after optimizing.

    python scripts/bench_derep_hotpaths.py
"""

from __future__ import annotations

import random
import time
from pathlib import Path

from repgenr.dereplicators.base import DerepResult
from repgenr.dereplicators.sourmash import _greedy_cluster
from repgenr.stages.dereplicate import _compose_two_stage

random.seed(0)


def _synthetic_matrix(n: int, frac_similar: float = 0.001) -> tuple[list[str], list[list[float]]]:
    """N x N similarity matrix; a small fraction of off-diagonal pairs are high."""
    labels = [f"g{i}.fasta" for i in range(n)]
    sim = [[0.0] * n for _ in range(n)]
    for i in range(n):
        sim[i][i] = 1.0
    n_pairs = int(frac_similar * n * n)
    for _ in range(n_pairs):
        i = random.randrange(n)
        j = random.randrange(n)
        if i != j:
            sim[i][j] = sim[j][i] = 0.999
    return labels, sim


def bench_greedy(n: int) -> float:
    labels, sim = _synthetic_matrix(n)
    name_by_label = {label: label for label in labels}
    start = time.perf_counter()
    clusters, status = _greedy_cluster(labels, sim, name_by_label, 0.99)
    elapsed = time.perf_counter() - start
    print(f"  _greedy_cluster      n={n:>5}: {elapsed:7.3f}s  ({len(clusters)} reps)")
    return elapsed


def _synthetic_two_stage(n: int, chunk: int) -> tuple[list[DerepResult], DerepResult]:
    """Stage-1: n genomes in chunks of `chunk`, each chunk keeps ~half as reps.
    Stage-2: folds the chunk reps down to ~n/4 final reps."""
    stage1: list[DerepResult] = []
    s1_reps: list[str] = []
    for c0 in range(0, n, chunk):
        names = [f"g{i}.fasta" for i in range(c0, min(c0 + chunk, n))]
        reps = names[::2]  # every other genome is a rep
        clusters = {r: [m] for r, m in zip(reps, names[1::2], strict=False)}
        status = {r: "representative" for r in reps}
        for members in clusters.values():
            for m in members:
                status[m] = "contained"
        stage1.append(DerepResult([Path(r) for r in reps], clusters, status))
        s1_reps.extend(reps)
    final_reps = s1_reps[::2]
    s2_clusters = {r: [c] for r, c in zip(final_reps, s1_reps[1::2], strict=False)}
    s2_status = {r: "representative" for r in final_reps}
    for members in s2_clusters.values():
        for m in members:
            s2_status[m] = "contained"
    stage2 = DerepResult([Path(r) for r in final_reps], s2_clusters, s2_status)
    return stage1, stage2


def bench_compose(n: int, chunk: int = 500) -> float:
    stage1, stage2 = _synthetic_two_stage(n, chunk)
    start = time.perf_counter()
    result = _compose_two_stage(stage1, stage2)
    elapsed = time.perf_counter() - start
    total = sum(len(v) for v in result.clusters.values()) + len(result.representatives)
    nreps = len(result.representatives)
    print(f"  _compose_two_stage   n={n:>5}: {elapsed:7.3f}s  ({nreps} reps, {total} accounted)")
    return elapsed


if __name__ == "__main__":
    print("sourmash greedy clustering (dense N x N):")
    for n in (1000, 2000, 4000):
        bench_greedy(n)
    print("two-stage chunk composition:")
    for n in (2000, 5000, 10000):
        bench_compose(n)
