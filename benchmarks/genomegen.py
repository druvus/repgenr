"""Synthetic genome sets with controlled similarity structure (audit tooling).

Generates seeded FASTA sets for the scaling/bias benchmarks: a random ancestor
sequence is mutated into cluster founders and cluster members, giving direct
control over intra- and inter-cluster identity. The in-scope tools (skani,
sourmash, mash, minimap2) measure k-mer/alignment identity, which this model
controls exactly; gene-model tools (parsnp) would not behave realistically on
random sequence and are outside the native benchmark scope.

Scenarios:

* ``balanced(n)``  -- n/50 clusters at ~96% inter- and ~99.6% intra-cluster ANI.
* ``clonal(n, f)`` -- one near-identical clone block (>=99.99% ANI) holding a
  fraction ``f`` of the set, plus a balanced-style diverse background.
* ``mixed(n)``     -- clusters over a divergence gradient.

``order`` controls how accessions (and therefore alphabetical sort position,
itself a bias mechanism under audit) map onto clusters: ``clustered`` gives the
clone block the lowest accessions (contiguous, sorts first), ``interleaved``
round-robins, ``random`` shuffles.

Usage::

    python -m benchmarks.genomegen --scenario clonal --n 1000 \
        --out /Volumes/sekvens2/repgenr/sets/clonal_1000_clustered --seed 1
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

_BASES = "ACGT"
_LINE = 70

# Per-genome divergence from its cluster founder (pairwise intra-cluster
# divergence is ~2x these). Chosen so intra-cluster ANI is ~99.6% (balanced)
# and >=99.99% (clone block), with ~96% between balanced cluster founders.
_FOUNDER_DIVERGENCE = 0.02
_MEMBER_DIVERGENCE = 0.002
_CLONE_DIVERGENCE = 0.00005


def pairwise_identity(a: str, b: str) -> float:
    """Fraction of positions identical between two equal-length sequences."""
    length = min(len(a), len(b))
    matches = sum(1 for x, y in zip(a[:length], b[:length], strict=False) if x == y)
    return matches / length


def _random_sequence(rng: random.Random, length: int) -> bytearray:
    return bytearray(rng.choices(_BASES.encode(), k=length))


def _mutate(rng: random.Random, seq: bytearray, divergence: float) -> bytearray:
    """Substitute ``round(divergence * len)`` distinct positions."""
    out = bytearray(seq)
    k = round(divergence * len(seq))
    for pos in rng.sample(range(len(seq)), k):
        current = chr(out[pos])
        out[pos] = ord(rng.choice([b for b in _BASES if b != current]))
    return out


def _cluster_plan(
    scenario: str, n: int, clone_fraction: float
) -> list[tuple[str, int, float, float]]:
    """(cluster_id, size, founder_divergence, member_divergence) per cluster."""
    if scenario == "clonal":
        clone_size = max(2, round(clone_fraction * n))
        rest = n - clone_size
        plan = [("clone", clone_size, _FOUNDER_DIVERGENCE, _CLONE_DIVERGENCE)]
        plan += _balanced_plan(rest, start=1)
        return plan
    if scenario == "balanced":
        return _balanced_plan(n, start=1)
    if scenario == "mixed":
        n_clusters = max(2, n // 50)
        sizes = _split_sizes(n, n_clusters)
        return [
            (f"m{i + 1}", size, 0.005 + 0.045 * i / max(1, n_clusters - 1),
             _MEMBER_DIVERGENCE)
            for i, size in enumerate(sizes)
        ]
    raise ValueError(f"Unknown scenario '{scenario}' (balanced|clonal|mixed).")


def _balanced_plan(n: int, start: int) -> list[tuple[str, int, float, float]]:
    if n <= 0:
        return []
    n_clusters = max(2, n // 50) if n >= 4 else 2
    n_clusters = min(n_clusters, n)
    sizes = _split_sizes(n, n_clusters)
    return [
        (f"sp{start + i}", size, _FOUNDER_DIVERGENCE, _MEMBER_DIVERGENCE)
        for i, size in enumerate(sizes)
    ]


def _split_sizes(n: int, parts: int) -> list[int]:
    base, extra = divmod(n, parts)
    return [base + (1 if i < extra else 0) for i in range(parts)]


def _assign_accessions(
    rng: random.Random, memberships: list[str], order: str
) -> list[int]:
    """Accession index per genome (position i in ``memberships``).

    clustered: input order (clone block first) -> contiguous low accessions.
    interleaved: round-robin across clusters. random: seeded shuffle.
    """
    indices = list(range(len(memberships)))
    if order == "clustered":
        return indices
    if order == "random":
        rng.shuffle(indices)
        return indices
    if order == "interleaved":
        by_cluster: dict[str, list[int]] = {}
        for i, cluster in enumerate(memberships):
            by_cluster.setdefault(cluster, []).append(i)
        rotated: list[int] = []
        queues = list(by_cluster.values())
        while any(queues):
            for queue in queues:
                if queue:
                    rotated.append(queue.pop(0))
        # rotated[k] = genome receiving accession k
        return rotated
    raise ValueError(f"Unknown order '{order}' (clustered|interleaved|random).")


def _write_fasta(path: Path, header: str, seq: bytearray) -> None:
    with open(path, "w", encoding="utf-8") as fo:
        fo.write(f">{header}\n")
        data = seq.decode()
        for pos in range(0, len(data), _LINE):
            fo.write(data[pos : pos + _LINE] + "\n")


def generate_set(
    out_dir: Path,
    *,
    scenario: str,
    n: int,
    seed: int,
    genome_length: int = 2_000_000,
    clone_fraction: float = 0.4,
    order: str = "clustered",
) -> dict:
    """Generate a genome set; return (and write) the truth record."""
    rng = random.Random(seed)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ancestor = _random_sequence(rng, genome_length)
    plan = _cluster_plan(scenario, n, clone_fraction)

    memberships: list[str] = []
    divergences: list[float] = []
    sequences: list[bytearray] = []
    for cluster_id, size, founder_div, member_div in plan:
        founder = _mutate(rng, ancestor, founder_div)
        for _ in range(size):
            memberships.append(cluster_id)
            divergences.append(founder_div + member_div)
            sequences.append(_mutate(rng, founder, member_div))

    # accession_of[i] is the accession index for genome i. The species token
    # repeats the index so alphabetical sort order == accession order, making
    # the `order` knob fully control where clone members land in the sorted
    # listing (the mechanism under audit). Truth clusters live in truth.json,
    # not in the filename taxonomy.
    accession_of = _assign_accessions(rng, memberships, order)
    clusters: dict[str, str] = {}
    divergence: dict[str, float] = {}
    for i, cluster in enumerate(memberships):
        idx = accession_of[i]
        name = f"Benchfam_Benchgen_g{idx:06d}_GCF9{idx:06d}.1.fasta"
        clusters[name] = cluster
        divergence[name] = divergences[i]
        _write_fasta(out_dir / name, f"GCF9{idx:06d}.1", sequences[i])

    truth = {
        "scenario": scenario,
        "n": n,
        "seed": seed,
        "order": order,
        "genome_length": genome_length,
        "clone_cluster": "clone" if scenario == "clonal" else None,
        "clusters": clusters,
        "divergence": divergence,
    }
    (out_dir / "truth.json").write_text(
        json.dumps(truth, indent=1, sort_keys=True), encoding="utf-8"
    )
    return truth


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--scenario", required=True, choices=["balanced", "clonal", "mixed"])
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--genome-length", type=int, default=2_000_000)
    parser.add_argument("--clone-fraction", type=float, default=0.4)
    parser.add_argument("--order", default="clustered",
                        choices=["clustered", "interleaved", "random"])
    args = parser.parse_args()
    truth = generate_set(
        args.out, scenario=args.scenario, n=args.n, seed=args.seed,
        genome_length=args.genome_length, clone_fraction=args.clone_fraction,
        order=args.order,
    )
    print(f"Wrote {len(truth['clusters'])} genomes to {args.out}")


if __name__ == "__main__":
    main()
