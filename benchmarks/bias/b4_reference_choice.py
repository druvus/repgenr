"""B4: effect of the reference-genome default on the MSA and tree.

All aligners project onto a single reference, defaulting to the alphabetically
first genome. This experiment builds the phylogeny for the same mixed
50-genome set three times: with the default reference, with another genome
from the same truth cluster, and with a genome from a different cluster. It
then compares the topologies (Robinson-Foulds distance against the
default-reference tree) and whether the truth-cluster bipartition survives.

Writes benchmarks/results/bias_b4.json.

Usage: python -m benchmarks.bias.b4_reference_choice [--aligner sibeliaz]
"""

from __future__ import annotations

import argparse
import json

from benchmarks.metrics import load_truth, newick_splits, robinson_foulds
from benchmarks.run_bench import RESULTS, STORAGE, _repgenr, run_group

SET = "mixed_50_clustered"


def run_one(label: str, reference: str | None, aligner: str, threads: int) -> dict:
    set_dir = STORAGE / "sets" / SET
    work = STORAGE / "work" / f"b4-{label}"
    work.mkdir(parents=True, exist_ok=True)
    out = work / "phylo_out"
    argv = [_repgenr(), "phylo-build", "--genomes-dir", str(set_dir),
            "-o", str(out), "--treebuilder", "fasttree", "--aligner", aligner,
            "--no-outgroup", "-t", str(threads)]
    if reference:
        argv += ["--reference", reference]
    proc = run_group(argv, timeout_s=5400)
    row: dict = {"label": label, "reference": reference, "status": "ok"}
    if proc.returncode != 0:
        return {"label": label, "reference": reference, "status": "failed",
                "stderr_tail": proc.stderr[-1500:]}
    tree = out / "tree" / "tree.nwk"
    if tree.exists():
        row["tree"] = tree.read_text(encoding="utf-8").strip()
    return row


def _truth_split_present(tree_text: str, truth: dict[str, str]) -> bool:
    """Whether the tree contains the bipartition separating the truth clusters."""
    leaves, splits = newick_splits(tree_text)
    sample = next(iter(leaves))
    # tree leaf names may drop the .fasta suffix of the truth keys
    key = (lambda n: n) if sample.endswith(".fasta") else (
        lambda n: n.rsplit(".fasta", 1)[0])
    clusters = sorted({c for c in truth.values()})
    side = frozenset(key(n) for n, c in truth.items() if c == clusters[0]) & leaves
    first = min(leaves)
    canonical = side if first in side else leaves - side
    return canonical in splits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aligner", default="sibeliaz")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    truth = load_truth(STORAGE / "sets" / SET)
    names = sorted(truth)
    default_ref = names[0]  # what the pipeline picks implicitly
    default_cluster = truth[default_ref]
    same_cluster = next(
        n for n in names if n != default_ref and truth[n] == default_cluster)
    cross_cluster = next(n for n in names if truth[n] != default_cluster)

    rows = [
        {"note": "reference variants for the same 50-genome mixed set",
         "aligner": args.aligner,
         "default_ref": default_ref, "default_ref_cluster": default_cluster,
         "samecluster_ref": same_cluster,
         "crosscluster_ref": cross_cluster,
         "crosscluster_ref_cluster": truth[cross_cluster]},
        run_one("default-ref", None, args.aligner, args.threads),
        run_one("samecluster-ref", same_cluster, args.aligner, args.threads),
        run_one("crosscluster-ref", cross_cluster, args.aligner, args.threads),
    ]

    baseline = next((r for r in rows[1:] if r.get("tree")), None)
    for row in rows[1:]:
        tree = row.get("tree")
        if not tree:
            continue
        row["truth_split_present"] = _truth_split_present(tree, truth)
        if baseline and row is not baseline:
            row["vs_default"] = robinson_foulds(baseline["tree"], tree)

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "bias_b4.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print(f"wrote {RESULTS / 'bias_b4.json'}")


if __name__ == "__main__":
    main()
