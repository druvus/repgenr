"""B1: which genome represents an over-represented clone block, and why.

Runs each dereplicator on the same 50-genome clonal set under three accession
orderings (clustered, interleaved, random). The genomes are identical across
orderings up to renaming, so a selection driven by content would pick the same
underlying genome; a selection driven by filename order tracks the ordering.

Writes benchmarks/results/bias_b1.json.

Usage: python -m benchmarks.bias.b1_clone_rep [--tools skder,galah,sourmash]
"""

from __future__ import annotations

import argparse
import json
import subprocess

from benchmarks.metrics import clone_representative, load_truth, partition_from_clusters_tsv
from benchmarks.run_bench import RESULTS, STORAGE, _env, _fofn_for, _repgenr

ORDERS = ("clustered", "interleaved", "random")


def run_one(tool: str, order: str) -> dict:
    set_dir = STORAGE / "sets" / f"clonal_50_{order}"
    work = STORAGE / "work" / f"b1-{tool}-{order}"
    work.mkdir(parents=True, exist_ok=True)
    out = work / "derep_out"
    fofn = _fofn_for(set_dir, work, None)
    argv = [
        _repgenr(),
        "dereplicate-chunk",
        "--genomes-fofn",
        str(fofn),
        "-o",
        str(out),
        "--tool",
        tool,
        "-t",
        "8",
    ]
    proc = subprocess.run(argv, capture_output=True, text=True, timeout=1800, env=_env())
    if proc.returncode != 0:
        return {
            "tool": tool,
            "order": order,
            "status": "failed",
            "stderr_tail": proc.stderr[-1500:],
        }
    truth = load_truth(set_dir)
    partition = partition_from_clusters_tsv(out / "clusters.tsv")
    winner = clone_representative(partition, truth)
    clone_members = sorted(g for g, c in truth.items() if c == "clone")
    return {
        "tool": tool,
        "order": order,
        "status": "ok",
        "n_representatives": len(set(partition.values())),
        "clone_representative": winner,
        "clone_rep_is_alphabetically_first_clone": winner == clone_members[0],
        "first_clone_member": clone_members[0],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tools", default="skder,galah,sourmash")
    args = parser.parse_args()

    rows = []
    for tool in args.tools.split(","):
        for order in ORDERS:
            row = run_one(tool.strip(), order)
            rows.append(row)
            print(json.dumps(row))
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "bias_b1.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print(f"wrote {RESULTS / 'bias_b1.json'}")


if __name__ == "__main__":
    main()
