"""B5: what happens on a (near-)pure clone set with almost no variation.

Generates a small all-clone set (clone divergence ~0.005 percent) and runs the
aligner -> fasttree path. The audit expects a degenerate star-like tree with no
warning anywhere except the `simple` snptyper's zero-SNP guard.

Writes benchmarks/results/bias_b5.json.

Usage: python -m benchmarks.bias.b5_zero_snp [--aligner sibeliaz]
"""

from __future__ import annotations

import argparse
import json
import sys

from benchmarks.genomegen import generate_set
from benchmarks.run_bench import RESULTS, STORAGE, _repgenr, run_group


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aligner", default="sibeliaz")
    args = parser.parse_args()

    set_dir = STORAGE / "sets" / "pureclone_20"
    if not (set_dir / "truth.json").exists():
        generate_set(set_dir, scenario="clonal", n=20, seed=1, clone_fraction=1.0)
        print(f"generated {set_dir}", file=sys.stderr)

    work = STORAGE / "work" / "b5-pureclone"
    work.mkdir(parents=True, exist_ok=True)
    out = work / "phylo_out"
    argv = [_repgenr(), "phylo-build", "--genomes-dir", str(set_dir),
            "-o", str(out), "--treebuilder", "fasttree", "--aligner", args.aligner,
            "--no-outgroup", "-t", "8"]
    proc = run_group(argv, timeout_s=5400)

    row: dict = {"aligner": args.aligner, "exit_code": proc.returncode}
    stderr = proc.stderr[-4000:]
    row["warned_low_diversity"] = "diverg" in stderr.lower() or "warning" in stderr.lower()
    tree = out / "tree" / "tree.nwk"
    if tree.exists():
        text = tree.read_text(encoding="utf-8").strip()
        row["tree"] = text
        # branch lengths: a degenerate tree is all ~zero
        lengths = [
            float(tok.split(":")[1].rstrip(");"))
            for tok in text.replace("(", ",").replace(")", ",").split(",")
            if ":" in tok
        ]
        row["max_branch_length"] = max(lengths) if lengths else None
    row["stderr_tail"] = stderr[-2000:]
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "bias_b5.json").write_text(json.dumps(row, indent=1), encoding="utf-8")
    print(f"wrote {RESULTS / 'bias_b5.json'}")


if __name__ == "__main__":
    main()
