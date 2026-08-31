"""B4: effect of the reference-genome default on the MSA and tree.

All aligners project onto a single reference, defaulting to the alphabetically
first genome — on a clonal set, a clone member. This experiment builds the
phylogeny for the same mixed 50-genome set twice: once with the default
reference and once with a diverse-background genome, then compares alignment
size, variant density, and the resulting topologies.

Writes benchmarks/results/bias_b4.json.

Usage: python -m benchmarks.bias.b4_reference_choice [--aligner sibeliaz]
"""

from __future__ import annotations

import argparse
import json

from benchmarks.metrics import load_truth
from benchmarks.run_bench import RESULTS, STORAGE, _repgenr, run_group

SET = "mixed_50_clustered"


def _tree_leafset_and_length(tree_path) -> dict:  # noqa: ANN001
    text = tree_path.read_text(encoding="utf-8").strip()
    leaves = sorted(
        tok.split(":")[0]
        for tok in text.replace("(", ",").replace(")", ",").split(",")
        if tok and not tok.startswith(";") and ":" in tok
    )
    return {"n_leaves": len(leaves), "newick_chars": len(text)}


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
        row.update(_tree_leafset_and_length(tree))
        row["tree"] = tree.read_text(encoding="utf-8").strip()
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aligner", default="sibeliaz")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    truth = load_truth(STORAGE / "sets" / SET)
    names = sorted(truth)
    default_ref = names[0]  # what the pipeline picks implicitly
    # a background (non-clone) genome as the alternative reference
    diverse = next(n for n in names if truth[n] != "clone" and n != default_ref)

    rows = [
        {"note": "default reference (alphabetically first)",
         "default_ref": default_ref, "default_ref_cluster": truth[default_ref],
         "alt_ref": diverse, "alt_ref_cluster": truth[diverse]},
        run_one("default-ref", None, args.aligner, args.threads),
        run_one("diverse-ref", diverse, args.aligner, args.threads),
    ]
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "bias_b4.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print(f"wrote {RESULTS / 'bias_b4.json'}")


if __name__ == "__main__":
    main()
