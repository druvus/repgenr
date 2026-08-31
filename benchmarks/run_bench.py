"""Resumable benchmark runner for the scaling/bias audit.

Runs the cell matrix from :mod:`benchmarks.cells` against generated sets,
measuring wall time and peak RSS via ``/usr/bin/time -l`` (macOS), and result
metrics (representative count, ARI vs truth, clone-block winner). One JSON per
completed cell under ``benchmarks/results/cells/``; ``--resume`` skips cells
already recorded as ok.

Usage::

    python -m benchmarks.run_bench --tier smoke --resume
    python -m benchmarks.run_bench --only 'derep-*-1000' --resume
    python -m benchmarks.run_bench --collect      # aggregate to summary.tsv
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from benchmarks.cells import Cell, all_cells, tiers
from benchmarks.metrics import (
    adjusted_rand_index,
    clone_representative,
    load_truth,
    partition_from_clusters_tsv,
)

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "benchmarks" / "results"
STORAGE = Path("/Volumes/sekvens2/repgenr")
_RSS_RE = re.compile(r"(\d+)\s+maximum resident set size")


def _python() -> str:
    return sys.executable


def _repgenr() -> str:
    """The repgenr console script installed next to the interpreter."""
    return str(Path(sys.executable).parent / "repgenr")


def _fastas(set_dir: Path) -> list[Path]:
    """Sorted genome FASTAs, skipping macOS AppleDouble ``._*`` files (exFAT)."""
    return sorted(p for p in set_dir.glob("*.fasta") if not p.name.startswith("."))


def _fofn_for(set_dir: Path, work: Path, subset: int | None) -> Path:
    fastas = _fastas(set_dir)
    if subset is not None:
        fastas = fastas[:subset]
    fofn = work / "genomes.fofn"
    fofn.write_text("".join(f"{p}\n" for p in fastas), encoding="utf-8")
    return fofn


def _subset_dir(set_dir: Path, work: Path, subset: int) -> Path:
    """Materialize the first N genomes as copies (sets live on exFAT: no links)."""
    out = work / "genomes_subset"
    out.mkdir(parents=True, exist_ok=True)
    for p in _fastas(set_dir)[:subset]:
        target = out / p.name
        if not target.exists():
            shutil.copy2(p, target)
    return out


def _command(cell: Cell, set_dir: Path, work: Path) -> tuple[list[str], Path]:
    """Build the argv for a cell; returns (argv, clusters_tsv_or_out_dir)."""
    py = _python()
    rg = _repgenr()
    if cell.kind == "derep_step":
        fofn = _fofn_for(set_dir, work, cell.subset)
        out = work / "derep_out"
        argv = [rg, "dereplicate-chunk",
                "--genomes-fofn", str(fofn), "-o", str(out),
                "--tool", cell.tool, "-t", "8", *cell.extra_args]
        return argv, out
    if cell.kind == "derep_dense":
        fofn = _fofn_for(set_dir, work, cell.subset)
        out = work / "derep_out"
        argv = [py, "-m", "benchmarks.dense_driver", str(fofn), str(out), "8"]
        return argv, out
    if cell.kind == "derep_stage":
        wd = work / "wd"
        genomes = wd / "genomes"
        if not genomes.exists():
            genomes.mkdir(parents=True)
            for p in _fastas(set_dir):
                shutil.copy2(p, genomes / p.name)
        argv = [rg, "--force", "dereplicate", "-wd", str(wd),
                "--tool", cell.tool, "-t", "8", *cell.extra_args]
        return argv, wd / "derep"
    if cell.kind == "tree_step":
        genomes_dir = (
            _subset_dir(set_dir, work, cell.subset) if cell.subset is not None else set_dir
        )
        out = work / "phylo_out"
        argv = [rg, "phylo-build",
                "--genomes-dir", str(genomes_dir), "-o", str(out),
                "--treebuilder", cell.tool, "--no-outgroup", "-t", "8",
                *cell.extra_args]
        return argv, out
    raise ValueError(f"unknown cell kind {cell.kind}")


def _measure(argv: list[str], timeout_s: int, log_path: Path) -> dict:
    start = time.monotonic()
    proc = subprocess.run(
        ["/usr/bin/time", "-l", *argv],
        capture_output=True, text=True, timeout=timeout_s,
    )
    wall = time.monotonic() - start
    log_path.write_text(proc.stdout[-20000:] + "\n---stderr---\n" + proc.stderr[-40000:],
                        encoding="utf-8")
    match = _RSS_RE.search(proc.stderr)
    return {
        "wall_s": round(wall, 2),
        "max_rss_mb": round(int(match.group(1)) / 1e6, 1) if match else None,
        "exit_code": proc.returncode,
    }


def _result_metrics(cell: Cell, out_dir: Path, set_dir: Path) -> dict:
    metrics: dict = {}
    clusters_tsv = out_dir / "clusters.tsv"
    if clusters_tsv.exists():
        truth = load_truth(set_dir)
        if cell.subset is not None:
            keep = set(sorted(truth)[: cell.subset])
            truth = {k: v for k, v in truth.items() if k in keep}
        partition = partition_from_clusters_tsv(clusters_tsv)
        metrics["n_representatives"] = len(set(partition.values()))
        metrics["ari_vs_truth"] = round(adjusted_rand_index(truth, partition), 4)
        metrics["clone_representative"] = clone_representative(partition, truth)
    tree = out_dir / "tree" / "tree.nwk"
    if tree.exists():
        content = tree.read_text(encoding="utf-8").strip()
        metrics["tree_ok"] = content.endswith(";")
        metrics["tree_leaves"] = content.count(",") + 1 if content else 0
    return metrics


def run_cell(cell: Cell, *, keep_work: bool = False) -> dict:
    set_dir = STORAGE / "sets" / cell.set_name
    work = STORAGE / "work" / cell.id
    work.mkdir(parents=True, exist_ok=True)
    log_path = STORAGE / "logs" / f"{cell.id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    argv, out_dir = _command(cell, set_dir, work)
    record: dict = {"id": cell.id, "kind": cell.kind, "tool": cell.tool,
                    "set": cell.set_name, "n": cell.n, "argv": argv}
    try:
        record.update(_measure(argv, cell.timeout_s, log_path))
    except subprocess.TimeoutExpired:
        record.update({"exit_code": None, "status": "timeout",
                       "wall_s": cell.timeout_s, "max_rss_mb": None})
    else:
        record["status"] = "ok" if record["exit_code"] == 0 else "failed"
        if record["status"] == "ok":
            record.update(_result_metrics(cell, out_dir, set_dir))
    (RESULTS / "cells").mkdir(parents=True, exist_ok=True)
    (RESULTS / "cells" / f"{cell.id}.json").write_text(
        json.dumps(record, indent=1, sort_keys=True), encoding="utf-8"
    )
    if not keep_work and cell.kind != "derep_stage":  # keep copied workdirs for reuse
        shutil.rmtree(work, ignore_errors=True)
    return record


def is_done(cell: Cell) -> bool:
    path = RESULTS / "cells" / f"{cell.id}.json"
    if not path.exists():
        return False
    return json.loads(path.read_text(encoding="utf-8")).get("status") == "ok"


def collect() -> Path:
    rows = []
    for path in sorted((RESULTS / "cells").glob("*.json")):
        rows.append(json.loads(path.read_text(encoding="utf-8")))
    columns = ["id", "kind", "tool", "set", "n", "status", "wall_s", "max_rss_mb",
               "n_representatives", "ari_vs_truth", "clone_representative",
               "tree_ok", "tree_leaves", "exit_code"]
    out = RESULTS / "summary.tsv"
    with open(out, "w", encoding="utf-8") as fo:
        fo.write("\t".join(columns) + "\n")
        for row in rows:
            fo.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--tier", choices=["smoke", "mid", "heavy"])
    parser.add_argument("--only", help="fnmatch pattern over cell ids")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--collect", action="store_true")
    args = parser.parse_args()

    if args.collect:
        print(f"wrote {collect()}")
        return

    cells = all_cells()
    if args.tier:
        cells = tiers(cells)[args.tier]
    if args.only:
        cells = [c for c in cells if fnmatch.fnmatch(c.id, args.only)]
    if args.resume:
        cells = [c for c in cells if not is_done(c)]
    if args.dry_run:
        for cell in cells:
            print(cell.id)
        print(f"{len(cells)} cell(s) scheduled")
        return
    for i, cell in enumerate(cells, 1):
        print(f"[{i}/{len(cells)}] {cell.id} ...", flush=True)
        record = run_cell(cell)
        print(f"    {record['status']} wall={record.get('wall_s')}s "
              f"rss={record.get('max_rss_mb')}MB "
              f"reps={record.get('n_representatives', '-')}", flush=True)


if __name__ == "__main__":
    main()
