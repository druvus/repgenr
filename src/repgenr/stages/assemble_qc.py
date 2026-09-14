"""CheckM2 quality assessment for the assemble stage.

One fixed tool (like the NCBI ``datasets`` client): a single ``checkm2 predict``
over a batch of assemblies returns the completeness and contamination the
quality keeper scores with. The reference database comes from ``--checkm2-db``
or ``CHECKM2DB`` (CheckM2's own variable).
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities, preflight

CHECKM2_CAPS = ToolCapabilities(
    name="checkm2",
    container="quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1",
    conda=("bioconda::checkm2",),
    required_binaries=(BinarySpec("checkm2", version_args=("--version",), min_version="1.0.1"),),
)
CHECKM2_DB_ENV = "CHECKM2DB"
# One predict call per batch: the DIAMOND database is loaded once per call.
_BATCH = 500


def checkm2_db_from_env() -> str | None:
    return os.environ.get(CHECKM2_DB_ENV) or None


def preflight_checkm2() -> dict[str, str]:
    return preflight(CHECKM2_CAPS)


def run_checkm2(
    genomes: list[Path], out_dir: Path, *, db: Path, threads: int, logger: logging.Logger
) -> dict[str, tuple[float, float]]:
    """Completeness and contamination per genome filename, over batches."""
    quality: dict[str, tuple[float, float]] = {}
    for i in range(0, len(genomes), _BATCH):
        batch = genomes[i : i + _BATCH]
        batch_dir = out_dir / f"batch{i // _BATCH}"
        inputs = batch_dir / "input"
        if batch_dir.exists():
            shutil.rmtree(batch_dir, ignore_errors=True)
        inputs.mkdir(parents=True)
        for g in batch:
            os.symlink(g.resolve(), inputs / g.name)
        cmd: list[str | Path] = [
            "checkm2",
            "predict",
            "--input",
            inputs,
            "-x",
            "fasta",
            "--output-directory",
            batch_dir / "out",
            "--threads",
            str(threads),
            "--database_path",
            db,
        ]
        run_tool(
            CHECKM2_CAPS,
            cmd,
            logger=logger,
            log_prefix="checkm2",
            extra_mounts=[
                str(Path(db).resolve().parent),
                *{str(g.resolve().parent) for g in batch},
            ],
        )
        report = batch_dir / "out" / "quality_report.tsv"
        if not report.exists():
            raise WorkdirError(f"CheckM2 wrote no quality_report.tsv under {batch_dir / 'out'}.")
        by_stem = parse_checkm2_report(report)
        for g in batch:
            if g.stem in by_stem:
                quality[g.name] = by_stem[g.stem]
    return quality


def parse_checkm2_report(path: Path) -> dict[str, tuple[float, float]]:
    """``quality_report.tsv`` -> genome stem -> (completeness, contamination)."""
    out: dict[str, tuple[float, float]] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        for rec in csv.DictReader(fo, delimiter="\t"):
            out[rec["Name"]] = (float(rec["Completeness"]), float(rec["Contamination"]))
    return out
