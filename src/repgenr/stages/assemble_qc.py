"""CheckM2 quality assessment for the assemble stage.

One fixed tool (like the NCBI ``datasets`` client): a single ``checkm2 predict``
over a batch of assemblies returns the completeness and contamination the
quality keeper scores with. The reference database comes from ``--checkm2-db``
or ``CHECKM2DB`` (CheckM2's own variable).
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import shutil
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import atomic_replace
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
# Per-run store of the scores, beside the run's done marker, so a change of
# only the quality gate reapplies them instead of running CheckM2 again.
CHECKM2_CACHE = "checkm2.json"


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


# --- stored scores ----------------------------------------------------------------


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fo:
        for block in iter(lambda: fo.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def quality_cache_key(contigs: Path, db: str | Path, version: str) -> dict[str, object]:
    """What a stored score depends on: the contigs, the database and the CheckM2 version.

    The database is identified by its resolved path and size, so a database
    replaced in place by another release also invalidates the stored scores.
    """
    db_path = Path(db).expanduser().resolve()
    try:
        db_size = db_path.stat().st_size
    except OSError:
        db_size = -1
    return {
        "contigs_sha256": _sha256(contigs),
        "db": str(db_path),
        "db_size": db_size,
        "checkm2": version,
    }


def load_cached_quality(path: Path, key: dict[str, object]) -> tuple[float, float] | None:
    """The stored (completeness, contamination) when ``path`` holds one for ``key``."""
    try:
        stored = json.loads(path.read_text(encoding="utf-8"))
        if stored.get("key") != key:
            return None
        return float(stored["completeness"]), float(stored["contamination"])
    except (OSError, ValueError, TypeError, KeyError, AttributeError):
        return None


def store_cached_quality(path: Path, key: dict[str, object], quality: tuple[float, float]) -> None:
    payload = {"key": key, "completeness": quality[0], "contamination": quality[1]}
    with atomic_replace(path) as fo:
        fo.write(json.dumps(payload, indent=1))
