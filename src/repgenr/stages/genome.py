"""genome stage: download and organize genomes selected by the metadata stage.

Ports ``genome.py``: read selected accessions from the manifest, download them
with NCBI ``datasets`` (dehydrated zip -> rehydrate), and store each as
``genomes/{family}_{genus}_{species}_{accession}.fasta``. The outgroup is
downloaded into ``outgroup/``. Filenames are written back to the manifest.
"""

from __future__ import annotations

import shutil
import zipfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.context import WorkdirContext
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities, preflight

_DATASETS = BinarySpec("datasets", version_args=("--version",))
_DATASETS_CAPS = ToolCapabilities(
    name="datasets",
    required_binaries=(_DATASETS,),
    conda=("conda-forge::ncbi-datasets-cli",),
)


def _run_cmd(cmd, **kwargs):
    """Run the datasets CLI, containerized when a backend is active."""
    return run_tool(_DATASETS_CAPS, cmd, **kwargs)


@dataclass
class GenomeParams:
    accession_list_only: bool = False
    keep_files: bool = False


def run(ctx: WorkdirContext, params: GenomeParams) -> int:
    logger = ctx.logger
    versions = preflight(_DATASETS_CAPS)

    manifest = ctx.manifest
    selected = [g for g in manifest.all_genomes(include_outgroup=False)]
    outgroup = [g for g in manifest.all_genomes(include_outgroup=True) if g.is_outgroup]
    if not selected:
        raise WorkdirError("No selected genomes in manifest. Run the metadata stage first.")

    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    filenames = {g.accession: _output_name(g) for g in selected}

    # Drop genome files that are no longer selected.
    _prune(ctx.genomes_dir, set(filenames.values()), logger)

    to_download = [
        acc for acc, name in filenames.items()
        if not (ctx.genomes_dir / name).exists()
        or (ctx.genomes_dir / name).stat().st_size == 0
    ]
    logger.info(
        "%d to download, %d already present", len(to_download), len(selected) - len(to_download)
    )

    acc_list = ctx.workdir / "ncbi_acc_download_list.txt"
    acc_list.write_text("\n".join(to_download))
    if params.accession_list_only:
        logger.info("Accession list written; stopping (--accession-list-only)")
        return 0

    if to_download:
        _download_batch(ctx, acc_list, filenames, logger, params.keep_files)

    if outgroup:
        _download_outgroup(ctx, outgroup[0], logger)

    # record filenames back into the manifest (one batched transaction)
    for g in selected:
        g.filename = filenames[g.accession]
    manifest.upsert_many(list(selected))

    ctx.config.record_stage(
        "genome",
        params={"downloaded": len(to_download), "total": len(selected)},
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return len(selected)


def _output_name(g) -> str:
    return f"{g.family}_{g.genus}_{g.species}_{g.accession}.fasta"


def _prune(genomes_dir: Path, keep: set[str], logger) -> None:
    for f in genomes_dir.iterdir():
        if f.name not in keep:
            f.unlink()
            logger.info("Removed %s (no longer selected)", f.name)


def _download_batch(ctx, acc_list, filenames, logger, keep_files) -> None:
    zip_path = ctx.workdir / "ncbi_download.zip"
    extract = ctx.workdir / "ncbi_extract"
    if extract.exists():
        shutil.rmtree(extract)

    _run_cmd(
        [
            "datasets", "download", "genome", "accession",
            "--dehydrated", "--inputfile", acc_list, "--filename", zip_path,
        ],
        logger=logger, log_prefix="datasets",
    )
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(extract)
    _run_cmd(
        ["datasets", "rehydrate", "--directory", extract],
        logger=logger, log_prefix="datasets",
    )

    accession_to_name = filenames
    for fna in extract.rglob("*.fna"):
        accession = fna.parent.name
        name = accession_to_name.get(accession)
        if name:
            shutil.move(str(fna), str(ctx.genomes_dir / name))

    if not keep_files:
        zip_path.unlink(missing_ok=True)
        shutil.rmtree(extract, ignore_errors=True)


def _download_outgroup(ctx, outgroup, logger) -> None:
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    zip_path = ctx.workdir / "ncbi_download_outgroup.zip"
    _run_cmd(
        [
            "datasets", "download", "genome", "accession", outgroup.accession,
            "--filename", zip_path,
        ],
        logger=logger, log_prefix="datasets",
    )
    name = _output_name(outgroup)
    with zipfile.ZipFile(zip_path) as zf:
        for item in zf.namelist():
            if item.endswith(".fna") and outgroup.accession in item:
                (ctx.outgroup_dir / name).write_bytes(zf.read(item))
                break
    zip_path.unlink(missing_ok=True)
