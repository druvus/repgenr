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

from ..core import process
from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.context import WorkdirContext
from ..core.contracts import FASTA_SUFFIXES, genome_filename
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities, preflight

_DATASETS = BinarySpec("datasets", version_args=("--version",))
DATASETS_CAPS = ToolCapabilities(
    name="datasets",
    required_binaries=(_DATASETS,),
    conda=("conda-forge::ncbi-datasets-cli",),
)
_DOWNLOAD_BATCH_SIZE = 5000  # accessions per datasets download/rehydrate call
_EST_BYTES_PER_GENOME = 5_000_000  # rough peak-disk estimate (zip + extract + final)
_MIN_FREE_BYTES = 1_000_000_000  # hard floor: refuse to start a batch under ~1 GB free


def _run_cmd(cmd, **kwargs):
    """Run the datasets CLI, containerized when a backend is active."""
    return run_tool(DATASETS_CAPS, cmd, **kwargs)


@dataclass
class GenomeParams:
    accession_list_only: bool = False
    keep_files: bool = False


def run(ctx: WorkdirContext, params: GenomeParams) -> int:
    logger = ctx.logger
    versions = preflight(DATASETS_CAPS)

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
        download_accessions(
            to_download, filenames, ctx.genomes_dir, ctx.workdir, logger, params.keep_files
        )

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
    return genome_filename(g.family, g.genus, g.species, g.accession)


def _prune(genomes_dir: Path, keep: set[str], logger) -> None:
    # Only ever remove genome FASTA files (never directories or user-placed
    # files): a non-selected FASTA is stale, anything else is left untouched.
    for f in genomes_dir.iterdir():
        if f.name in keep or not f.is_file() or not f.name.endswith(FASTA_SUFFIXES):
            continue
        f.unlink()
        logger.info("Removed %s (no longer selected)", f.name)


def _check_disk(scratch_dir: Path, n_accessions: int, logger) -> None:
    """Refuse to start a download with almost no free disk; warn when tight.

    Genome sizes are unknown ahead of time, so this is a coarse guard against the
    "filled the volume mid-run" failure, not a precise reservation.
    """
    free = shutil.disk_usage(scratch_dir).free
    if free < _MIN_FREE_BYTES:
        raise WorkdirError(
            f"Only {free / 1e9:.1f} GB free under {scratch_dir}; refusing to download "
            f"{n_accessions} genomes. Free disk space or point --outdir at a larger volume."
        )
    estimate = n_accessions * _EST_BYTES_PER_GENOME * 2
    if free < estimate:
        logger.warning(
            "Low disk: ~%.1f GB free, up to ~%.1f GB may be needed for %d genomes.",
            free / 1e9, estimate / 1e9, n_accessions,
        )


def _assert_fasta(path: Path) -> None:
    """Cheap sanity check that a downloaded file is FASTA, not an error page.

    NCBI/GTDB occasionally serve an HTML error body with a 200 status; catching
    it here names the offending file instead of failing cryptically inside an
    aligner far downstream.
    """
    with open(path, "rb") as fo:
        head = fo.read(1)
    if head != b">":
        raise WorkdirError(
            f"Downloaded genome is not FASTA (does not start with '>'): {path.name}. "
            "The download may be an error page; re-run to retry."
        )


def download_accessions(
    accessions: list[str],
    filenames: dict[str, str],
    dest_dir: Path,
    scratch_dir: Path,
    logger,
    keep_files: bool = False,
) -> None:
    """Download genomes in fixed-size sub-batches into ``dest_dir``.

    Each sub-batch is downloaded, rehydrated, and moved into ``dest_dir``
    independently, so a failure loses only one batch and a re-run resumes (the
    caller recomputes the still-missing accessions). This also bounds the size of
    any single dehydrated zip / rehydrate at 1000s-100000s of accessions.
    ``filenames`` maps accession -> output filename; scratch (zips, extracts)
    lives under ``scratch_dir``. ctx-free so both the genome stage and the
    stateless ``genome-fetch`` step reuse it.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    total = len(accessions)
    _check_disk(scratch_dir, total, logger)
    n_batches = (total + _DOWNLOAD_BATCH_SIZE - 1) // _DOWNLOAD_BATCH_SIZE
    for bi, start in enumerate(range(0, total, _DOWNLOAD_BATCH_SIZE)):
        batch = accessions[start : start + _DOWNLOAD_BATCH_SIZE]
        logger.info("Download batch %d/%d (%d accessions)", bi + 1, n_batches, len(batch))
        _download_one_batch(batch, filenames, dest_dir, scratch_dir, logger, keep_files, bi)


def _download_one_batch(
    batch: list[str], filenames: dict[str, str], dest_dir: Path, scratch_dir: Path,
    logger, keep_files: bool, bi: int,
) -> None:
    acc_file = scratch_dir / f"ncbi_acc_batch{bi}.txt"
    acc_file.write_text("\n".join(batch))
    zip_path = scratch_dir / f"ncbi_download_{bi}.zip"
    extract = scratch_dir / f"ncbi_extract_{bi}"
    if extract.exists():
        shutil.rmtree(extract)

    _run_cmd(
        [
            "datasets", "download", "genome", "accession",
            "--dehydrated", "--inputfile", acc_file, "--filename", zip_path,
        ],
        logger=logger, log_prefix="datasets",
    )
    process.unzip(zip_path, extract)
    _run_cmd(
        ["datasets", "rehydrate", "--directory", extract],
        logger=logger, log_prefix="datasets",
    )

    produced: set[str] = set()
    for fna in extract.rglob("*.fna"):
        name = filenames.get(fna.parent.name)
        if name:
            target = dest_dir / name
            shutil.move(str(fna), str(target))
            _assert_fasta(target)
            produced.add(fna.parent.name)

    missing = [acc for acc in batch if acc not in produced]
    if missing:
        logger.warning(
            "NCBI returned no genome for %d of %d accessions in batch %d (e.g. %s)",
            len(missing), len(batch), bi, ", ".join(missing[:3]),
        )

    if not keep_files:
        zip_path.unlink(missing_ok=True)
        acc_file.unlink(missing_ok=True)
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
    try:
        with zipfile.ZipFile(zip_path) as zf:
            for item in zf.namelist():
                if item.endswith(".fna") and outgroup.accession in item:
                    (ctx.outgroup_dir / name).write_bytes(zf.read(item))
                    break
    except zipfile.BadZipFile as exc:
        raise WorkdirError(
            f"Corrupt outgroup download: {zip_path} ({exc}). Re-run to retry."
        ) from exc
    out_path = ctx.outgroup_dir / name
    if out_path.exists():
        _assert_fasta(out_path)
    zip_path.unlink(missing_ok=True)
