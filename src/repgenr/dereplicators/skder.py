"""skDER dereplication adapter (skani-based).

skDER scales to tens of thousands of genomes natively (no chunking needed). It
estimates pairwise ANI/AF with skani and picks representatives in one round.

CLI used::

    skder -g <files> -o <out_dir> -i <ANI%> -f <AF%> -c <threads> \
          -d dynamic|greedy|low_mem_greedy

Two robustness notes:

* skDER is run in a local scratch directory and the results are copied into the
  working directory afterwards. skDER's ``determineN50`` does
  ``os.listdir(N50/)`` and parses every file as a result row, which breaks on
  exFAT/NTFS volumes where macOS leaves ``._*`` AppleDouble companions. Keeping
  skDER's output on the local (APFS/ext4) temp filesystem avoids that.
* Plain skDER reports representatives only (membership is a CiDDER feature), so
  non-representatives are assigned to their closest representative using skDER's
  ``Skani_Triangle_Edge_Output.txt`` pairwise ANI/AF table.
"""

from __future__ import annotations

import logging
import shutil
import tempfile
from collections.abc import Iterator, Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import ToolExecutionError, WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import link_or_copy
from .base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
)

_FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fas")


class SkderDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="skder",
        required_binaries=(
            BinarySpec("skder", version_args=("-v",), min_version="1.0"),
            BinarySpec("skani", version_args=("-V",), min_version="0.2"),
        ),
        default_params={"mode": "greedy"},
        recommended_max_genomes=None,  # scales natively
        supports_native_scaling=True,
        threads_param="-c",
        conda=("bioconda::skder",),
    )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        out_dir.mkdir(parents=True, exist_ok=True)

        # Run skDER on a local temp filesystem to avoid exFAT/NTFS ._* breakage,
        # then copy its results back next to the working directory.
        local_tmp = Path(tempfile.mkdtemp(prefix="repgenr_skder_"))
        result_dir = local_tmp / "skder_out"
        mode = params.extra.get("mode", self.capabilities.default_params["mode"])
        ani_pct = _as_percent(params.secondary_ani)
        af_pct = _as_percent(params.aligned_fraction)
        cmd = [
            "skder",
            "-g", *[str(g) for g in genomes],
            "-o", result_dir,
            "-i", ani_pct,
            "-f", af_pct,
            "-c", params.threads,
            "-d", mode,
        ]
        try:
            run_tool(self.capabilities, cmd, logger=logger, log_prefix="skder")
            staged = out_dir / "skder_out"
            if staged.exists():
                shutil.rmtree(staged)
            # Hardlink the result files when possible (no extra disk; the temp dir
            # is removed below, but hardlinks keep the data alive). Falls back to a
            # copy across filesystems (temp vs workdir) or on exFAT.
            shutil.copytree(result_dir, staged, copy_function=link_or_copy)
        finally:
            shutil.rmtree(local_tmp, ignore_errors=True)

        return _parse_skder_output(staged, genomes, float(ani_pct), float(af_pct), logger)


def _as_percent(value: float) -> str:
    """skDER expects identity/fraction as a percentage (e.g. 99.0, 50.0)."""
    pct = value * 100 if value <= 1.0 else value
    return f"{pct:g}"


def _parse_skder_output(
    result_dir: Path,
    genomes: Sequence[Path],
    ani_cutoff: float,
    af_cutoff: float,
    logger: logging.Logger,
) -> DerepResult:
    rep_dir = _find_representatives_dir(result_dir)
    if rep_dir is None:
        raise WorkdirError(
            f"Could not locate skDER representative genomes under {result_dir}. "
            "Confirm skDER ran successfully."
        )
    representatives = sorted(
        p for p in rep_dir.iterdir()
        if not p.name.startswith(".") and p.suffix in _FASTA_SUFFIXES
    )
    if not representatives:
        raise ToolExecutionError(["skder"], 1, "skDER produced no representative genomes")

    rep_names = {p.name for p in representatives}
    clusters: dict[str, list[str]] = {p.name: [] for p in representatives}
    status: dict[str, str] = {p.name: STATUS_REPRESENTATIVE for p in representatives}

    # Assign each non-representative to its closest representative via the skani
    # edge table (best ANI above the cutoffs). Streamed, so the full edge list is
    # never materialized -- the table can be large for big inputs.
    best: dict[str, tuple[str, float]] = {}  # member -> (rep, ani)
    for a, b, ani, af in _iter_edges(result_dir):
        if ani < ani_cutoff or af < af_cutoff:
            continue
        for member, rep in ((a, b), (b, a)):
            if rep in rep_names and member not in rep_names:
                if member not in best or ani > best[member][1]:
                    best[member] = (rep, ani)
    for member, (rep, _ani) in best.items():
        clusters[rep].append(member)
        status[member] = STATUS_CONTAINED

    # any remaining input genome that is neither a representative nor assigned
    for g in genomes:
        status.setdefault(g.name, STATUS_CONTAINED)

    return DerepResult(
        representatives=representatives,
        clusters=clusters,
        genome_status=status,
    )


def _find_representatives_dir(result_dir: Path) -> Path | None:
    candidates = [
        result_dir / "Dereplicated_Representative_Genomes",
        result_dir / "skDER_Dereplicated_Representative_Genomes",
    ]
    for cand in candidates:
        if cand.is_dir():
            return cand
    for sub in result_dir.rglob("*"):
        if sub.is_dir() and ("epresentative" in sub.name or "Dereplicated" in sub.name):
            if any(c.suffix in _FASTA_SUFFIXES for c in sub.iterdir()):
                return sub
    return None


def _iter_edges(result_dir: Path) -> Iterator[tuple[str, str, float, float]]:
    """Stream skDER's skani edge table as (a, b, ANI, min_AF) basename tuples."""
    edge_file = result_dir / "Skani_Triangle_Edge_Output.txt"
    if not edge_file.exists():
        matches = sorted(result_dir.rglob("*Edge_Output*.txt"))
        if not matches:
            return
        edge_file = matches[0]

    with open(edge_file) as fo:
        for ln, line in enumerate(fo):
            if ln == 0:
                continue  # header
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                ani = float(fields[2])
                af = min(float(fields[3]), float(fields[4]))
            except ValueError:
                continue
            yield (Path(fields[0]).name, Path(fields[1]).name, ani, af)
