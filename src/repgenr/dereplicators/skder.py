"""skDER dereplication adapter (skani-based).

skDER scales to tens of thousands of genomes natively (no chunking needed). It
estimates pairwise ANI/AF with skani and picks representatives in one round.

CLI used::

    skder -g <fofn or files> -o <out_dir> -i <ANI%> -f <AF%> -c <threads> \
          -d dynamic|greedy|low_mem_greedy

With ``-n`` skDER also assigns non-representatives to their closest
representative, which we parse into cluster membership.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.errors import ToolExecutionError, WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import run, write_fofn
from .base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
)


class SkderDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="skder",
        required_binaries=(
            BinarySpec("skder", version_args=("-v",)),
            BinarySpec("skani", version_args=("-V",)),
        ),
        default_params={"mode": "greedy"},
        recommended_max_genomes=None,  # scales natively
        supports_native_scaling=True,
        threads_param="-c",
    )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        fofn = write_fofn(genomes, out_dir / "genomes.fofn")
        result_dir = out_dir / "skder_out"

        mode = params.extra.get("mode", self.capabilities.default_params["mode"])
        cmd = [
            "skder",
            "-l", fofn,                       # genome list file
            "-o", result_dir,
            "-i", _as_percent(params.secondary_ani),
            "-f", _as_percent(params.aligned_fraction),
            "-c", params.threads,
            "-d", mode,
            "-n",                              # emit cluster membership
        ]
        run(cmd, logger=logger, log_prefix="skder")

        return _parse_skder_output(result_dir, genomes, logger)


def _as_percent(value: float) -> str:
    """skDER expects identity/fraction as a percentage (e.g. 99.0, 50.0)."""
    pct = value * 100 if value <= 1.0 else value
    return f"{pct:g}"


def _parse_skder_output(
    result_dir: Path, genomes: Sequence[Path], logger: logging.Logger
) -> DerepResult:
    rep_dir = _find_representatives_dir(result_dir)
    if rep_dir is None:
        raise WorkdirError(
            f"Could not locate skDER representative genomes under {result_dir}. "
            "Confirm skDER ran successfully."
        )
    representatives = sorted(
        p for p in rep_dir.iterdir() if p.suffix in (".fasta", ".fa", ".fna", ".fas")
    )
    if not representatives:
        raise ToolExecutionError(["skder"], 1, "skDER produced no representative genomes")

    by_basename = {p.name: p.name for p in genomes}

    clusters: dict[str, list[str]] = {p.name: [] for p in representatives}
    status: dict[str, str] = {p.name: STATUS_REPRESENTATIVE for p in representatives}

    cluster_file = _find_clustering_file(result_dir)
    if cluster_file is not None:
        for member, rep in _read_membership(cluster_file):
            member_name = by_basename.get(member, member)
            rep_name = by_basename.get(rep, rep)
            if rep_name not in clusters:
                clusters[rep_name] = []
            if member_name != rep_name:
                clusters[rep_name].append(member_name)
                status[member_name] = STATUS_CONTAINED
    else:
        logger.warning(
            "skDER clustering file not found; cluster membership will be representatives only"
        )

    # any input genome not seen is contained-but-unassigned; mark contained
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
    # fall back: a single subdir that looks like it holds FASTAs
    for sub in result_dir.rglob("*"):
        if sub.is_dir() and any(
            c.suffix in (".fasta", ".fa", ".fna", ".fas") for c in sub.iterdir()
        ):
            if "epresentative" in sub.name or "Dereplicated" in sub.name:
                return sub
    return None


def _find_clustering_file(result_dir: Path) -> Path | None:
    for pattern in ("*Clustering*.txt", "*Clustering*.tsv", "*clustering*.txt"):
        matches = sorted(result_dir.rglob(pattern))
        if matches:
            return matches[0]
    return None


def _read_membership(path: Path) -> list[tuple[str, str]]:
    """Return (member_basename, representative_basename) pairs.

    skDER clustering rows pair a genome path with its representative path; we key
    on basenames so they match the genome file inventory.
    """
    pairs: list[tuple[str, str]] = []
    with open(path) as fo:
        for ln, line in enumerate(fo):
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if ln == 0 and not fields[0].endswith((".fasta", ".fa", ".fna", ".fas")):
                continue  # header
            if len(fields) >= 2:
                member = Path(fields[0]).name
                rep = Path(fields[1]).name
                pairs.append((member, rep))
    return pairs
