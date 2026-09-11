"""Built-in 'simple' SNP typer (minimap2 + samtools + bcftools).

A lightweight reference-based core-SNP pipeline that needs no dedicated typing
tool. For each genome: align to the reference (minimap2), call SNP-only variants
(bcftools), and build a SNP-only consensus that preserves reference length. The
per-genome consensuses (plus the reference) are stacked into a whole-genome
alignment, then reduced to variable columns to form the core-SNP alignment.

Genomes are independent, so they run concurrently: the thread budget becomes
``workers x threads-per-worker``, each worker driving its own chain of tools.
The per-genome intermediates (SAM, BAM, pileup, calls) are scratch. They are
written compressed and deleted as soon as that genome's consensus has been
read, which keeps peak scratch at a few hundred megabytes instead of growing
with the genome count.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ..core.binaries import BinarySpec
from ..core.containers import run_chain, run_tool
from ..core.errors import WorkdirError
from ..core.executors import parallel_map
from ..core.plugins import ToolCapabilities
from .base import SnpParams, SnpResult, SnpTyper

_CAPABILITIES = ToolCapabilities(
    name="simple",
    required_binaries=(
        BinarySpec("minimap2", version_args=("--version",), min_version="2.17"),
        # samtools/bcftools >= 1.10: an ancient 0.1.x (pulled in by some perl
        # deps) lacks `bcftools mpileup` and breaks the SNP pipeline. strict_version
        # so an ancient build that does not answer --version is rejected, not skipped.
        BinarySpec(
            "samtools", version_args=("--version",), min_version="1.10", strict_version=True
        ),
        BinarySpec(
            "bcftools", version_args=("--version",), min_version="1.10", strict_version=True
        ),
    ),
    recommended_max_genomes=2000,
    # multi-tool: resolved to one image via Wave (or pin an explicit container)
    conda=("bioconda::minimap2", "bioconda::samtools", "bioconda::bcftools"),
)


class SimpleSnpTyper(SnpTyper):
    capabilities = _CAPABILITIES
    requires_reference = True

    def call(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: SnpParams,
        logger: logging.Logger,
    ) -> SnpResult:
        genomes = list(genomes)
        if reference is None:
            reference = genomes[0]
        out_dir.mkdir(parents=True, exist_ok=True)

        ref = out_dir / "reference.fasta"
        ref.write_text(reference.read_text())
        run_tool(_CAPABILITIES, ["samtools", "faidx", ref], logger=logger, log_prefix="samtools")

        consensuses: dict[str, str] = {reference.stem: _concat_fasta(ref)}
        per_genome_dir = out_dir / "per_genome"
        per_genome_dir.mkdir(exist_ok=True)

        others = [g for g in genomes if g.resolve() != reference.resolve()]
        workers, per_worker = _split_threads(params.threads, len(others))
        if workers > 1:
            logger.info(
                "Calling %d genomes on %d worker(s), %d thread(s) each",
                len(others),
                workers,
                per_worker,
            )
        called = parallel_map(
            lambda genome: _call_one(genome, ref, per_genome_dir, per_worker, params, logger),
            others,
            workers,
            logger=logger,
        )
        consensuses.update(zip((g.stem for g in others), called, strict=True))

        full_fasta = out_dir / "full_alignment.fasta"
        with open(full_fasta, "w", encoding="utf-8") as fo:
            for name, seq in consensuses.items():
                fo.write(f">{name}\n{seq}\n")

        core_fasta = out_dir / "core_snp.fasta"
        snp_matrix = out_dir / "snp_distance_matrix.tsv"
        n_sites = _write_core_snps(consensuses, core_fasta, snp_matrix)
        logger.info(
            "simple SNP typer: %d core SNP sites across %d genomes", n_sites, len(consensuses)
        )
        if n_sites == 0:
            raise WorkdirError("No variable sites found; cannot build a SNP tree.")

        return SnpResult(
            core_snp_fasta=core_fasta,
            snp_distance_matrix=snp_matrix,
            masked=False,
            full_alignment=full_fasta,
        )


def _split_threads(threads: int, genomes: int) -> tuple[int, int]:
    """Workers and threads each, for ``genomes`` independent per-genome chains.

    One genome's chain is dominated by single-threaded steps (mpileup above
    all), so the budget buys far more as concurrent genomes than as threads
    inside one chain. Threads go to the tools only once there are more than
    enough workers for the genomes at hand.
    """
    if threads < 1 or genomes < 1:
        return 1, 1
    workers = max(1, min(threads, genomes))
    return workers, max(1, threads // workers)


def _call_one(genome: Path, ref: Path, work: Path, threads: int, params: SnpParams, logger) -> str:
    stem = genome.stem
    sam = work / f"{stem}.sam"
    bam = work / f"{stem}.bam"
    # Compressed BCF, not plain VCF: the pileup of a bacterial genome is a
    # record per reference base, and it is read once by `bcftools call`.
    pileup = work / f"{stem}.pileup.bcf"
    calls = work / f"{stem}.calls.bcf"
    snps = work / f"{stem}.snps.vcf.gz"
    cons = work / f"{stem}.consensus.fasta"

    nt = str(threads)
    log = "bcftools"

    # One unit of work: these run in order, and in a single container when one
    # is in use. minimap2 writes with -o rather than to stdout so that every
    # step is a plain argument vector.
    run_chain(
        _CAPABILITIES,
        [
            ("minimap2", ["minimap2", "-a", "-t", nt, "-o", sam, ref, genome.resolve()]),
            ("samtools", ["samtools", "sort", "-@", nt, "-o", bam, sam]),
            ("samtools", ["samtools", "index", "-@", nt, bam]),
            (
                "bcftools",
                ["bcftools", "mpileup", "--threads", nt, "-Ob", "-f", ref, "-o", pileup, bam],
            ),
            # Haploid calling: the diploid default emits heterozygous genotypes
            # on bacterial genomes, which `bcftools consensus` renders as IUPAC
            # codes that Gubbins (and most alignment tools) reject.
            (
                log,
                [
                    "bcftools",
                    "call",
                    "--threads",
                    nt,
                    "-mv",
                    "--ploidy",
                    "1",
                    "-Ob",
                    "-o",
                    calls,
                    pileup,
                ],
            ),
            (log, ["bcftools", "view", "--threads", nt, "-v", "snps", "-Oz", "-o", snps, calls]),
            (log, ["bcftools", "index", snps]),
            (log, ["bcftools", "consensus", "-f", ref, "-o", cons, snps]),
        ],
        logger=logger,
        extra_mounts=[genome.resolve().parent, work, ref.parent],
    )
    consensus = _concat_fasta(cons)
    # The consensus is in memory now; nothing downstream reads this genome's
    # intermediates. Keep them on failure, where they are the evidence.
    for leftover in (sam, bam, Path(f"{bam}.bai"), pileup, calls, snps, Path(f"{snps}.csi"), cons):
        leftover.unlink(missing_ok=True)
    return consensus


def _concat_fasta(path: Path) -> str:
    parts: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(">"):
            parts.append(line.strip())
    return "".join(parts)


# Columns per pass when reducing the stacked consensuses. One block of every
# genome is held as bytes at a time, so peak memory follows the genome count
# rather than the reference length.
_COLUMN_BLOCK = 1 << 20


def _core_columns(seqs: list[str], length: int) -> npt.NDArray[np.uint8]:
    """Stack the variable columns of ``seqs`` as a (genomes x sites) byte array.

    A column is variable when any genome differs from the first, which is what
    comparing every pair in the column amounts to. Characters are compared as
    bytes, so an N or a gap counts as a difference exactly as before.
    """
    kept: list[npt.NDArray[np.uint8]] = []
    for start in range(0, length, _COLUMN_BLOCK):
        end = min(start + _COLUMN_BLOCK, length)
        block = np.frombuffer(
            b"".join(s[start:end].encode("ascii", "replace") for s in seqs), dtype=np.uint8
        ).reshape(len(seqs), end - start)
        varying = (block != block[0]).any(axis=0)
        if varying.any():
            kept.append(block[:, varying])
    if not kept:
        return np.empty((len(seqs), 0), dtype=np.uint8)
    return np.concatenate(kept, axis=1)


def _write_core_snps(consensuses: dict[str, str], core_fasta: Path, snp_matrix: Path) -> int:
    names = list(consensuses)
    seqs = [consensuses[n] for n in names]
    length = min(len(s) for s in seqs) if seqs else 0

    core = _core_columns(seqs, length)
    n_sites = int(core.shape[1])
    with open(core_fasta, "w", encoding="utf-8") as fo:
        for name, row in zip(names, core, strict=True):
            snp_seq = row.tobytes().decode("ascii")
            fo.write(f">{name}\n")
            for pos in range(0, n_sites, 80):
                fo.write(snp_seq[pos : pos + 80] + "\n")

    # Pairwise SNP distances over those columns. Each row is compared against
    # the whole matrix at once; the work is still quadratic in genomes, but
    # each pair costs a vector comparison instead of a Python loop.
    with open(snp_matrix, "w", encoding="utf-8") as fo:
        fo.write("\t" + "\t".join(names) + "\n")
        for name, row in zip(names, core, strict=True):
            dists = (core != row).sum(axis=1)
            fo.write(name + "\t" + "\t".join(str(int(d)) for d in dists) + "\n")
    return n_sites
