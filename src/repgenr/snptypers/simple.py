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

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
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
    caps = _CAPABILITIES

    def tool(cmd, **kw):
        run_tool(caps, cmd, logger=logger, **kw)

    tool(
        ["minimap2", "-a", "-t", nt, ref, genome.resolve()],
        log_prefix="minimap2",
        stdout_path=sam,
    )
    tool(["samtools", "sort", "-@", nt, "-o", bam, sam], log_prefix="samtools")
    sam.unlink(missing_ok=True)
    tool(["samtools", "index", "-@", nt, bam], log_prefix="samtools")
    tool(
        ["bcftools", "mpileup", "--threads", nt, "-Ob", "-f", ref, "-o", pileup, bam],
        log_prefix=log,
    )
    # Haploid calling: the diploid default emits heterozygous genotypes on
    # bacterial genomes, which `bcftools consensus` renders as IUPAC codes
    # that Gubbins (and most alignment tools) reject.
    tool(
        ["bcftools", "call", "--threads", nt, "-mv", "--ploidy", "1", "-Ob", "-o", calls, pileup],
        log_prefix=log,
    )
    pileup.unlink(missing_ok=True)
    tool(
        ["bcftools", "view", "--threads", nt, "-v", "snps", "-Oz", "-o", snps, calls],
        log_prefix=log,
    )
    tool(["bcftools", "index", snps], log_prefix=log)
    tool(["bcftools", "consensus", "-f", ref, "-o", cons, snps], log_prefix=log)
    consensus = _concat_fasta(cons)
    # The consensus is in memory now; nothing downstream reads this genome's
    # intermediates. Keep them on failure, where they are the evidence.
    for leftover in (bam, Path(f"{bam}.bai"), calls, snps, Path(f"{snps}.csi"), cons):
        leftover.unlink(missing_ok=True)
    return consensus


def _concat_fasta(path: Path) -> str:
    parts: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(">"):
            parts.append(line.strip())
    return "".join(parts)


def _write_core_snps(consensuses: dict[str, str], core_fasta: Path, snp_matrix: Path) -> int:
    names = list(consensuses)
    seqs = [consensuses[n] for n in names]
    length = min(len(s) for s in seqs) if seqs else 0

    variable_cols = [col for col in range(length) if len({s[col] for s in seqs}) > 1]
    with open(core_fasta, "w", encoding="utf-8") as fo:
        for name, seq in zip(names, seqs, strict=True):
            snp_seq = "".join(seq[c] for c in variable_cols)
            fo.write(f">{name}\n")
            for pos in range(0, len(snp_seq), 80):
                fo.write(snp_seq[pos : pos + 80] + "\n")

    # pairwise SNP distance matrix
    snp_rows = {name: "".join(seqs[i][c] for c in variable_cols) for i, name in enumerate(names)}
    with open(snp_matrix, "w", encoding="utf-8") as fo:
        fo.write("\t" + "\t".join(names) + "\n")
        for a in names:
            dists = [
                str(sum(1 for x, y in zip(snp_rows[a], snp_rows[b], strict=True) if x != y))
                for b in names
            ]
            fo.write(a + "\t" + "\t".join(dists) + "\n")
    return len(variable_cols)
