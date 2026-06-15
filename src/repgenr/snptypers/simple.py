"""Built-in 'simple' SNP typer (minimap2 + samtools + bcftools).

A lightweight reference-based core-SNP pipeline that needs no dedicated typing
tool. For each genome: align to the reference (minimap2), call SNP-only variants
(bcftools), and build a SNP-only consensus that preserves reference length. The
per-genome consensuses (plus the reference) are stacked into a whole-genome
alignment, then reduced to variable columns to form the core-SNP alignment.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import run
from .base import SnpParams, SnpResult, SnpTyper


class SimpleSnpTyper(SnpTyper):
    capabilities = ToolCapabilities(
        name="simple",
        required_binaries=(
            BinarySpec("minimap2", version_args=("--version",)),
            BinarySpec("samtools", version_args=("--version",)),
            BinarySpec("bcftools", version_args=("--version",)),
        ),
        recommended_max_genomes=2000,
        threads_param=None,
    )
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
        run(["samtools", "faidx", ref], logger=logger, log_prefix="samtools")

        consensuses: dict[str, str] = {reference.stem: _concat_fasta(ref)}
        per_genome_dir = out_dir / "per_genome"
        per_genome_dir.mkdir(exist_ok=True)

        for genome in genomes:
            if genome.resolve() == reference.resolve():
                continue
            consensuses[genome.stem] = _call_one(genome, ref, per_genome_dir, params, logger)

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
        )


def _call_one(genome: Path, ref: Path, work: Path, params: SnpParams, logger) -> str:
    stem = genome.stem
    sam = work / f"{stem}.sam"
    bam = work / f"{stem}.bam"
    calls = work / f"{stem}.calls.vcf"
    snps = work / f"{stem}.snps.vcf.gz"
    cons = work / f"{stem}.consensus.fasta"

    pileup = work / f"{stem}.pileup.vcf"
    log = "bcftools"
    run(
        ["minimap2", "-a", ref, genome.resolve()],
        logger=logger, log_prefix="minimap2", stdout_path=sam,
    )
    run(["samtools", "sort", "-o", bam, sam], logger=logger, log_prefix="samtools")
    run(["samtools", "index", bam], logger=logger, log_prefix="samtools")
    run(["bcftools", "mpileup", "-f", ref, "-o", pileup, bam], logger=logger, log_prefix=log)
    run(["bcftools", "call", "-mv", "-Ov", "-o", calls, pileup], logger=logger, log_prefix=log)
    run(
        ["bcftools", "view", "-v", "snps", "-Oz", "-o", snps, calls],
        logger=logger, log_prefix=log,
    )
    run(["bcftools", "index", snps], logger=logger, log_prefix=log)
    run(["bcftools", "consensus", "-f", ref, "-o", cons, snps], logger=logger, log_prefix=log)
    return _concat_fasta(cons)


def _concat_fasta(path: Path) -> str:
    parts: list[str] = []
    for line in path.read_text().splitlines():
        if not line.startswith(">"):
            parts.append(line.strip())
    return "".join(parts)


def _write_core_snps(consensuses: dict[str, str], core_fasta: Path, snp_matrix: Path) -> int:
    names = list(consensuses)
    seqs = [consensuses[n] for n in names]
    length = min(len(s) for s in seqs) if seqs else 0

    variable_cols = [
        col for col in range(length) if len({s[col] for s in seqs}) > 1
    ]
    with open(core_fasta, "w") as fo:
        for name, seq in zip(names, seqs, strict=True):
            snp_seq = "".join(seq[c] for c in variable_cols)
            fo.write(f">{name}\n")
            for pos in range(0, len(snp_seq), 80):
                fo.write(snp_seq[pos : pos + 80] + "\n")

    # pairwise SNP distance matrix
    snp_rows = {
        name: "".join(seqs[i][c] for c in variable_cols) for i, name in enumerate(names)
    }
    with open(snp_matrix, "w") as fo:
        fo.write("\t" + "\t".join(names) + "\n")
        for a in names:
            dists = [
                str(sum(1 for x, y in zip(snp_rows[a], snp_rows[b], strict=True) if x != y))
                for b in names
            ]
            fo.write(a + "\t" + "\t".join(dists) + "\n")
    return len(variable_cols)
