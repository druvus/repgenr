"""Built-in 'simple' SNP typer (minimap2 + samtools + bcftools).

A lightweight reference-based core-SNP pipeline that needs no dedicated typing
tool. For each genome: align to the reference (minimap2, assembly preset asm20
by default), call SNP-only variants (bcftools), and build a SNP-only consensus
that preserves reference length. Reference positions that no primary or
supplementary alignment of the genome covers are set to N in its consensus, so
sequence a genome lacks is recorded as missing rather than as the reference
base. The per-genome consensuses (plus the reference) are stacked into a
whole-genome alignment, then reduced to the columns where at least two of A, C,
G and T occur to form the core-SNP alignment.

Genomes are independent, so they run concurrently: the thread budget becomes
``workers x threads-per-worker``, each worker driving its own chain of tools.
The per-genome intermediates (SAM, BAM, pileup, calls) are scratch. They are
written compressed and deleted as soon as that genome's consensus has been
read, which keeps peak scratch at a few hundred megabytes instead of growing
with the genome count.
"""

from __future__ import annotations

import gzip
import logging
import re
import shutil
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
    # minimap2 preset for assembly-to-reference mapping. asm20 tolerates the
    # divergence seen within a bacterial species (up to several percent);
    # ``--tool-arg preset=asm5`` suits near-identical genomes, and
    # ``preset=none`` restores minimap2's default (no -x).
    default_params={"preset": "asm20"},
    accepted_extras=frozenset({"preset"}),
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
        _copy_plain_fasta(reference, ref)
        run_tool(_CAPABILITIES, ["samtools", "faidx", ref], logger=logger, log_prefix="samtools")

        consensuses: dict[str, str] = {reference.stem: _concat_fasta(ref)}
        contigs = _reference_contigs(ref)
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
            lambda genome: _call_one(
                genome, ref, per_genome_dir, per_worker, params, logger, contigs=contigs
            ),
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
            raise WorkdirError(
                f"No variable sites found among {len(consensuses)} genomes against the "
                f"reference {reference.stem}, so no SNP tree can be built. The genomes "
                "either are identical at every position they share or did not align to "
                "the reference (unaligned positions are N and do not count). Try a "
                "closer reference (--reference), more divergent genomes, or an "
                "alignment-free tree (phylo --treebuilder mashtree)."
            )

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


def _call_one(
    genome: Path,
    ref: Path,
    work: Path,
    threads: int,
    params: SnpParams,
    logger,
    *,
    contigs: list[tuple[str, int]] | None = None,
) -> str:
    """Consensus of ``genome`` in reference coordinates, N where it does not align."""
    stem = genome.stem
    if contigs is None:
        contigs = _reference_contigs(ref)
    preset = str(params.extra.get("preset", _CAPABILITIES.default_params["preset"])).strip()
    preset_args = [] if preset.lower() in ("", "none") else ["-x", preset]
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
            (
                "minimap2",
                ["minimap2", "-a", *preset_args, "-t", nt, "-o", sam, ref, genome.resolve()],
            ),
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
    covered = _covered_positions(sam, contigs)
    consensus = _mask_uncovered(_concat_fasta(cons), covered, stem)
    # The consensus is in memory now; nothing downstream reads this genome's
    # intermediates. Keep them on failure, where they are the evidence.
    for leftover in (sam, bam, Path(f"{bam}.bai"), pileup, calls, snps, Path(f"{snps}.csi"), cons):
        leftover.unlink(missing_ok=True)
    return consensus


def _copy_plain_fasta(src: Path, dest: Path) -> None:
    """Copy ``src`` to ``dest`` uncompressed.

    The reference copy is indexed with ``samtools faidx`` and read here as
    text, so a gzipped genome (judged by its magic bytes, not its name) is
    decompressed. Query genomes are passed to minimap2 as they are, since
    minimap2 reads gzipped FASTA itself.
    """
    with open(src, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    opener = gzip.open if gzipped else open
    with opener(src, "rb") as fi, open(dest, "wb") as fo:
        shutil.copyfileobj(fi, fo, 1 << 20)


def _reference_contigs(ref: Path) -> list[tuple[str, int]]:
    """(name, length) of each reference record, in file order.

    This is the order in which ``bcftools consensus`` writes the records and
    hence the order of the concatenated consensus.
    """
    contigs: list[tuple[str, int]] = []
    name: str | None = None
    length = 0
    with open(ref, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    contigs.append((name, length))
                name = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        contigs.append((name, length))
    return contigs


# CIGAR operations that consume the reference: M, D, N, = and X.
_CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
_REF_CONSUMING = frozenset("MDN=X")
# SAM flags of records that do not describe where the genome aligns:
# unmapped (0x4) and secondary (0x100).
_SKIP_FLAGS = 0x4 | 0x100


def _covered_positions(sam: Path, contigs: list[tuple[str, int]]) -> npt.NDArray[np.bool_]:
    """Reference positions spanned by the genome's alignments, over the concatenated contigs.

    Primary and supplementary records count; secondary and unmapped records do
    not. A record covers its reference span from POS to the end given by the
    CIGAR, deletions included, since the consensus there is defined by the
    alignment.
    """
    offsets: dict[str, int] = {}
    total = 0
    for name, length in contigs:
        offsets[name] = total
        total += length
    covered = np.zeros(total, dtype=np.bool_)
    with open(sam, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.split("\t", 6)
            if len(fields) < 6 or int(fields[1]) & _SKIP_FLAGS:
                continue
            offset = offsets.get(fields[2])
            if offset is None or fields[5] == "*":
                continue
            span = sum(int(n) for n, op in _CIGAR_OP.findall(fields[5]) if op in _REF_CONSUMING)
            start = offset + int(fields[3]) - 1
            covered[start : start + span] = True
    return covered


def _mask_uncovered(consensus: str, covered: npt.NDArray[np.bool_], name: str) -> str:
    if len(consensus) != len(covered):
        raise WorkdirError(
            f"The consensus of {name} has {len(consensus)} positions but the reference has "
            f"{len(covered)}; the simple typer expects a SNP-only consensus of reference length."
        )
    seq = np.frombuffer(consensus.encode("ascii", "replace"), dtype=np.uint8).copy()
    seq[~covered] = ord("N")
    return seq.tobytes().decode("ascii")


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


# Byte to base code: A, C, G, T (either case) to 0..3, everything else
# (N, IUPAC codes, gaps) to 4, which marks a position without a base.
_NO_BASE = 4
_BASE_CODE = np.full(256, _NO_BASE, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _BASE_CODE[ord(_b)] = _i
    _BASE_CODE[ord(_b.lower())] = _i
# The same codes with "no base" at 0, so a column maximum ignores missing data.
_BASE_CODE_LOW = np.where(_BASE_CODE == _NO_BASE, 0, _BASE_CODE).astype(np.uint8)
_CODE_BASE = np.frombuffer(b"ACGTN", dtype=np.uint8)


def _core_columns(seqs: list[str], length: int) -> npt.NDArray[np.uint8]:
    """Base codes (0..3 for A, C, G, T; 4 for no base) of the variable columns of ``seqs``.

    A column is variable when at least two of A, C, G and T occur in it. N,
    other ambiguity codes and gaps are missing data: a column holding them and
    one base is not variable. Columns where every genome holds the same byte
    are dropped first with one comparison, so the base-level test runs only on
    the few columns that differ at all.
    """
    kept: list[npt.NDArray[np.uint8]] = []
    for start in range(0, length, _COLUMN_BLOCK):
        end = min(start + _COLUMN_BLOCK, length)
        block = np.frombuffer(
            b"".join(s[start:end].encode("ascii", "replace") for s in seqs), dtype=np.uint8
        ).reshape(len(seqs), end - start)
        differing = block[:, (block != block[0]).any(axis=0)]
        if differing.shape[1] == 0:
            continue
        # Smallest base code with missing data high, largest with it low: the
        # two differ only when the column holds two different bases.
        lowest = _BASE_CODE[differing].min(axis=0)
        highest = _BASE_CODE_LOW[differing].max(axis=0)
        variable = highest > lowest
        if variable.any():
            kept.append(_BASE_CODE[differing[:, variable]])
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
            snp_seq = _CODE_BASE[row].tobytes().decode("ascii")
            fo.write(f">{name}\n")
            for pos in range(0, n_sites, 80):
                fo.write(snp_seq[pos : pos + 80] + "\n")

    # Pairwise SNP distances over those columns, counting only sites where both
    # genomes have a base. Each row is compared against the whole matrix at
    # once; the work is still quadratic in genomes, but each pair costs a
    # vector comparison instead of a Python loop.
    has_base = core != _NO_BASE
    with open(snp_matrix, "w", encoding="utf-8") as fo:
        fo.write("\t" + "\t".join(names) + "\n")
        for name, row, row_base in zip(names, core, has_base, strict=True):
            dists = ((core != row) & has_base & row_base).sum(axis=1)
            fo.write(name + "\t" + "\t".join(str(int(d)) for d in dists) + "\n")
    return n_sites
