"""Gzipped genomes reach tools that cannot read gzip as decompressed copies.

SibeliaZ (TwoPaCo) and ParSNP fail on gzipped FASTA (an empty MAF, and a
UnicodeDecodeError); the stages give them plain copies from scratch. Records
and leaves carry the record name (no '.fasta'), and the copies are removed.
"""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest
from live_helpers import newick_leaves

from repgenr.core.contracts import MSA_FASTA, TREE_NWK, record_name

pytestmark = [pytest.mark.live]


@pytest.fixture
def gz_genomes(synthetic_set, tmp_path: Path) -> Path:
    """Four gzipped genomes (phylo needs three; ParSNP's RAxML step needs four)."""
    plain = synthetic_set("balanced", n=4, length=50_000)
    out = tmp_path / "gz"
    out.mkdir()
    for genome in sorted(plain.glob("*.fasta")):
        with open(genome, "rb") as fi, gzip.open(out / f"{genome.name}.gz", "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return out


def _names(genomes: Path) -> set[str]:
    return {record_name(p) for p in genomes.glob("*.fasta.gz")}


def _headers(fasta: Path) -> set[str]:
    return {
        line[1:].split()[0]
        for line in fasta.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    }


@pytest.mark.requires_binary("sibeliaz", "FastTree")
def test_sibeliaz_aligns_gzipped_genomes(run_repgenr, gz_genomes: Path, tmp_path: Path) -> None:
    out = tmp_path / "sib"
    run_repgenr(
        "phylo-build",
        "--genomes-dir",
        gz_genomes,
        "-o",
        out,
        "--no-outgroup",
        "--msa-source",
        "aligner",
        "--aligner",
        "sibeliaz",
        "--treebuilder",
        "fasttree",
        "-t",
        "2",
    )
    assert _headers(out / "align" / MSA_FASTA) == _names(gz_genomes)
    assert newick_leaves((out / "tree" / TREE_NWK).read_text(encoding="utf-8")) == _names(
        gz_genomes
    )
    assert not (out / "scratch" / "phylo_inputs").exists()


@pytest.mark.requires_binary("parsnp", "harvesttools", "FastTree")
def test_parsnp_types_gzipped_genomes(run_repgenr, gz_genomes: Path, tmp_path: Path) -> None:
    out = tmp_path / "parsnp"
    run_repgenr(
        "phylo-build",
        "--genomes-dir",
        gz_genomes,
        "-o",
        out,
        "--no-outgroup",
        "--msa-source",
        "snptype",
        "--snptyper",
        "parsnp",
        "--treebuilder",
        "fasttree",
        "-t",
        "2",
    )
    assert newick_leaves((out / "tree" / TREE_NWK).read_text(encoding="utf-8")) == _names(
        gz_genomes
    )
    assert not list((out / "scratch").rglob("inputs/*"))
