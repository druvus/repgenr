"""The per-stage input specs feed the v2 fingerprint (STAGE_INPUTS table)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from repgenr.cli import base as cli
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import CLUSTERS_TSV, CORE_SNP_FASTA, TREE_NWK
from repgenr.core.inputs import ABSENT


def _digests(ctx: WorkdirContext, stage: str, params) -> dict[str, str]:
    return cli._stage_input_digests(ctx, stage, params)


def test_dereplicate_inputs_track_genomes_dir(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    ctx.genomes_dir.mkdir()
    (ctx.genomes_dir / "g1.fasta").write_text(">g1\nACGT\n", encoding="utf-8")
    params = SimpleNamespace()
    before = _digests(ctx, "dereplicate", params)
    assert set(before) == {"genomes"}

    (ctx.genomes_dir / "g2.fasta").write_text(">g2\nACGT\n", encoding="utf-8")
    after = _digests(ctx, "dereplicate", params)
    assert after["genomes"] != before["genomes"]


def test_phylo_inputs_conditional_on_msa_source(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    aligner = _digests(ctx, "phylo", SimpleNamespace(all_genomes=False, msa_source="aligner"))
    snptype = _digests(ctx, "phylo", SimpleNamespace(all_genomes=False, msa_source="snptype"))
    core_snp_key = f"snp/{CORE_SNP_FASTA}"
    assert core_snp_key not in aligner
    assert core_snp_key in snptype
    assert "derep/representatives" in aligner and "outgroup" in aligner


def test_phylo_all_genomes_switches_input_dir(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    reps = _digests(ctx, "phylo", SimpleNamespace(all_genomes=False, msa_source="aligner"))
    allg = _digests(ctx, "phylo", SimpleNamespace(all_genomes=True, msa_source="aligner"))
    assert "derep/representatives" in reps and "genomes" not in reps
    assert "genomes" in allg and "derep/representatives" not in allg


def test_tree2tax_inputs_track_tree_manifest_and_clusters(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    ctx.tree_dir.mkdir()
    tree = ctx.tree_dir / TREE_NWK
    tree.write_text("(a,b);\n", encoding="utf-8")

    plain = _digests(ctx, "tree2tax", SimpleNamespace(include_dereplicated=False))
    assert set(plain) == {f"tree/{TREE_NWK}", "manifest"}

    withderep = _digests(ctx, "tree2tax", SimpleNamespace(include_dereplicated=True))
    clusters_key = f"derep/{CLUSTERS_TSV}"
    assert clusters_key in withderep
    assert withderep[clusters_key] == ABSENT  # not written yet -> stable sentinel

    tree.write_text("(a,b,c);\n", encoding="utf-8")
    changed = _digests(ctx, "tree2tax", SimpleNamespace(include_dereplicated=False))
    assert changed[f"tree/{TREE_NWK}"] != plain[f"tree/{TREE_NWK}"]


def test_genome_inputs_track_selection_tsv(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    sel = ctx.workdir / "selection.tsv"
    sel.write_text("accession\nGCF_1\n", encoding="utf-8")
    before = _digests(ctx, "genome", SimpleNamespace())
    sel.write_text("accession\nGCF_1\nGCF_2\n", encoding="utf-8")
    after = _digests(ctx, "genome", SimpleNamespace())
    assert after["selection.tsv"] != before["selection.tsv"]


def test_metadata_has_no_local_inputs(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    assert _digests(ctx, "metadata", SimpleNamespace()) == {}
