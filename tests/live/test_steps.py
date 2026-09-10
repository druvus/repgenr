"""Stateless data-channel steps run as processes: chunk, merge, phylo-build,
tree2tax-relations, and --versions-out on each."""

from __future__ import annotations

from pathlib import Path

import pytest
from live_helpers import newick_leaves, truth_of, truth_partition

from repgenr.core.contracts import (
    CLUSTERS_TSV,
    GENOMES_MAP_TSV,
    TREE2TAX_TSV,
    TREE_NWK,
    read_clusters,
)

pytestmark = [pytest.mark.live, pytest.mark.requires_binary("sourmash", "mashtree")]


def _partition(clusters_tsv: Path):
    return {frozenset({rep, *m}) for rep, m in read_clusters(clusters_tsv).items()}


@pytest.fixture
def chunks(synthetic_set, run_repgenr, tmp_path: Path):
    """Three chunk results over a clonal n=12 set, plus the set itself."""
    genomes = synthetic_set("clonal", n=12, length=50_000)
    files = sorted(genomes.glob("*.fasta"))
    outs = []
    for i in range(3):
        fofn = tmp_path / f"chunk{i}.fofn"
        fofn.write_text("".join(f"{f}\n" for f in files[i * 4 : (i + 1) * 4]), encoding="utf-8")
        out = tmp_path / f"chunk{i}"
        run_repgenr(
            "dereplicate-chunk",
            "--genomes-fofn",
            fofn,
            "-o",
            out,
            "--tool",
            "sourmash",
            "-t",
            "2",
            "--versions-out",
            out / "versions.yml",
        )
        outs.append(out)
    return genomes, outs


def test_chunk_results_carry_the_contract_and_versions(chunks) -> None:
    _, outs = chunks
    for out in outs:
        assert (out / CLUSTERS_TSV).is_file() and (out / "representatives").is_dir()
        assert "sourmash:" in (out / "versions.yml").read_text(encoding="utf-8")


def test_merge_by_chunk_dir_recovers_the_partition(run_repgenr, chunks, tmp_path: Path) -> None:
    genomes, outs = chunks
    merged = tmp_path / "merged"
    run_repgenr(
        "dereplicate-merge",
        "-o",
        merged,
        "--tool",
        "sourmash",
        "-t",
        "2",
        *[a for out in outs for a in ("--chunk-dir", str(out))],
    )
    reps = set(read_clusters(merged / CLUSTERS_TSV))
    expected = truth_partition(truth_of(genomes))
    # The merge sees only chunk representatives; every truth cluster must be
    # represented exactly once among them.
    assert len(reps) == len(expected)
    by_cluster = {name: c for name, c in truth_of(genomes)["clusters"].items()}
    assert len({by_cluster[r] for r in reps}) == len(expected)


def test_merge_by_chunk_fofn(run_repgenr, chunks, tmp_path: Path) -> None:
    _, outs = chunks
    fofn = tmp_path / "chunks.fofn"
    fofn.write_text("".join(f"{o}\n" for o in outs), encoding="utf-8")
    merged = tmp_path / "merged_fofn"
    run_repgenr(
        "dereplicate-merge",
        "-o",
        merged,
        "--chunk-fofn",
        fofn,
        "--tool",
        "sourmash",
        "-t",
        "2",
        "--versions-out",
        merged / "versions.yml",
    )
    assert (merged / CLUSTERS_TSV).is_file()
    assert "sourmash:" in (merged / "versions.yml").read_text(encoding="utf-8")


def test_phylo_build_and_tree2tax_relations_with_outgroup(
    run_repgenr, synthetic_set, tmp_path: Path
) -> None:
    genomes = synthetic_set("balanced", n=5, length=50_000)
    files = sorted(genomes.glob("*.fasta"))
    outgroup_dir = tmp_path / "outgroup"
    outgroup_dir.mkdir()
    og = files[-1]
    (outgroup_dir / og.name).symlink_to(og)
    ingroup_dir = tmp_path / "ingroup"
    ingroup_dir.mkdir()
    for f in files[:-1]:
        (ingroup_dir / f.name).symlink_to(f)
    acc_file = tmp_path / "outgroup_accession.txt"
    acc_file.write_text(og.stem.rsplit("_", 1)[-1] + "\n", encoding="utf-8")

    out = tmp_path / "phylo"
    run_repgenr(
        "phylo-build",
        "--genomes-dir",
        ingroup_dir,
        "-o",
        out,
        "--outgroup-dir",
        outgroup_dir,
        "--outgroup-accession",
        acc_file,
        "--treebuilder",
        "mashtree",
        "-t",
        "2",
        "--versions-out",
        out / "versions.yml",
    )
    tree = (out / "tree" / TREE_NWK).read_text(encoding="utf-8")
    assert newick_leaves(tree) == {f.stem for f in files}, "outgroup leaf joins the ingroup"
    assert "mashtree:" in (out / "versions.yml").read_text(encoding="utf-8")

    clusters = tmp_path / CLUSTERS_TSV
    rows = ["representative\tmember"]
    for f in files[:-1]:
        rows.append(f"{f.name}\t{f.name}")
    rows.append(f"{files[0].name}\tredundant_GCF9999999.1.fasta")
    clusters.write_text("\n".join(rows) + "\n", encoding="utf-8")
    rel = tmp_path / "relations"
    run_repgenr(
        "tree2tax-relations",
        "--tree",
        out / "tree" / TREE_NWK,
        "-o",
        rel,
        "--clusters",
        clusters,
        "--include-dereplicated",
        "--outgroup-dir",
        outgroup_dir,
        "--outgroup-accession",
        acc_file,
        "--remove-outgroup",
        "--node-basename",
        "N",
        "-r",
        "top",
        "--versions-out",
        rel / "versions.yml",
    )
    relations = (rel / TREE2TAX_TSV).read_text(encoding="utf-8")
    assert "\ttop\n" in relations and og.stem not in relations, "outgroup dropped, root renamed"
    assert "\tN" in relations or "N1" in relations, "internal nodes carry the basename"
    mapped = (rel / GENOMES_MAP_TSV).read_text(encoding="utf-8")
    assert "GCF9999999.1" in mapped, "--include-dereplicated maps the redundant genome"


def test_tree2tax_relations_collapse_flags(run_repgenr, tmp_path: Path) -> None:
    """Weak splits collapse by support or by branch length (pure Python step)."""
    tree = tmp_path / "t.nwk"
    tree.write_text(
        "(((A:0.1,B:0.1)40:0.001,C:0.1)95:0.2,(D:0.1,E:0.1)99:0.2);\n", encoding="utf-8"
    )

    def rows(out: Path) -> list[str]:
        return (out / TREE2TAX_TSV).read_text(encoding="utf-8").splitlines()[1:]

    plain = tmp_path / "plain"
    run_repgenr("tree2tax-relations", "--tree", tree, "-o", plain)
    by_support = tmp_path / "support"
    run_repgenr("tree2tax-relations", "--tree", tree, "-o", by_support, "--collapse-support", "0.5")
    by_length = tmp_path / "length"
    run_repgenr("tree2tax-relations", "--tree", tree, "-o", by_length, "--collapse-length", "0.01")
    assert len(rows(by_support)) == len(rows(plain)) - 1, (
        "the support-40 node merges into its parent"
    )
    assert len(rows(by_length)) == len(rows(plain)) - 1, "the 0.001 branch merges into its parent"
