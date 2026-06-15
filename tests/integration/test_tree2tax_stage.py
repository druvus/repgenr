"""tree2tax stage test against a synthetic rooted tree (uses real ete3)."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.context import WorkdirContext
from repgenr.stages.tree2tax import Tree2taxParams, run


def _setup(workdir: Path) -> WorkdirContext:
    # tree leaves are genome stems; include an outgroup leaf
    tree_dir = workdir / "tree"
    tree_dir.mkdir(parents=True)
    (tree_dir / "tree.nwk").write_text(
        "((Fam_gen_sp_GCA_000001:0.1,Fam_gen_sp_GCA_000002:0.1):0.2,"
        "Out_gen_sp_GCA_000099:0.5);\n"
    )
    # outgroup metadata
    (workdir / "outgroup_accession.txt").write_text("GCA_000099\n")
    og = workdir / "outgroup"
    og.mkdir()
    (og / "Out_gen_sp_GCA_000099.fasta").write_text(">x\nACGT\n")

    # derep clusters: rep GCA_000001 contains a redundant GCA_000003
    derep = workdir / "derep"
    derep.mkdir()
    (derep / "clusters.tsv").write_text(
        "representative\tmember\n"
        "Fam_gen_sp_GCA_000001.fasta\tFam_gen_sp_GCA_000001.fasta\n"
        "Fam_gen_sp_GCA_000001.fasta\tFam_gen_sp_GCA_000003.fasta\n"
        "Fam_gen_sp_GCA_000002.fasta\tFam_gen_sp_GCA_000002.fasta\n"
    )
    return WorkdirContext(workdir, create=True)


def test_tree2tax_outputs(workdir: Path) -> None:
    ctx = _setup(workdir)
    t2t, gmap = run(ctx, Tree2taxParams(include_dereplicated=True, remove_outgroup=True))

    edges = [line.split("\t") for line in t2t.read_text().splitlines()[1:]]
    children = {c for c, _ in edges}
    assert "Fam_gen_sp_GCA_000001" in children
    assert "Fam_gen_sp_GCA_000002" in children
    # outgroup removed from relations
    assert "Out_gen_sp_GCA_000099" not in children

    rows = [line.split("\t") for line in gmap.read_text().splitlines()]
    mapping = {acc: leaf for acc, leaf in rows}
    # representative maps to its own leaf
    assert mapping["GCA_000001"] == "Fam_gen_sp_GCA_000001"
    # redundant genome maps to the representative leaf
    assert mapping["GCA_000003"] == "Fam_gen_sp_GCA_000001"
