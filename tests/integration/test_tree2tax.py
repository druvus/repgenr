"""tree2tax on a small tree (dendropy port): edges, genomes_map, outgroup rooting."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.context import WorkdirContext
from repgenr.stages.tree2tax import Tree2taxParams
from repgenr.stages.tree2tax import run as tree2tax_run

_NWK = (
    "((Fam_Gen_sp_GCA_000001.1,Fam_Gen_sp_GCA_000002.1),"
    "(Fam_Gen_sp_GCA_000003.1,Fam_Gen_sp_GCA_000004.1));"
)


def _edges(path: Path) -> set[tuple[str, str]]:
    out = set()
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        c, p = line.split("\t")
        out.add((c, p))
    return out


def _setup(workdir: Path, *, with_outgroup: bool = False, clusters: dict | None = None) -> None:
    (workdir / "tree").mkdir(parents=True)
    (workdir / "tree" / "tree.nwk").write_text(_NWK + "\n")
    if clusters is not None:
        (workdir / "derep").mkdir(parents=True)
        from repgenr.core.contracts import write_clusters
        write_clusters(workdir / "derep" / "clusters.tsv", clusters)
    if with_outgroup:
        og = workdir / "outgroup"
        og.mkdir(parents=True)
        (og / "Fam_Gen_sp_GCA_000004.1.fasta").write_text(">x\nACGT\n")
        (workdir / "outgroup_accession.txt").write_text("GCA_000004.1\n")


def test_tree2tax_unrooted_edges_and_map(workdir: Path) -> None:
    _setup(workdir)
    ctx = WorkdirContext(workdir, create=True)
    t2t, gmap = tree2tax_run(ctx, Tree2taxParams())
    edges = _edges(t2t)
    # every leaf has a path to root; two cherries -> two distinct internal nodes
    leaves = {c for c, _ in edges if c.startswith("Fam_")}
    assert leaves == {
        "Fam_Gen_sp_GCA_000001.1", "Fam_Gen_sp_GCA_000002.1",
        "Fam_Gen_sp_GCA_000003.1", "Fam_Gen_sp_GCA_000004.1",
    }
    assert any(p == "root" for _, p in edges)
    # genomes_map: each accession -> its leaf
    lines = gmap.read_text().splitlines()
    accs = {ln.split("\t")[0] for ln in lines}
    assert accs == {"GCA_000001.1", "GCA_000002.1", "GCA_000003.1", "GCA_000004.1"}


def test_tree2tax_outgroup_rooting(workdir: Path) -> None:
    _setup(workdir, with_outgroup=True)
    ctx = WorkdirContext(workdir, create=True)
    t2t, _ = tree2tax_run(ctx, Tree2taxParams(remove_outgroup=True))
    edges = _edges(t2t)
    # outgroup removed from the taxonomy after rooting
    assert not any(c == "Fam_Gen_sp_GCA_000004.1" for c, _ in edges)
    assert any(p == "root" for _, p in edges)


def test_tree2tax_include_dereplicated(workdir: Path) -> None:
    # a redundant genome listed under its representative leaf appears in genomes_map
    clusters = {"Fam_Gen_sp_GCA_000001.1": ["Fam_Gen_sp_GCA_000009.1"]}
    _setup(workdir, clusters=clusters)
    ctx = WorkdirContext(workdir, create=True)
    _, gmap = tree2tax_run(ctx, Tree2taxParams(include_dereplicated=True))
    accs = {ln.split("\t")[0] for ln in gmap.read_text().splitlines()}
    assert "GCA_000009.1" in accs  # the redundant member is mapped too


# An unrooted tree as mashtree/fasttree emit it: the seed node is trifurcating and
# the outgroup sits deep inside, next to one ingroup leaf.
_UNROOTED_NWK = (
    "(Fam_Gen_sp_GCA_000001.1:0.01,Fam_Gen_sp_GCA_000002.1:0.01,"
    "((Fam_Gen_sp_GCA_000003.1:0.02,Fam_Out_sp_GCA_000009.1:0.2):0.01,"
    "Fam_Gen_sp_GCA_000004.1:0.01):0.01);"
)


def test_tree2tax_roots_on_outgroup_edge_for_unrooted_tree(workdir: Path) -> None:
    (workdir / "tree").mkdir(parents=True)
    (workdir / "tree" / "tree.nwk").write_text(_UNROOTED_NWK + "\n")
    og = workdir / "outgroup"
    og.mkdir(parents=True)
    (og / "Fam_Out_sp_GCA_000009.1.fasta").write_text(">x\nACGT\n")
    (workdir / "outgroup_accession.txt").write_text("GCA_000009.1\n")
    ctx = WorkdirContext(workdir, create=True)
    t2t, _ = tree2tax_run(ctx, Tree2taxParams())
    edges = _edges(t2t)
    root_children = {c for c, p in edges if p == "root"}
    # The root splits the outgroup from one ingroup clade; no ingroup leaf may be
    # placed beside the outgroup because the input tree happened to be unrooted.
    assert "Fam_Out_sp_GCA_000009.1" in root_children
    assert len(root_children) == 2
    (ingroup,) = root_children - {"Fam_Out_sp_GCA_000009.1"}
    # The outgroup's former neighbour is inside the ingroup clade, not at the root.
    assert not ingroup.startswith("Fam_")
    assert ("Fam_Gen_sp_GCA_000003.1", "root") not in edges
