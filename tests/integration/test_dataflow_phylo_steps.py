"""Stateless data-channel steps for phylo/tree2tax (no shared workdir).

Covers the pure, binary-free parts: the tree2tax-relations step end to end
(dendropy only) and the phylo input resolvers (genome listing, reference and
outgroup resolution) that let build_tree run without a WorkdirContext.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.stages.phylo import _resolve_reference, list_fasta, resolve_outgroup_files
from repgenr.stages.tree2tax import Tree2taxStepParams, tree2tax_relations

_LOG = logging.getLogger("test")

_NWK = (
    "((Fam_Gen_sp_GCA_000001.1,Fam_Gen_sp_GCA_000002.1),"
    "(Fam_Gen_sp_GCA_000003.1,Fam_Gen_sp_GCA_000004.1));"
)


def _edges(path: Path) -> set[tuple[str, str]]:
    return {tuple(line.split("\t")) for line in path.read_text().splitlines()[1:]}  # type: ignore[misc]


def test_tree2tax_relations_step(tmp_path: Path) -> None:
    tree = tmp_path / "tree.nwk"
    tree.write_text(_NWK + "\n")
    out = tmp_path / "out"
    t2t, gmap = tree2tax_relations(
        Tree2taxStepParams(tree=tree, out_dir=out), _LOG
    )
    assert t2t == out / "tree2tax.tsv"
    edges = _edges(t2t)
    leaves = {c for c, _ in edges if c.startswith("Fam_")}
    assert leaves == {
        "Fam_Gen_sp_GCA_000001.1", "Fam_Gen_sp_GCA_000002.1",
        "Fam_Gen_sp_GCA_000003.1", "Fam_Gen_sp_GCA_000004.1",
    }
    assert any(p == "root" for _, p in edges)
    accs = {ln.split("\t")[0] for ln in gmap.read_text().splitlines()}
    assert accs == {"GCA_000001.1", "GCA_000002.1", "GCA_000003.1", "GCA_000004.1"}


def test_tree2tax_relations_step_outgroup_and_dereplicated(tmp_path: Path) -> None:
    tree = tmp_path / "tree.nwk"
    tree.write_text(_NWK + "\n")
    og = tmp_path / "outgroup"
    og.mkdir()
    (og / "Fam_Gen_sp_GCA_000004.1.fasta").write_text(">x\nACGT\n")
    acc = tmp_path / "outgroup_accession.txt"
    acc.write_text("GCA_000004.1\n")

    from repgenr.core.contracts import write_clusters

    clusters = tmp_path / "clusters.tsv"
    write_clusters(clusters, {"Fam_Gen_sp_GCA_000001.1": ["Fam_Gen_sp_GCA_000009.1"]})

    out = tmp_path / "out"
    t2t, gmap = tree2tax_relations(
        Tree2taxStepParams(
            tree=tree, out_dir=out, clusters=clusters,
            outgroup_dir=og, outgroup_accession=acc,
            remove_outgroup=True, include_dereplicated=True,
        ),
        _LOG,
    )
    edges = _edges(t2t)
    assert not any(c == "Fam_Gen_sp_GCA_000004.1" for c, _ in edges)  # outgroup dropped
    accs = {ln.split("\t")[0] for ln in gmap.read_text().splitlines()}
    assert "GCA_000009.1" in accs  # redundant member mapped under its representative


def test_list_fasta(tmp_path: Path) -> None:
    (tmp_path / "a.fasta").write_text(">a\nAC\n")
    (tmp_path / "b.fna").write_text(">b\nGT\n")
    (tmp_path / ".hidden.fasta").write_text(">h\nNN\n")
    (tmp_path / "notes.txt").write_text("ignore")
    names = [p.name for p in list_fasta(tmp_path)]
    assert names == ["a.fasta", "b.fna"]
    assert list_fasta(tmp_path / "missing") == []


def test_resolve_reference(tmp_path: Path) -> None:
    g1 = tmp_path / "Fam_Gen_sp_GCA_1.fasta"
    g2 = tmp_path / "Fam_Gen_sp_GCA_2.fasta"
    for g in (g1, g2):
        g.write_text(">x\nAC\n")
    assert _resolve_reference(None, [g1, g2], None) == g1  # default: first
    assert _resolve_reference("Fam_Gen_sp_GCA_2.fasta", [g1, g2], None) == g2
    with pytest.raises(UserInputError):
        _resolve_reference("nope.fasta", [g1, g2], None)


def test_resolve_outgroup_files(tmp_path: Path) -> None:
    og = tmp_path / "outgroup"
    og.mkdir()
    f = og / "Fam_Gen_sp_GCA_000004.1.fasta"
    f.write_text(">x\nAC\n")
    acc = tmp_path / "outgroup_accession.txt"
    acc.write_text("GCA_000004.1\n")
    found, leaf = resolve_outgroup_files(og, acc, _LOG)
    assert found == f
    assert leaf == "Fam_Gen_sp_GCA_000004.1"
    # missing dir / accession -> no outgroup, no crash
    assert resolve_outgroup_files(tmp_path / "nope", acc, _LOG) == (None, None)
