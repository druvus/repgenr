"""tree2tax stage test against a synthetic rooted tree (uses real ete3)."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.stages.tree2tax import Tree2taxParams, run


def _setup(workdir: Path) -> WorkdirContext:
    # tree leaves are genome stems; include an outgroup leaf
    tree_dir = workdir / "tree"
    tree_dir.mkdir(parents=True)
    (tree_dir / "tree.nwk").write_text(
        "((Fam_gen_sp_GCA_000001:0.1,Fam_gen_sp_GCA_000002:0.1):0.2,Out_gen_sp_GCA_000099:0.5);\n"
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


@pytest.mark.parametrize("unreachable", [False, True])
def test_missing_workdir_exits_3_without_creating_it(tmp_path: Path, unreachable: bool) -> None:
    """A nonexistent -wd is a workdir error (exit 3), not a traceback or a new directory."""
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    if unreachable:
        blocker = tmp_path / "a_file"
        blocker.write_text("")
        missing = blocker / "wd"  # mkdir under a regular file fails
    else:
        missing = tmp_path / "missing_wd"
    result = CliRunner().invoke(app, ["tree2tax", "-wd", str(missing)])
    assert result.exit_code == 3, result.output
    assert "Traceback" not in result.output
    assert not missing.exists()


def test_flag_change_is_not_reported_as_a_changed_input(workdir: Path) -> None:
    """--no-include-dereplicated drops derep/clusters.tsv from the stage inputs;
    the rerun must not claim that file changed when it did not."""
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    _setup(workdir).close()
    runner = CliRunner()
    first = runner.invoke(app, ["tree2tax", "-wd", str(workdir)])
    assert first.exit_code == 0, first.output
    for flag in ("--no-include-dereplicated", "--include-dereplicated"):
        result = runner.invoke(app, ["tree2tax", "-wd", str(workdir), flag])
        assert result.exit_code == 0, result.output
    log = (workdir / "repgenr.log").read_text()
    assert log.count("Wrote tree2tax.tsv") == 3  # each flag change re-ran
    assert "'derep/clusters.tsv' changed" not in log


def test_member_that_is_a_leaf_maps_only_to_its_own_leaf(workdir: Path) -> None:
    """A tree built with --all-genomes has every genome as a leaf; a contained
    genome must not be listed a second time under its representative."""
    ctx = _setup(workdir)
    (workdir / "tree" / "tree.nwk").write_text(
        "(((Fam_gen_sp_GCA_000001:0.1,Fam_gen_sp_GCA_000003:0.1):0.1,"
        "Fam_gen_sp_GCA_000002:0.1):0.2,Out_gen_sp_GCA_000099:0.5);\n"
    )
    _t2t, gmap = run(ctx, Tree2taxParams(include_dereplicated=True))
    rows = [tuple(line.split("\t")) for line in gmap.read_text().splitlines()]
    accessions = [acc for acc, _ in rows]
    assert len(accessions) == len(set(accessions)), rows
    assert ("GCA_000003", "Fam_gen_sp_GCA_000003") in rows
