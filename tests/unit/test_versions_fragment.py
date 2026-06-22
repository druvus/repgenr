"""The versions.yml fragment writer used by the data-channel steps."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.core.versions import write_versions_fragment
from repgenr.stages.tree2tax import Tree2taxStepParams, tree2tax_relations

_NWK = "(Fam_Gen_sp_GCA_1.1,Fam_Gen_sp_GCA_2.1);"


def test_fragment_is_sorted_and_indented(tmp_path: Path) -> None:
    out = tmp_path / "v.yml"
    write_versions_fragment(out, {"skder": "1.3.0", "skani": "0.2.1"})
    # 4-space indent (slots under a process key) and sorted by tool name
    assert out.read_text() == "    skani: 0.2.1\n    skder: 1.3.0\n"


def test_fragment_empty_is_empty_file(tmp_path: Path) -> None:
    out = tmp_path / "v.yml"
    write_versions_fragment(out, {})
    assert out.read_text() == ""


def test_tree2tax_step_writes_empty_versions(tmp_path: Path) -> None:
    tree = tmp_path / "tree.nwk"
    tree.write_text(_NWK + "\n")
    versions = tmp_path / "tool_versions.yml"
    tree2tax_relations(
        Tree2taxStepParams(tree=tree, out_dir=tmp_path / "out", versions_out=versions),
        logging.getLogger("test"),
    )
    # tree2tax uses no external tool -> empty fragment (module still adds repgenr)
    assert versions.exists() and versions.read_text() == ""
