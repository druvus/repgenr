"""Outgroup resolution matches the accession exactly; stale files are pruned."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.stages.phylo import resolve_outgroup_files
from repgenr.stages.tree2tax import _resolve_outgroup_leaf_from

_LOG = logging.getLogger("test")


def _setup(tmp_path: Path, accession: str, files: list[str]) -> tuple[Path, Path]:
    og_dir = tmp_path / "outgroup"
    og_dir.mkdir()
    for name in files:
        (og_dir / name).write_text(">og\nACGT\n", encoding="utf-8")
    acc_file = tmp_path / "outgroup_accession.txt"
    acc_file.write_text(accession + "\n", encoding="utf-8")
    return og_dir, acc_file


def test_exact_accession_match_wins_over_substring(tmp_path: Path) -> None:
    """A stale outgroup file whose name merely CONTAINS the accession must not
    shadow the file whose parsed accession IS the accession."""
    og_dir, acc_file = _setup(
        tmp_path, "GCF_1.1",
        # the stale file sorts first AND contains "GCF_1.1" as a substring,
        # but its parsed accession is "1.1_backup", not "GCF_1.1"
        ["AAA_Old_GCF_1.1_backup.fasta", "Fam_Out_grp_GCF_1.1.fasta"],
    )
    file, leaf = resolve_outgroup_files(og_dir, acc_file, _LOG)
    assert file is not None and file.name == "Fam_Out_grp_GCF_1.1.fasta"
    assert leaf == "Fam_Out_grp_GCF_1.1"

    assert _resolve_outgroup_leaf_from(og_dir, acc_file, _LOG) == "Fam_Out_grp_GCF_1.1"


def test_substring_fallback_still_resolves_with_warning(tmp_path: Path, caplog) -> None:
    """Non-canonical outgroup filenames (viral: '<accession>.fasta' variants)
    keep resolving through the substring fallback."""
    og_dir, acc_file = _setup(tmp_path, "MN908947.3", ["outgroup_MN908947.3.fa"])
    with caplog.at_level(logging.INFO):
        file, leaf = resolve_outgroup_files(og_dir, acc_file, _LOG)
    assert file is not None and file.name == "outgroup_MN908947.3.fa"


def test_unresolvable_outgroup_warns_and_returns_none(tmp_path: Path) -> None:
    og_dir, acc_file = _setup(tmp_path, "GCF_9.9", ["Fam_Out_grp_GCF_1.1.fasta"])
    file, leaf = resolve_outgroup_files(og_dir, acc_file, _LOG)
    assert file is None and leaf is None
