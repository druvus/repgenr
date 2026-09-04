"""Input-completeness guards: a partial upstream output is refused, not processed."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.contracts import (
    MISSING_ACCESSIONS_TXT,
    SelectionRow,
    write_clusters,
    write_selection,
)
from repgenr.core.errors import WorkdirError
from repgenr.core.integrity import (
    check_genome_completeness,
    check_representatives_consistency,
    looks_like_fasta,
)

_LOG = logging.getLogger("test")


def _selection(workdir: Path, names: list[str], outgroup: str | None = None) -> None:
    rows = [
        SelectionRow(f"ACC_{i}", "Fam", "Gen", f"sp{i}", False, name)
        for i, name in enumerate(names)
    ]
    if outgroup:
        rows.append(SelectionRow("ACC_OG", "Fam", "Out", "grp", True, outgroup))
    write_selection(workdir / "selection.tsv", rows)


def _fill(genomes_dir: Path, names: list[str]) -> None:
    genomes_dir.mkdir(parents=True, exist_ok=True)
    for name in names:
        (genomes_dir / name).write_text(">x\nACGT\n", encoding="utf-8")


def test_complete_set_passes(tmp_path: Path) -> None:
    _selection(tmp_path, ["a.fasta", "b.fasta"], outgroup="og.fasta")
    _fill(tmp_path / "genomes", ["a.fasta", "b.fasta"])
    check_genome_completeness(tmp_path / "genomes", tmp_path, logger=_LOG)


def test_partial_set_refused_naming_shortfall(tmp_path: Path) -> None:
    _selection(tmp_path, ["a.fasta", "b.fasta", "c.fasta"])
    _fill(tmp_path / "genomes", ["a.fasta"])
    with pytest.raises(WorkdirError, match="b.fasta"):
        check_genome_completeness(tmp_path / "genomes", tmp_path, logger=_LOG)


def test_allow_incomplete_downgrades_to_warning(tmp_path: Path, caplog) -> None:
    _selection(tmp_path, ["a.fasta", "b.fasta"])
    _fill(tmp_path / "genomes", ["a.fasta"])
    with caplog.at_level(logging.WARNING):
        check_genome_completeness(
            tmp_path / "genomes", tmp_path, logger=_LOG, allow_incomplete=True
        )
    assert any("b.fasta" in r.getMessage() for r in caplog.records)


def test_known_missing_accessions_are_excused(tmp_path: Path) -> None:
    """NCBI legitimately returns nothing for some accessions; the genome stage
    records them and the completeness check must not fail on them."""
    _selection(tmp_path, ["a.fasta", "b.fasta"])
    _fill(tmp_path / "genomes", ["a.fasta"])
    (tmp_path / MISSING_ACCESSIONS_TXT).write_text("ACC_1\n", encoding="utf-8")
    check_genome_completeness(tmp_path / "genomes", tmp_path, logger=_LOG)


def test_no_selection_file_is_a_noop(tmp_path: Path) -> None:
    _fill(tmp_path / "genomes", ["whatever.fasta"])
    check_genome_completeness(tmp_path / "genomes", tmp_path, logger=_LOG)


def test_representatives_must_match_clusters(tmp_path: Path) -> None:
    reps = tmp_path / "derep" / "representatives"
    _fill(reps, ["r1.fasta"])  # crashed _write_contract: prefix of reps
    write_clusters(tmp_path / "derep" / "clusters.tsv", {"r1.fasta": [], "r2.fasta": ["m.fasta"]})
    with pytest.raises(WorkdirError, match="r2.fasta"):
        check_representatives_consistency(reps, tmp_path / "derep" / "clusters.tsv", logger=_LOG)
    _fill(reps, ["r2.fasta"])
    check_representatives_consistency(reps, tmp_path / "derep" / "clusters.tsv", logger=_LOG)


def test_representatives_check_noop_without_clusters(tmp_path: Path) -> None:
    reps = tmp_path / "derep" / "representatives"
    _fill(reps, ["r1.fasta"])
    check_representatives_consistency(reps, tmp_path / "derep" / "clusters.tsv", logger=_LOG)


def test_looks_like_fasta(tmp_path: Path) -> None:
    good = tmp_path / "g.fasta"
    good.write_text(">seq\nACGT\n", encoding="utf-8")
    bad = tmp_path / "b.fasta"
    bad.write_text("<html>Service unavailable</html>", encoding="utf-8")
    assert looks_like_fasta(good)
    assert not looks_like_fasta(bad)
    assert not looks_like_fasta(tmp_path / "absent.fasta")
