"""Every input genome must leave dereplication with a status and a home."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    DerepResult,
    check_result_complete,
)


def _ok() -> DerepResult:
    return DerepResult(
        representatives=[Path("a.fasta")],
        clusters={"a.fasta": ["b.fasta"]},
        genome_status={"a.fasta": STATUS_REPRESENTATIVE, "b.fasta": STATUS_CONTAINED},
    )


def test_complete_result_passes() -> None:
    check_result_complete(_ok(), ["a.fasta", "b.fasta"])


def test_missing_status_raises() -> None:
    with pytest.raises(WorkdirError, match="c.fasta"):
        check_result_complete(_ok(), ["a.fasta", "b.fasta", "c.fasta"])


def test_contained_without_cluster_raises() -> None:
    r = _ok()
    r.clusters = {"a.fasta": []}
    with pytest.raises(WorkdirError, match="no representative"):
        check_result_complete(r, ["a.fasta", "b.fasta"])


def test_representative_without_status_raises() -> None:
    r = _ok()
    r.genome_status = {"b.fasta": STATUS_CONTAINED}
    with pytest.raises(WorkdirError, match="a.fasta"):
        check_result_complete(r, ["a.fasta", "b.fasta"])


def test_unknown_status_value_raises() -> None:
    r = _ok()
    r.genome_status["b.fasta"] = "weird"
    with pytest.raises(WorkdirError, match="weird"):
        check_result_complete(r, ["a.fasta", "b.fasta"])


def test_fail_qc_needs_no_cluster() -> None:
    r = DerepResult(
        representatives=[Path("a.fasta")],
        clusters={"a.fasta": []},
        genome_status={"a.fasta": STATUS_REPRESENTATIVE, "b.fasta": "fail_qc"},
    )
    check_result_complete(r, ["a.fasta", "b.fasta"])
