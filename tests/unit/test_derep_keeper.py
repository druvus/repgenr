"""Quality-aware keeper: the best-scoring cluster member becomes the representative."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    DerepResult,
    check_result_complete,
)
from repgenr.stages.derep_keeper import quality_score, rescore_representatives

_LOG = logging.getLogger("keeper")


def _result() -> DerepResult:
    return DerepResult(
        representatives=[Path("/g/rep.fasta"), Path("/g/solo.fasta")],
        clusters={"rep.fasta": ["m1.fasta", "m2.fasta"], "solo.fasta": []},
        genome_status={
            "rep.fasta": STATUS_REPRESENTATIVE, "solo.fasta": STATUS_REPRESENTATIVE,
            "m1.fasta": STATUS_CONTAINED, "m2.fasta": STATUS_CONTAINED,
        },
    )


def test_score_penalises_contamination() -> None:
    assert quality_score(100.0, 0.0) > quality_score(100.0, 2.0)
    assert quality_score(95.0, 0.0) == 95.0


def test_better_member_replaces_representative() -> None:
    quality = {"rep.fasta": (90.0, 3.0), "m1.fasta": (99.0, 0.2), "m2.fasta": (95.0, 1.0)}
    out, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 1
    assert sorted(p.name for p in out.representatives) == ["m1.fasta", "solo.fasta"]
    assert out.clusters["m1.fasta"] == ["m2.fasta", "rep.fasta"]
    assert out.genome_status["m1.fasta"] == STATUS_REPRESENTATIVE
    assert out.genome_status["rep.fasta"] == STATUS_CONTAINED
    assert out.representatives[0].parent == Path("/g")
    check_result_complete(out, ["rep.fasta", "solo.fasta", "m1.fasta", "m2.fasta"])


def test_scored_member_wins_over_unscored_representative() -> None:
    """An unscored representative has no quality to defend; any scored member
    -- even a mediocre one -- must replace it."""
    quality = {"m1.fasta": (60.0, 2.0)}  # rep.fasta and m2.fasta carry no quality
    out, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 1
    assert sorted(p.name for p in out.representatives) == ["m1.fasta", "solo.fasta"]
    assert out.clusters["m1.fasta"] == ["m2.fasta", "rep.fasta"]
    assert out.genome_status["m1.fasta"] == STATUS_REPRESENTATIVE
    assert out.genome_status["rep.fasta"] == STATUS_CONTAINED


def test_unscored_member_never_wins() -> None:
    quality = {"rep.fasta": (90.0, 3.0)}
    out, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 0
    assert out.clusters == _result().clusters


def test_no_quality_keeps_adapter_choice() -> None:
    out, swaps = rescore_representatives(_result(), {}, _LOG)
    assert swaps == 0
    assert out == _result()


def test_tie_keeps_current_representative() -> None:
    quality = {"rep.fasta": (99.0, 0.0), "m1.fasta": (99.0, 0.0), "m2.fasta": (50.0, 0.0)}
    _, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 0


def test_fail_qc_genome_survives_rescore_untouched() -> None:
    """A genome with no cluster (rejected on QC upstream) keeps its status and
    still lets the rescored result pass check_result_complete."""
    result = _result()
    result.genome_status["bad.fasta"] = STATUS_FAIL_QC
    quality = {"rep.fasta": (90.0, 3.0), "m1.fasta": (99.0, 0.2), "m2.fasta": (95.0, 1.0)}

    out, swaps = rescore_representatives(result, quality, _LOG)

    assert swaps == 1
    assert out.genome_status["bad.fasta"] == STATUS_FAIL_QC
    check_result_complete(out, ["rep.fasta", "solo.fasta", "m1.fasta", "m2.fasta", "bad.fasta"])
