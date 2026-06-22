"""The resume fingerprint ignores non-result params (threads/num_processes)."""

from __future__ import annotations

from dataclasses import dataclass

from repgenr.cli.base import _stage_fingerprint


@dataclass
class _P:
    secondary_ani: float = 0.99
    threads: int = 16
    num_processes: int = 1


def test_threads_do_not_change_fingerprint() -> None:
    base = _stage_fingerprint("dereplicate", _P(threads=16, num_processes=1))
    more = _stage_fingerprint("dereplicate", _P(threads=64, num_processes=8))
    assert base == more  # only scheduling differs -> same fingerprint, resume skips


def test_result_param_changes_fingerprint() -> None:
    a = _stage_fingerprint("dereplicate", _P(secondary_ani=0.99))
    b = _stage_fingerprint("dereplicate", _P(secondary_ani=0.95))
    assert a != b


def test_stage_name_changes_fingerprint() -> None:
    assert _stage_fingerprint("phylo", _P()) != _stage_fingerprint("dereplicate", _P())
