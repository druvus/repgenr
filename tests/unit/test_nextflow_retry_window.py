"""The Nextflow retry window covers resource kills, not a user's interrupt."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _retry_exit_codes() -> set[int]:
    config = (ROOT / "nextflow" / "nextflow.config").read_text(encoding="utf-8")
    line = next(ln for ln in config.splitlines() if "errorStrategy" in ln and "exitStatus" in ln)
    codes: set[int] = set()
    for lo, hi in re.findall(r"\((\d+)\.\.(\d+)\)", line):
        codes |= set(range(int(lo), int(hi) + 1))
    for group in re.findall(r"\[([\d,\s]+)\]", line):
        codes |= {int(x) for x in group.split(",") if x.strip()}
    assert codes, f"could not parse the retry window from: {line.strip()}"
    return codes


def test_oom_and_signal_kills_are_retried() -> None:
    codes = _retry_exit_codes()
    assert {137, 139, 143} <= codes  # SIGKILL (OOM), SIGSEGV, SIGTERM (scheduler)


def test_user_interrupt_is_not_retried() -> None:
    # 130 = 128 + SIGINT: someone pressed Ctrl-C or cancelled the task on purpose.
    assert 130 not in _retry_exit_codes()
