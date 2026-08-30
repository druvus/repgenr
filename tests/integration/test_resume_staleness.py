"""The resume guard warns when a skipped stage is older than an upstream re-run."""

from __future__ import annotations

import logging

from repgenr.cli import base as cli
from repgenr.core.config import Config


def test_warn_if_stale_detects_newer_upstream(caplog) -> None:
    cfg = Config()
    # 'phylo' completed before 'dereplicate' was re-run -> phylo is stale
    cfg.record_stage("metadata", completed="2026-01-01T00:00:00")
    cfg.record_stage("phylo", completed="2026-01-01T01:00:00")
    cfg.record_stage("dereplicate", completed="2026-01-01T02:00:00")

    logger = logging.getLogger("repgenr.test.stale")
    with caplog.at_level(logging.WARNING, logger="repgenr.test.stale"):
        cli._warn_if_stale("phylo", cfg, logger)
    assert any("may be stale" in r.message for r in caplog.records)


def test_warn_if_stale_quiet_when_current(caplog) -> None:
    cfg = Config()
    cfg.record_stage("dereplicate", completed="2026-01-01T00:00:00")
    cfg.record_stage("phylo", completed="2026-01-01T01:00:00")  # newer than upstream: fine

    logger = logging.getLogger("repgenr.test.stale2")
    with caplog.at_level(logging.WARNING, logger="repgenr.test.stale2"):
        cli._warn_if_stale("phylo", cfg, logger)
    assert not any("may be stale" in r.message for r in caplog.records)


def test_warn_if_stale_compares_parsed_datetimes(caplog) -> None:
    """Mixed ISO offset forms compare by instant, not lexicographically."""
    cfg = Config()
    # +00:00 sorts after Z lexicographically although it is the same instant
    # format family; the naive form here is strictly older by instant.
    cfg.record_stage("dereplicate", completed="2026-01-01T02:00:00+00:00")
    cfg.record_stage("phylo", completed="2026-01-01T03:00:00+00:00")

    logger = logging.getLogger("repgenr.test.stale3")
    with caplog.at_level(logging.WARNING, logger="repgenr.test.stale3"):
        cli._warn_if_stale("phylo", cfg, logger)
    assert not any("may be stale" in r.message for r in caplog.records)


def test_warn_if_stale_ignores_malformed_timestamp(caplog) -> None:
    cfg = Config()
    cfg.record_stage("dereplicate", completed="not-a-timestamp")
    cfg.record_stage("phylo", completed="2026-01-01T01:00:00")

    logger = logging.getLogger("repgenr.test.stale4")
    with caplog.at_level(logging.WARNING, logger="repgenr.test.stale4"):
        cli._warn_if_stale("phylo", cfg, logger)  # must not raise
    assert not any("may be stale" in r.message for r in caplog.records)
