"""scripts/live_report.py turns a junit XML into the verification results table."""

from __future__ import annotations

from pathlib import Path

from scripts.live_report import END, START, render, splice

_JUNIT = """<?xml version="1.0" encoding="utf-8"?>
<testsuites><testsuite name="pytest" tests="3">
<testcase classname="tests.live.test_smoke" name="test_offline_chain" time="4.9"/>
<testcase classname="tests.live.test_network" name="test_api_genus" time="12.0">
<failure message="x"/></testcase>
<testcase classname="tests.live.test_container_runs" name="test_cactus" time="600">
<skipped message="docker"/></testcase>
</testsuite></testsuites>
"""


def test_render_counts_and_rows(tmp_path: Path) -> None:
    junit = tmp_path / "junit.xml"
    junit.write_text(_JUNIT, encoding="utf-8")
    table = render(junit)
    assert "1 passed, 1 failed, 0 errors, 1 skipped, 10 min in total." in table
    assert "| test_smoke | test_offline_chain | passed | 5 |" in table
    assert "| test_network | test_api_genus | failed | 12 |" in table
    assert "| test_container_runs | test_cactus | skipped | 600 |" in table


def test_splice_replaces_only_the_marked_block() -> None:
    doc = f"# Title\n\nintro\n\n{START}\nold\n{END}\n\ntail\n"
    out = splice(doc, "new table\n")
    assert out == f"# Title\n\nintro\n\n{START}\nnew table\n{END}\n\ntail\n"


def test_later_file_overrides_earlier_rows(tmp_path: Path) -> None:
    first = tmp_path / "a.xml"
    first.write_text(_JUNIT, encoding="utf-8")
    second = tmp_path / "b.xml"
    second.write_text(
        '<testsuites><testsuite name="pytest" tests="1">'
        '<testcase classname="tests.live.test_network" name="test_api_genus" time="3"/>'
        "</testsuite></testsuites>",
        encoding="utf-8",
    )
    table = render(first, second)
    assert "| test_network | test_api_genus | passed | 3 |" in table
    assert "2 passed, 0 failed, 0 errors, 1 skipped" in table
