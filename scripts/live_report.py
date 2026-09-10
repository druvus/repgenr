"""Render a pytest junit XML from the live suite as a Markdown results table.

Usage: python scripts/live_report.py <junit.xml> [more.xml ...] [--out docs/verification.md]

Several files merge by test; a later file's row wins, so a rerun of one
module can update the rows of an earlier full run.

Without --out the table is printed. With --out, the block between the
markers ``<!-- live-results:start -->`` and ``<!-- live-results:end -->`` in
that file is replaced, so the surrounding prose is kept.
"""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from datetime import UTC, datetime
from pathlib import Path

START = "<!-- live-results:start -->"
END = "<!-- live-results:end -->"


def _status(case: ET.Element) -> str:
    for tag, label in (("failure", "failed"), ("error", "error"), ("skipped", "skipped")):
        if case.find(tag) is not None:
            return label
    return "passed"


def _rows(root: ET.Element) -> list[tuple[str, str, str, float]]:
    rows = []
    for case in root.iter("testcase"):
        module = case.get("classname", "").split(".")[-1]
        rows.append((module, case.get("name", ""), _status(case), float(case.get("time", 0))))
    return rows


def render(*junits: Path) -> str:
    """Rows from every file; a later file overrides an earlier row of the same test."""
    merged: dict[tuple[str, str], tuple[str, str, str, float]] = {}
    for junit in junits:
        for row in _rows(ET.parse(junit).getroot()):
            merged[(row[0], row[1])] = row
    rows = list(merged.values())
    junit = max(junits, key=lambda j: j.stat().st_mtime)
    counts = {
        s: sum(1 for r in rows if r[2] == s) for s in ("passed", "failed", "error", "skipped")
    }
    total = sum(r[3] for r in rows)
    stamp = datetime.fromtimestamp(junit.stat().st_mtime, UTC).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        f"Last run {stamp}: {counts['passed']} passed, {counts['failed']} failed, "
        f"{counts['error']} errors, {counts['skipped']} skipped, {total / 60:.0f} min in total.",
        "",
        "| module | test | result | seconds |",
        "|---|---|---|---|",
    ]
    for module, name, status, seconds in sorted(rows):
        lines.append(f"| {module} | {name} | {status} | {seconds:.0f} |")
    return "\n".join(lines) + "\n"


def splice(text: str, table: str) -> str:
    if START not in text or END not in text:
        raise SystemExit(f"markers {START} / {END} not found")
    head, _, rest = text.partition(START)
    _, _, tail = rest.partition(END)
    return f"{head}{START}\n{table}{END}{tail}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("junit", type=Path, nargs="+", help="junit files; later ones override")
    parser.add_argument("--out", type=Path, help="Markdown file with the result markers to update.")
    args = parser.parse_args(argv)
    table = render(*args.junit)
    if args.out is None:
        sys.stdout.write(table)
        return 0
    args.out.write_text(splice(args.out.read_text(encoding="utf-8"), table), encoding="utf-8")
    print(f"updated {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
