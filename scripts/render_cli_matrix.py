"""Render tests/audit/cli_matrix.yaml to docs/audit/cli-matrix.md.

Usage: python scripts/render_cli_matrix.py   (from the repository root)
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "tests" / "audit" / "cli_matrix.yaml"
OUT = ROOT / "docs" / "audit" / "cli-matrix.md"
REFERENCE = ROOT / "docs" / "cli-reference.md"
SKIP_OPTS = {"--help", "--install-completion", "--show-completion"}


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return ", ".join(str(v) for v in value)
    return str(value).replace("|", "\\|")


def _table(flags: dict) -> list[str]:
    rows = [
        "| flag | aliases | param | validated | nextflow | live | docs |",
        "|---|---|---|---|---|---|---|",
    ]
    for flag, rec in flags.items():
        rows.append(
            f"| `{flag}` | {_cell(rec.get('aliases'))} | {_cell(rec['param'])} | "
            f"{_cell(rec['validated'])} | {_cell(rec.get('nextflow'))} | {_cell(rec['live'])} | "
            f"{_cell(rec['docs'])} |"
        )
    return rows


def render(matrix: dict) -> str:
    n_flags = len(matrix["global_flags"]) + sum(
        len(c["flags"]) for c in matrix["commands"].values()
    )
    todo = sum(
        1
        for c in matrix["commands"].values()
        for r in c["flags"].values()
        if str(r["live"]).startswith("todo:")
    ) + sum(1 for r in matrix["global_flags"].values() if str(r["live"]).startswith("todo:"))
    lines = [
        "# CLI matrix",
        "",
        "Generated from `tests/audit/cli_matrix.yaml` by `scripts/render_cli_matrix.py`;",
        "`tests/unit/test_cli_matrix.py` keeps both in step with the command tree.",
        "",
        f"{len(matrix['commands'])} commands, {n_flags} flags "
        f"({n_flags - todo} with a live test or an n/a reason, {todo} pending).",
        "",
        "## Global flags",
        "",
        *_table(matrix["global_flags"]),
        "",
    ]
    for name, cmd in matrix["commands"].items():
        lines += [f"## {name}", "", f"dispatch: `{cmd['dispatch']}`", ""]
        lines += _table(cmd["flags"]) if cmd["flags"] else ["(no flags)"]
        lines.append("")
    lines += [
        "## Short-alias collisions",
        "",
        *(
            f"- `{a}`: {', '.join(f'`{f}`' for f in fl)}"
            for a, fl in matrix["known_short_alias_collisions"].items()
        ),
        "",
    ]
    return "\n".join(lines)


def _click_tree():
    import typer

    from repgenr.cli.main import app

    return typer.main.get_command(app)


def _option_rows(cmd) -> list[str]:
    rows = ["| option | default | description |", "|---|---|---|"]
    for p in cmd.params:
        opts = getattr(p, "opts", None)
        if not opts or opts[0] in SKIP_OPTS:
            continue
        names = ", ".join(f"`{o}`" for o in [*opts, *(getattr(p, "secondary_opts", None) or [])])
        if getattr(p, "required", False):
            default = "required"
        elif getattr(p, "is_flag", False):
            default = "off" if not p.default else "on"
        else:
            default = "" if p.default in (None, [], ()) else f"`{p.default}`"
        rows.append(f"| {names} | {default} | {_cell(p.help or '')} |")
    return rows


def render_reference(cli=None) -> str:
    """The command reference: every command with its options, defaults and help."""
    cli = cli or _click_tree()
    lines = [
        "# Command reference",
        "",
        "Generated from the command tree by `scripts/render_cli_matrix.py`; the",
        "matrix test keeps it in sync. Global options go before the command name",
        "(`repgenr --container docker dereplicate ...`).",
        "",
        "## Global options",
        "",
        *_option_rows(cli),
        "",
    ]
    for name in sorted(cli.commands):
        cmd = cli.commands[name]
        lines += [f"## {name}", "", (cmd.help or "").strip(), ""]
        lines += (
            _option_rows(cmd)
            if any(getattr(p, "opts", None) for p in cmd.params)
            else ["(no options)"]
        )
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    with open(MATRIX, encoding="utf-8") as fh:
        matrix = yaml.safe_load(fh)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(render(matrix), encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    REFERENCE.write_text(render_reference(), encoding="utf-8")
    print(f"wrote {REFERENCE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
