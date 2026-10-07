"""Grouped --help: every command sits in one panel, in pipeline order."""

from __future__ import annotations

import re

import typer
from typer.testing import CliRunner

from repgenr.cli.base import COMMAND_ORDER, COMMAND_PANELS
from repgenr.cli.main import app

CLICK = typer.main.get_command(app)


def _help() -> str:
    out = CliRunner().invoke(app, ["--help"], terminal_width=120).output
    return re.sub(r"\x1b\[[0-9;]*m", "", out)


def test_every_command_is_in_exactly_one_panel() -> None:
    assigned = [c for cmds in COMMAND_PANELS.values() for c in cmds]
    assert sorted(assigned) == sorted(CLICK.commands)
    assert len(assigned) == len(set(assigned))


def test_registered_panel_matches_mapping() -> None:
    for panel, cmds in COMMAND_PANELS.items():
        for name in cmds:
            assert CLICK.commands[name].rich_help_panel == panel, name


def test_command_order_follows_mapping() -> None:
    assert CLICK.list_commands(None) == list(COMMAND_ORDER)


def test_help_shows_panels_in_order() -> None:
    out = _help()
    positions = [out.index(panel) for panel in COMMAND_PANELS]
    assert positions == sorted(positions)
    rows = [m.group(1) for m in re.finditer(r"^\u2502 (\S+) ", out, re.MULTILINE)]
    listed = [r for r in rows if r in COMMAND_ORDER]
    assert listed == list(COMMAND_ORDER)
    assert "\u2500 Commands " not in out
