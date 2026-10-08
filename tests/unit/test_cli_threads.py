"""-t/--threads defaults to the CPU limit of the process when that is lower."""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from repgenr.cli import cmd_misc
from repgenr.cli.main import app
from repgenr.core import resources

_runner = CliRunner()


@pytest.fixture()
def glance(monkeypatch, tmp_path):
    """Run `glance` with the stage replaced; return its threads and the output."""
    calls: list = []
    monkeypatch.setattr(
        cmd_misc, "_run", lambda stage, workdir, build, *, create=False: calls.append(build())
    )

    def invoke(*extra: str, cpus: int, quiet: bool = False):
        monkeypatch.setattr(resources, "usable_cpus", lambda: cpus)
        args = ["-q"] if quiet else []
        result = _runner.invoke(app, [*args, "glance", "-wd", str(tmp_path), *extra])
        assert result.exit_code == 0, result.output
        return calls.pop().threads, result

    return invoke


def test_the_default_is_lowered_to_the_cpu_limit(glance) -> None:
    threads, result = glance(cpus=4)
    assert threads == 4
    # The note goes to stderr; stdout stays free for a command's own output.
    assert "Using 4 threads, the CPU limit of this process" in result.stderr
    assert result.stdout == ""


def test_the_default_stays_16_with_enough_cpus(glance) -> None:
    threads, result = glance(cpus=64)
    assert threads == 16
    assert "Using" not in result.output and "WARNING -t/--threads" not in result.output


def test_quiet_hides_the_note(glance) -> None:
    threads, result = glance(cpus=4, quiet=True)
    assert threads == 4 and "Using 4 threads" not in result.output


def test_an_explicit_count_above_the_limit_is_kept_with_a_warning(glance) -> None:
    threads, result = glance("-t", "32", cpus=4)
    assert threads == 32
    assert result.stderr.count("WARNING -t/--threads 32 exceeds the 4 CPU(s)") == 1
    assert result.stdout == ""


def test_an_explicit_count_within_the_limit_is_silent(glance) -> None:
    threads, result = glance("-t", "4", cpus=4)
    assert threads == 4 and "WARNING -t/--threads" not in result.output


def _commands(typer_app):  # noqa: ANN001, ANN202
    for info in typer_app.registered_commands:
        yield info.callback
    for group in typer_app.registered_groups:
        yield from _commands(group.typer_instance)


def test_every_threads_option_uses_the_shared_callback() -> None:
    import inspect

    from repgenr.cli.base import resolve_threads

    found = []
    for function in _commands(app):
        for parameter in inspect.signature(function).parameters.values():
            default = parameter.default
            if "--threads" in (getattr(default, "param_decls", None) or ()):
                found.append(function.__name__)
                assert default.callback is resolve_threads, function.__name__
    assert len(found) >= 11, found
