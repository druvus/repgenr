"""Global container options that have no effect are reported, not ignored."""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.containers import configure_container

_runner = CliRunner()


@pytest.fixture(autouse=True)
def _native_after():
    yield
    configure_container("none")


@pytest.mark.parametrize(
    ("args", "flag"),
    [
        (["--platform", "linux/amd64"], "--platform"),
        (["--wave"], "--wave"),
        (["--container-engine", "podman"], "--container-engine"),
        (["--container-cache", "/tmp/sif"], "--container-cache"),
    ],
)
def test_container_options_without_a_backend_are_reported(args, flag) -> None:
    result = _runner.invoke(app, [*args, "list-tools"])
    assert result.exit_code == 0
    assert f"{flag} has no effect without --container" in result.stderr


def test_container_cache_with_docker_is_reported() -> None:
    result = _runner.invoke(
        app, ["--container", "docker", "--container-cache", "/tmp/x", "list-tools"]
    )
    assert result.exit_code == 0
    assert "--container-cache is used only by --container singularity" in result.stderr


def test_effective_options_are_not_reported() -> None:
    result = _runner.invoke(
        app, ["--container", "docker", "--platform", "linux/amd64", "list-tools"]
    )
    assert result.exit_code == 0
    assert "no effect" not in result.stderr


@pytest.mark.parametrize("via_env", [False, True])
def test_an_unknown_backend_is_a_usage_error(monkeypatch, via_env) -> None:
    # It ended in a traceback with exit 1.
    if via_env:
        monkeypatch.setenv("REPGENR_CONTAINER", "bogus")
        result = _runner.invoke(app, ["list-tools"])
    else:
        result = _runner.invoke(app, ["--container", "bogus", "list-tools"])
    assert result.exit_code == 2
    assert "Unknown container backend 'bogus'" in result.output
    assert "Traceback" not in result.output
    assert result.exception is None or isinstance(result.exception, SystemExit)
