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
