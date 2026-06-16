"""Unit tests for the container execution backend (no daemon needed)."""

from __future__ import annotations

import logging

import pytest

from repgenr.core import containers
from repgenr.core.containers import (
    ContainerConfig,
    configure_container,
    resolve_image,
    wrap_command,
)
from repgenr.core.plugins import ToolCapabilities

_LOG = logging.getLogger("test")


@pytest.fixture(autouse=True)
def _reset_backend():
    yield
    configure_container("none")  # restore native default after each test


def test_resolve_prefers_explicit_image() -> None:
    caps = ToolCapabilities(name="cactus", container="quay.io/x/cactus:1", conda=("bioconda::x",))
    assert resolve_image(caps, ContainerConfig(backend="docker")) == "quay.io/x/cactus:1"


def test_resolve_none_without_wave() -> None:
    caps = ToolCapabilities(name="skder", conda=("bioconda::skder",))
    # docker backend but no explicit image and Wave disabled -> run native
    assert resolve_image(caps, ContainerConfig(backend="docker", wave_enabled=False)) is None


def test_docker_wrap_command() -> None:
    cfg = ContainerConfig(backend="docker", platform="linux/amd64")
    argv = ["skder", "-g", "/wd/a.fasta", "-o", "/wd/out"]
    cmd = wrap_command("img:1", argv, config=cfg, cwd="/wd", logger=_LOG)
    assert cmd[0] == "docker" and "run" in cmd and "--rm" in cmd
    assert "--platform" in cmd and "linux/amd64" in cmd
    assert "-w" in cmd and "/wd" in cmd
    # the image precedes the tool argv, which is preserved in order
    i = cmd.index("img:1")
    assert cmd[i + 1 :] == argv
    # workdir is bind-mounted
    assert any(c == "/wd:/wd" for c in cmd)
    # HOME points at the (mounted, writable) workdir, not the unwritable "/"
    assert any(c == "HOME=/wd" for c in cmd)


def test_docker_extra_mounts(tmp_path) -> None:
    # A directory referenced indirectly (not in argv) is mounted when declared.
    genomes = tmp_path / "reps"
    genomes.mkdir()
    cfg = ContainerConfig(backend="docker")
    cmd = wrap_command(
        "img:1", ["cactus", "seqfile.txt"], config=cfg, cwd="/wd",
        logger=_LOG, extra_mounts=[str(genomes)],
    )
    assert any(c == f"{genomes}:{genomes}" for c in cmd)


def test_singularity_wrap_command_no_cache() -> None:
    cfg = ContainerConfig(backend="singularity")  # no cache_dir -> docker:// ref
    argv = ["mashtree", "x.fasta"]
    cmd = wrap_command("img:1", argv, config=cfg, cwd="/wd", logger=_LOG)
    assert cmd[0] == "singularity" and cmd[1] == "exec"
    assert "--bind" in cmd and "--pwd" in cmd
    assert "docker://img:1" in cmd
    assert cmd[-len(argv):] == argv


def test_run_tool_native_when_backend_none(monkeypatch) -> None:
    configure_container("none")
    captured = {}

    def fake_run(command, **kw):
        captured["cmd"] = list(command)
        return 0

    monkeypatch.setattr(containers.process, "run", fake_run)
    caps = ToolCapabilities(name="skder", container="img:1")
    containers.run_tool(caps, ["skder", "-h"], logger=_LOG)
    # backend none -> not wrapped
    assert captured["cmd"] == ["skder", "-h"]


def test_run_tool_wraps_when_backend_active(monkeypatch) -> None:
    configure_container("docker")
    captured = {}

    def fake_run(command, **kw):
        captured["cmd"] = list(command)
        return 0

    monkeypatch.setattr(containers.process, "run", fake_run)
    caps = ToolCapabilities(name="skder", container="img:1")
    containers.run_tool(caps, ["skder", "-h"], logger=_LOG, cwd="/wd")
    assert captured["cmd"][0] == "docker"
    assert "img:1" in captured["cmd"]
    assert captured["cmd"][-2:] == ["skder", "-h"]
