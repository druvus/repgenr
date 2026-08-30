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


def test_docker_extra_mounts(tmp_path, monkeypatch) -> None:
    # A directory referenced indirectly (not in argv) is mounted when declared.
    # Isolate the temp dir so `genomes` is not nested under the default tempdir
    # mount (which would correctly dedup it away and mask what we're testing).
    sys_tmp = tmp_path / "systmp"
    sys_tmp.mkdir()
    monkeypatch.setattr(containers.tempfile, "gettempdir", lambda: str(sys_tmp))
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


@pytest.mark.parametrize("backend", ["none", "docker"])
def test_run_tool_forwards_timeout(monkeypatch, backend) -> None:
    configure_container(backend)
    captured = {}

    def fake_run(command, **kw):
        captured["timeout"] = kw.get("timeout")
        return 0

    monkeypatch.setattr(containers.process, "run", fake_run)
    caps = ToolCapabilities(name="skder", container="img:1")
    containers.run_tool(caps, ["skder", "-h"], logger=_LOG, timeout=42.0)
    assert captured["timeout"] == 42.0


# --- Wave image minting robustness --------------------------------------------


def _wave_caps() -> ToolCapabilities:
    return ToolCapabilities(name="sourmash", conda=("bioconda::sourmash",))


@pytest.fixture()
def _wave_env(monkeypatch):
    containers._WAVE_CACHE.clear()
    monkeypatch.setattr(containers.shutil, "which", lambda name: "/usr/bin/wave")
    yield
    containers._WAVE_CACHE.clear()


def test_wave_timeout_raises_tool_error(monkeypatch, _wave_env) -> None:
    import subprocess

    def fake_run(cmd, **kwargs):
        raise subprocess.TimeoutExpired(cmd, 600)

    monkeypatch.setattr(subprocess, "run", fake_run)
    cfg = ContainerConfig(backend="docker", wave_enabled=True)
    with pytest.raises(containers.ToolExecutionError, match="timed out"):
        resolve_image(_wave_caps(), cfg)


def test_wave_empty_stdout_raises_tool_error(monkeypatch, _wave_env) -> None:
    import subprocess
    from types import SimpleNamespace

    def fake_run(cmd, **kwargs):
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    cfg = ContainerConfig(backend="docker", wave_enabled=True)
    with pytest.raises(containers.ToolExecutionError, match="no image"):
        resolve_image(_wave_caps(), cfg)


# --- retrying tool runner -----------------------------------------------------


def test_run_tool_with_retries_recovers(monkeypatch) -> None:
    calls = {"n": 0}

    def flaky(caps, cmd, **kwargs):
        calls["n"] += 1
        if calls["n"] < 3:
            raise containers.ToolExecutionError(list(cmd), 1, output="net down")
        return 0

    monkeypatch.setattr(containers, "run_tool", flaky)
    monkeypatch.setattr(containers.time, "sleep", lambda s: None)
    rc = containers.run_tool_with_retries(
        _wave_caps(), ["datasets", "download"], logger=_LOG, attempts=3
    )
    assert rc == 0 and calls["n"] == 3


def test_run_tool_with_retries_exhausts(monkeypatch) -> None:
    def always_fail(caps, cmd, **kwargs):
        raise containers.ToolExecutionError(list(cmd), 1, output="net down")

    monkeypatch.setattr(containers, "run_tool", always_fail)
    monkeypatch.setattr(containers.time, "sleep", lambda s: None)
    with pytest.raises(containers.ToolExecutionError):
        containers.run_tool_with_retries(
            _wave_caps(), ["datasets", "download"], logger=_LOG, attempts=2
        )
