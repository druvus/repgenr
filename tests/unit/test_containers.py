"""Unit tests for the container execution backend (no daemon needed)."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pytest

from repgenr.core import containers
from repgenr.core.containers import (
    ContainerConfig,
    configure_container,
    resolve_image,
    wrap_command,
)
from repgenr.core.errors import MissingBinaryError
from repgenr.core.plugins import ToolCapabilities

_LOG = logging.getLogger("test")


@pytest.fixture(autouse=True)
def _reset_backend():
    yield
    configure_container("none")  # restore native default after each test


def test_resolve_prefers_explicit_image() -> None:
    caps = ToolCapabilities(name="cactus", container="quay.io/x/cactus:1", conda=("bioconda::x",))
    assert resolve_image(caps, ContainerConfig(backend="docker")) == "quay.io/x/cactus:1"


def test_resolve_wave_wins_over_the_pin_when_enabled(monkeypatch) -> None:
    # --wave asks for an image minted from the conda spec (native architecture);
    # the pinned BioContainer is the default when Wave is off.
    caps = ToolCapabilities(name="skder", container="quay.io/x/skder:1", conda=("bioconda::skder",))
    monkeypatch.setattr(containers, "_wave_image", lambda spec, cfg: "wave.seqera.io/x/skder")
    assert (
        resolve_image(caps, ContainerConfig(backend="docker", wave_enabled=True))
        == "wave.seqera.io/x/skder"
    )
    assert (
        resolve_image(caps, ContainerConfig(backend="docker", wave_enabled=False))
        == "quay.io/x/skder:1"
    )


def test_resolve_pin_when_wave_has_no_conda_spec(monkeypatch) -> None:
    caps = ToolCapabilities(name="cactus", container="quay.io/x/cactus:1")
    monkeypatch.setattr(containers, "_wave_image", lambda spec, cfg: "must-not-be-called")
    assert (
        resolve_image(caps, ContainerConfig(backend="docker", wave_enabled=True))
        == "quay.io/x/cactus:1"
    )


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


def test_docker_run_uses_an_init_process() -> None:
    """With --init the tool is not process 1 in the container, so the SIGTERM
    the docker client forwards when repgenr stops a tool ends the container."""
    cfg = ContainerConfig(backend="docker")
    cmd = wrap_command("img:1", ["fasttree", "x"], config=cfg, cwd="/wd", logger=_LOG)
    assert "--init" in cmd
    assert cmd.index("--init") < cmd.index("img:1")


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
        "img:1",
        ["cactus", "seqfile.txt"],
        config=cfg,
        cwd="/wd",
        logger=_LOG,
        extra_mounts=[str(genomes)],
    )
    assert any(c == f"{genomes}:{genomes}" for c in cmd)


def test_singularity_wrap_command_no_cache() -> None:
    cfg = ContainerConfig(backend="singularity")  # no cache_dir -> docker:// ref
    argv = ["mashtree", "x.fasta"]
    cmd = wrap_command("img:1", argv, config=cfg, cwd="/wd", logger=_LOG)
    assert cmd[0] == "singularity" and cmd[1] == "exec"
    assert "--bind" in cmd and "--pwd" in cmd
    assert "docker://img:1" in cmd
    assert cmd[-len(argv) :] == argv


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
    with pytest.raises(containers.ToolExecutionError) as ei:
        resolve_image(_wave_caps(), cfg)
    assert "timed out" in ei.value.details()


def test_wave_empty_stdout_raises_tool_error(monkeypatch, _wave_env) -> None:
    import subprocess
    from types import SimpleNamespace

    def fake_run(cmd, **kwargs):
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    cfg = ContainerConfig(backend="docker", wave_enabled=True)
    with pytest.raises(containers.ToolExecutionError) as ei:
        resolve_image(_wave_caps(), cfg)
    assert "no image" in ei.value.details()


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


def test_run_tool_with_retries_stops_on_a_permanent_failure(monkeypatch) -> None:
    calls = {"n": 0}

    def no_match(caps, cmd, **kwargs):
        calls["n"] += 1
        raise containers.ToolExecutionError(list(cmd), 1, output="Error: no match")

    monkeypatch.setattr(containers, "run_tool", no_match)
    monkeypatch.setattr(containers.time, "sleep", lambda s: None)
    with pytest.raises(containers.ToolExecutionError):
        containers.run_tool_with_retries(
            _wave_caps(),
            ["datasets", "download"],
            logger=_LOG,
            attempts=3,
            permanent=lambda exc: "no match" in (exc.output or ""),
        )
    assert calls["n"] == 1


# --- backlog fixes from the 2026-09-01 audit self-review ----------------------


def test_unknown_backend_is_a_user_input_error() -> None:
    from repgenr.core.errors import UserInputError

    try:
        with pytest.raises(UserInputError, match="podman"):
            containers.configure_container("podman")
    finally:
        containers.configure_container("none")


def test_wave_cache_is_keyed_by_platform(monkeypatch, _wave_env) -> None:
    import subprocess
    from types import SimpleNamespace

    minted: list[list[str]] = []

    def fake_run(cmd, **kwargs):
        minted.append(list(cmd))
        platform = cmd[cmd.index("--platform") + 1] if "--platform" in cmd else "host"
        return SimpleNamespace(returncode=0, stdout=f"wave.io/img:{platform}\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    amd = resolve_image(
        _wave_caps(), ContainerConfig(backend="docker", wave_enabled=True, platform="linux/amd64")
    )
    arm = resolve_image(
        _wave_caps(), ContainerConfig(backend="docker", wave_enabled=True, platform="linux/arm64")
    )
    assert amd == "wave.io/img:linux/amd64"
    assert arm == "wave.io/img:linux/arm64"
    assert len(minted) == 2  # one mint per platform, not a shared entry
    # The same platform twice is served from the cache.
    resolve_image(
        _wave_caps(), ContainerConfig(backend="docker", wave_enabled=True, platform="linux/amd64")
    )
    assert len(minted) == 2


def test_symlinked_genome_targets_are_bound(tmp_path, monkeypatch) -> None:
    """A genomes/ directory of symlinks (repgenr ingest) needs the targets'
    directory bound too, or the links dangle inside the container."""
    from repgenr.core import containers as c

    source = tmp_path / "elsewhere" / "set"
    source.mkdir(parents=True)
    (source / "a.fasta").write_text(">a\nACGT\n", encoding="utf-8")
    wd = tmp_path / "wd"
    genomes = wd / "genomes"
    genomes.mkdir(parents=True)
    (genomes / "a.fasta").symlink_to(source / "a.fasta")
    monkeypatch.setattr(c.tempfile, "gettempdir", lambda: str(tmp_path / "tmp"))
    (tmp_path / "tmp").mkdir()

    mounts = c._default_mounts(c.ContainerConfig(backend="docker"), wd, [str(genomes)])
    assert any(m == source or str(source).startswith(str(m)) for m in mounts), mounts


def test_container_workdir_defaults_to_the_host_cwd(tmp_path, monkeypatch) -> None:
    """A stateless step run with `-o .` hands the tool relative paths; the
    container must start in the host cwd, not in the temp mount (live audit:
    progressiveMauve under the docker profile could not open align/xmfa/...)."""
    from repgenr.core import containers as c

    sys_tmp = tmp_path / "systmp"
    sys_tmp.mkdir()
    monkeypatch.setattr(c.tempfile, "gettempdir", lambda: str(sys_tmp))
    task = tmp_path / "task"
    task.mkdir()
    monkeypatch.chdir(task)
    cmd = c.wrap_command(
        "img:1",
        ["tool", "align/out.xmfa"],
        config=c.ContainerConfig(backend="docker"),
        cwd=None,
        logger=_LOG,
    )
    assert cmd[cmd.index("-w") + 1] == str(task.resolve()) or cmd[cmd.index("-w") + 1] == str(task)
    assert any(x.startswith(f"{task}:") or x.startswith(f"{task.resolve()}:") for x in cmd)


def test_run_chain_runs_one_container_for_the_whole_sequence(tmp_path, monkeypatch, caplog):
    """Containerized, a per-genome chain is one engine invocation, not one per tool."""
    import logging

    from repgenr.core import containers as mod

    runs: list[list[str]] = []
    monkeypatch.setattr(mod.process, "run", lambda cmd, **kw: runs.append([str(c) for c in cmd]))
    monkeypatch.setattr(mod, "resolve_image", lambda caps, config=None: "example/image:1")
    mod.configure_container(backend="docker")
    try:
        caps = ToolCapabilities(name="chained", container="example/image:1")
        with caplog.at_level(logging.INFO):
            mod.run_chain(
                caps,
                [
                    ("minimap2", ["minimap2", "-o", tmp_path / "a.sam", tmp_path / "ref.fa"]),
                    ("samtools", ["samtools", "sort", "-o", tmp_path / "a.bam"]),
                ],
                logger=logging.getLogger("t"),
            )
    finally:
        mod.configure_container(backend="none")

    assert len(runs) == 1, "one container for the chain"
    script = runs[0][-1]
    assert runs[0][-2] == "-c" and runs[0][-3] == "sh"
    assert script.startswith("set -e\n"), "the chain stops at the first failure"
    assert "minimap2" in script and "samtools sort" in script
    mounts = [runs[0][i + 1].split(":")[0] for i, tok in enumerate(runs[0]) if tok == "-v"]
    covered = any(tmp_path == Path(m) or Path(m) in tmp_path.parents for m in mounts)
    assert covered, f"the paths in the script are mounted: {mounts}"
    assert "[minimap2] $" in caplog.text and "[samtools] $" in caplog.text


def test_run_chain_on_the_host_runs_each_command_itself(tmp_path, monkeypatch):
    """Without a container the chain is just a loop, keeping per-tool logging."""
    import logging

    from repgenr.core import containers as mod

    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(
        mod,
        "run_tool",
        lambda caps, cmd, **kw: calls.append((kw.get("log_prefix"), [str(c) for c in cmd])),
    )
    mod.run_chain(
        ToolCapabilities(name="plain"),
        [("one", ["echo", "1"]), ("two", ["echo", "2"])],
        logger=logging.getLogger("t"),
    )
    assert [prefix for prefix, _ in calls] == ["one", "two"]


def test_mounts_follow_a_symlinked_directory_on_the_path(tmp_path, monkeypatch) -> None:
    """Nextflow stages an input directory as a symlink inside the task dir
    (chunks/c0 -> /work/ab/c0). The genome path a tool receives goes through
    that link; inside the container the link resolves to the real directory,
    which must therefore be bound at its own path as well."""
    import os

    from repgenr.core import containers as c

    sys_tmp = tmp_path / "systmp"
    sys_tmp.mkdir()
    monkeypatch.setattr(c.tempfile, "gettempdir", lambda: str(sys_tmp))
    real = tmp_path / "work" / "ab" / "c0" / "representatives"
    real.mkdir(parents=True)
    (real / "g.fasta").write_text(">g\nACGT\n", encoding="utf-8")
    task = tmp_path / "work" / "cd"
    (task / "chunks").mkdir(parents=True)
    os.symlink(real.parent, task / "chunks" / "c0")
    via_link = task / "chunks" / "c0" / "representatives"

    mounts = c._default_mounts(c.ContainerConfig(backend="docker"), task, [], [str(via_link)])
    assert Path(os.path.realpath(real)) in mounts


def test_containerized_failure_names_the_adapter_tool_not_the_engine(monkeypatch) -> None:
    caps = ToolCapabilities(name="spades.py", container="quay.io/x/spades:1")
    monkeypatch.setattr(containers, "_CONFIG", ContainerConfig(backend="docker"))
    monkeypatch.setattr(
        containers,
        "wrap_command",
        lambda image, argv, **kw: ["sh", "-c", "echo oops >&2; exit 3"],
    )
    with pytest.raises(containers.ToolExecutionError) as ei:
        containers.run_tool(caps, ["spades.py", "-o", "x"], logger=_LOG)
    assert str(ei.value) == "spades.py failed (exit 3)"


def _engine_env(monkeypatch, returncode: int, stderr: str = "") -> list[list[str]]:
    import subprocess

    from repgenr.core import binaries

    calls: list[list[str]] = []

    def run(argv, **kwargs):
        calls.append(list(argv))
        return subprocess.CompletedProcess(argv, returncode, stdout="", stderr=stderr)

    monkeypatch.setattr(containers.subprocess, "run", run)
    monkeypatch.setattr(binaries.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(binaries, "_query_version", lambda name, args: "29.5.3")
    containers._ENGINE_READY.clear()
    return calls


def test_preflight_reports_an_unreachable_docker_daemon(monkeypatch) -> None:
    # `docker --version` answers without a daemon, so list-tools --check said
    # "ok" and the stage failed later with exit 6 and the cause only in the log.
    from repgenr.core.plugins import ToolCapabilities, preflight

    calls = _engine_env(monkeypatch, 1, "failed to connect to the docker API at unix:///x.sock\n")
    caps = ToolCapabilities(name="tool", container="quay.io/x/tool:1")
    try:
        containers.configure_container("docker")
        with pytest.raises(MissingBinaryError, match="daemon is not reachable") as exc:
            preflight(caps)
        assert "failed to connect to the docker API" in str(exc.value)
    finally:
        containers.configure_container("none")
    assert calls == [["docker", "info"]]


def test_engine_readiness_is_checked_once_per_configuration(monkeypatch) -> None:
    from repgenr.core.plugins import ToolCapabilities, preflight

    calls = _engine_env(monkeypatch, 0)
    caps = ToolCapabilities(name="tool", container="quay.io/x/tool:1")
    try:
        containers.configure_container("docker")
        assert preflight(caps) == {"tool": "quay.io/x/tool:1", "docker": "29.5.3"}
        assert preflight(caps) == {"tool": "quay.io/x/tool:1", "docker": "29.5.3"}
        containers.configure_container("singularity")
        preflight(caps)  # no daemon to ask
    finally:
        containers.configure_container("none")
    assert calls == [["docker", "info"]]


def test_preflight_records_the_engine_and_its_version(monkeypatch) -> None:
    # An image-run tool recorded only its image reference, so the engine that
    # ran it was not in repgenr.yaml or versions.yml.
    from repgenr.core.plugins import ToolCapabilities, preflight

    _engine_env(monkeypatch, 0)
    image_caps = ToolCapabilities(name="tool", container="quay.io/x/tool:1")
    try:
        containers.configure_container("docker", engine="/opt/bin/podman")
        assert preflight(image_caps) == {"tool": "quay.io/x/tool:1", "podman": "29.5.3"}
        containers.configure_container("singularity")
        assert preflight(image_caps) == {"tool": "quay.io/x/tool:1", "singularity": "29.5.3"}
    finally:
        containers.configure_container("none")


def test_versions_lists_the_engine_once_for_two_image_stages(tmp_path) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app
    from repgenr.core.config import Config

    cfg = Config()
    cfg.record_stage(
        "dereplicate",
        tool_versions={"skder": "quay.io/x/skder:1", "docker": "29.5.3"},
        completed="t",
    )
    cfg.record_stage(
        "phylo", tool_versions={"iqtree": "quay.io/x/iqtree:2", "docker": "29.5.3"}, completed="t"
    )
    cfg.save(tmp_path)
    result = CliRunner().invoke(app, ["versions", "-wd", str(tmp_path)])
    assert result.exit_code == 0, result.output
    assert result.stdout.splitlines().count("docker: 29.5.3") == 1


def test_an_image_docker_cannot_start_is_named_as_an_engine_failure(monkeypatch) -> None:
    # A pin with a missing tag: docker exits 125 before the tool starts. The
    # message said "sourmash failed (exit 125)", as if sourmash had run.
    caps = ToolCapabilities(name="sourmash", container="quay.io/x/sourmash:0.0.0")
    monkeypatch.setattr(containers, "_CONFIG", ContainerConfig(backend="docker"))
    monkeypatch.setattr(
        containers,
        "wrap_command",
        lambda image, argv, **kw: [
            "sh",
            "-c",
            "echo \"Unable to find image 'quay.io/x/sourmash:0.0.0' locally\" >&2; "
            "echo 'docker: Error response from daemon: manifest unknown' >&2; exit 125",
        ],
    )
    with pytest.raises(containers.ToolExecutionError) as ei:
        containers.run_tool(caps, ["sourmash", "compare"], logger=_LOG)
    msg = str(ei.value)
    assert "docker could not start image quay.io/x/sourmash:0.0.0 for sourmash" in msg
    assert "manifest unknown" in msg
    assert ei.value.exit_code == 6 and ei.value.returncode == 125


@pytest.mark.parametrize("backend", ["docker", "singularity"])
def test_checkm_data_path_reaches_the_container(tmp_path, monkeypatch, backend) -> None:
    # install.md tells dRep users to point CHECKM_DATA_PATH at the CheckM data;
    # docker passed only HOME into the container and mounted nothing there.
    sys_tmp = tmp_path / "systmp"
    sys_tmp.mkdir()
    monkeypatch.setattr(containers.tempfile, "gettempdir", lambda: str(sys_tmp))
    data = tmp_path / "checkm_data"
    data.mkdir()
    monkeypatch.setenv("CHECKM_DATA_PATH", str(data))
    cmd = wrap_command(
        "img:1",
        ["dRep", "dereplicate"],
        config=ContainerConfig(backend=backend),
        cwd="/wd",
        logger=_LOG,
    )
    if backend == "docker":
        assert f"CHECKM_DATA_PATH={data}" in cmd
        assert f"{data}:{data}" in cmd
    else:
        # Singularity passes the host environment; the directory needs a bind.
        assert str(data) in cmd and cmd[cmd.index(str(data)) - 1] == "--bind"


def test_an_unset_checkm_data_path_adds_nothing(monkeypatch) -> None:
    monkeypatch.delenv("CHECKM_DATA_PATH", raising=False)
    cmd = wrap_command(
        "img:1", ["dRep"], config=ContainerConfig(backend="docker"), cwd="/wd", logger=_LOG
    )
    assert not any(c.startswith("CHECKM_DATA_PATH") for c in cmd)


def test_docker_containers_are_named_for_cleanup() -> None:
    cmd = wrap_command(
        "img:1", ["tool"], config=ContainerConfig(backend="docker"), cwd="/wd", logger=_LOG
    )
    name = cmd[cmd.index("--name") + 1]
    assert name.startswith("repgenr-")
    assert cmd.index("--name") < cmd.index("img:1")


@pytest.mark.parametrize("interrupt", [SystemExit(143), KeyboardInterrupt()])
def test_a_stopped_repgenr_stops_the_container_too(monkeypatch, interrupt) -> None:
    # A tool that ignores the SIGTERM forwarded through --init kept its
    # container running after repgenr exited; the engine is asked to stop it.
    caps = ToolCapabilities(name="tool", container="quay.io/x/tool:1")
    monkeypatch.setattr(containers, "_CONFIG", ContainerConfig(backend="docker"))
    calls: list[list[str]] = []

    def run(cmd, **kw):
        raise interrupt

    def engine(argv, **kw):
        import subprocess

        calls.append(list(argv))
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(containers.process, "run", run)
    monkeypatch.setattr(containers.subprocess, "run", engine)
    with pytest.raises(type(interrupt)):
        containers.run_tool(caps, ["tool"], logger=_LOG)
    assert len(calls) == 1
    assert calls[0][:4] == ["docker", "stop", "--time", "0"] and calls[0][-1].startswith("repgenr-")


def test_a_tool_failure_does_not_stop_a_container(monkeypatch) -> None:
    caps = ToolCapabilities(name="tool", container="quay.io/x/tool:1")
    monkeypatch.setattr(containers, "_CONFIG", ContainerConfig(backend="docker"))
    calls: list[list[str]] = []

    def run(cmd, **kw):
        raise containers.ToolExecutionError(cmd, 2, output="bad")

    monkeypatch.setattr(containers.process, "run", run)
    monkeypatch.setattr(containers.subprocess, "run", lambda argv, **kw: calls.append(argv))
    with pytest.raises(containers.ToolExecutionError):
        containers.run_tool(caps, ["tool"], logger=_LOG)
    assert calls == []


def test_a_tool_that_itself_exits_125_is_reported_as_the_tool(monkeypatch) -> None:
    # Only docker's own error marks an engine failure; a tool may exit 125.
    caps = ToolCapabilities(name="sourmash", container="quay.io/x/sourmash:1")
    monkeypatch.setattr(containers, "_CONFIG", ContainerConfig(backend="docker"))
    monkeypatch.setattr(
        containers,
        "wrap_command",
        lambda image, argv, **kw: ["sh", "-c", "echo 'bad k-mer size' >&2; exit 125"],
    )
    with pytest.raises(containers.ToolExecutionError) as ei:
        containers.run_tool(caps, ["sourmash", "compare"], logger=_LOG)
    assert str(ei.value) == "sourmash failed (exit 125)"


def test_a_relative_checkm_data_path_is_made_absolute(tmp_path, monkeypatch) -> None:
    sys_tmp = tmp_path / "systmp"
    sys_tmp.mkdir()
    monkeypatch.setattr(containers.tempfile, "gettempdir", lambda: str(sys_tmp))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "checkm").mkdir()
    monkeypatch.setenv("CHECKM_DATA_PATH", "checkm")
    cmd = wrap_command(
        "img:1", ["dRep"], config=ContainerConfig(backend="docker"), cwd="/wd", logger=_LOG
    )
    data = os.path.abspath(tmp_path / "checkm")
    assert f"CHECKM_DATA_PATH={data}" in cmd
    assert f"{data}:{data}" in cmd


def test_available_cpus_on_the_host_follows_the_affinity_mask(monkeypatch) -> None:
    # A scheduler or taskset may confine repgenr to fewer CPUs than the host has.
    caps = ToolCapabilities(name="gubbins")
    monkeypatch.setattr(containers.os, "sched_getaffinity", lambda pid: {0, 1, 2}, raising=False)
    monkeypatch.setattr(containers.os, "cpu_count", lambda: 64)
    assert containers.available_cpus(caps) == 3
    monkeypatch.delattr(containers.os, "sched_getaffinity", raising=False)
    assert containers.available_cpus(caps) == 64
