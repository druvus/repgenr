"""Container execution backend for external tools.

RepGenR drives bioinformatics tools via :func:`repgenr.core.process.run`. When a
container backend is active, each tool's argument vector is wrapped so it runs
inside a pinned image (a BioContainer, the Cactus image, or a Wave-built image),
which unblocks Linux-only tools on any host and pins tool versions.

Only the external-tool subprocess is containerized; RepGenR's own Python (e.g.
the format converters) keeps running on the host and reads tool outputs from the
bind-mounted working directory.

The backend is process-global: configure it once (CLI/env), then adapters call
:func:`run_tool` with their capabilities.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from . import process
from .errors import MissingBinaryError, ToolExecutionError
from .plugins import ToolCapabilities

NATIVE = "none"
DOCKER = "docker"
SINGULARITY = "singularity"


@dataclass
class ContainerConfig:
    backend: str = NATIVE
    engine: str | None = None  # explicit engine binary; defaults per backend
    platform: str | None = None  # e.g. "linux/amd64" for emulated BioContainers
    cache_dir: Path | None = None  # where Singularity .sif / Wave cache live
    wave_enabled: bool = False
    extra_mounts: tuple[Path, ...] = field(default_factory=tuple)

    @property
    def active(self) -> bool:
        return self.backend != NATIVE

    def engine_binary(self) -> str:
        if self.engine:
            return self.engine
        return "docker" if self.backend == DOCKER else "singularity"


_CONFIG = ContainerConfig()
_WAVE_CACHE: dict[tuple[str, ...], str] = {}


def configure_container(
    backend: str = NATIVE,
    *,
    engine: str | None = None,
    platform: str | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
    wave_enabled: bool = False,
) -> ContainerConfig:
    """Set the process-global container backend. ``backend='none'`` = native."""
    global _CONFIG
    if backend not in (NATIVE, DOCKER, SINGULARITY):
        raise ToolExecutionError([backend], 2, f"unknown container backend '{backend}'")
    _CONFIG = ContainerConfig(
        backend=backend,
        engine=engine,
        platform=platform,
        cache_dir=Path(cache_dir).resolve() if cache_dir else None,
        wave_enabled=wave_enabled,
    )
    return _CONFIG


def get_config() -> ContainerConfig:
    return _CONFIG


def resolve_image(caps: ToolCapabilities, config: ContainerConfig | None = None) -> str | None:
    """Return the image URI for an adapter, or None to run natively.

    Prefers an explicit ``caps.container``; otherwise, when Wave is enabled and a
    ``caps.conda`` spec exists, mints (and caches) an image via the Wave CLI.
    """
    config = config or _CONFIG
    if caps.container:
        return caps.container
    if config.wave_enabled and caps.conda:
        return _wave_image(caps.conda, config)
    return None


def _wave_image(conda_spec: tuple[str, ...], config: ContainerConfig) -> str:
    if conda_spec in _WAVE_CACHE:
        return _WAVE_CACHE[conda_spec]
    if shutil.which("wave") is None:
        raise MissingBinaryError(
            "Wave is enabled but the 'wave' CLI is not on PATH. Install it or pin "
            f"an explicit image for conda spec {' '.join(conda_spec)}."
        )
    cmd = ["wave"]
    for pkg in conda_spec:
        cmd += ["--conda-package", pkg]
    if config.platform:
        cmd += ["--platform", config.platform]
    import subprocess

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        raise ToolExecutionError(cmd, proc.returncode, output=proc.stderr.strip())
    image = proc.stdout.strip().splitlines()[-1].strip()
    _WAVE_CACHE[conda_spec] = image
    return image


def _engine_env(config: ContainerConfig) -> dict[str, str]:
    """Env additions so Singularity/Apptainer caches live in cache_dir."""
    if config.backend != SINGULARITY or config.cache_dir is None:
        return {}
    cache = str(config.cache_dir)
    tmp = str(config.cache_dir / "tmp")
    Path(tmp).mkdir(parents=True, exist_ok=True)
    return {
        "SINGULARITY_CACHEDIR": cache,
        "APPTAINER_CACHEDIR": cache,
        "SINGULARITY_TMPDIR": tmp,
        "APPTAINER_TMPDIR": tmp,
    }


def _default_mounts(
    config: ContainerConfig,
    cwd: str | os.PathLike[str] | None,
    argv: Sequence[str] = (),
) -> list[Path]:
    # Use absolute (NOT symlink-resolved) paths: tools are given paths like
    # macOS /var/folders/... and /Users/...; resolving to /private/var or
    # /System/Volumes/Data would mount at paths that do not exist inside the
    # Linux container. Docker resolves host-side symlinks itself.
    def absp(p: str | os.PathLike[str]) -> Path:
        return Path(os.path.abspath(p))

    mounts: list[Path] = []
    if cwd is not None:
        mounts.append(absp(cwd))
    mounts.append(absp(tempfile.gettempdir()))
    mounts.extend(absp(m) for m in config.extra_mounts)
    # Bind the directories referenced by absolute-path arguments (genome inputs,
    # output dirs, references), so the tool sees its files at identical paths.
    for token in argv:
        if token.startswith("/"):
            p = Path(token)
            d = p if p.is_dir() else p.parent
            if d.exists():
                mounts.append(absp(d))
    # de-duplicate, dropping any mount nested under another
    unique: list[Path] = []
    for m in sorted(set(mounts), key=lambda p: len(str(p))):
        if not any(_is_relative_to(m, u) for u in unique):
            unique.append(m)
    return unique


def _is_relative_to(child: Path, parent: Path) -> bool:
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False


def _sanitize(image: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", image)


def _singularity_source(image: str, config: ContainerConfig, logger: logging.Logger) -> str:
    """Return a Singularity image reference, pre-pulling a .sif into cache_dir.

    A local ``.sif`` is used as-is. Otherwise the image is treated as a registry
    reference (``docker://`` is assumed when no scheme is present). With a cache
    dir set, the image is pulled once to ``<cache>/<name>.sif`` and reused.
    """
    if image.endswith(".sif"):
        return image
    source = image if "://" in image else f"docker://{image}"
    if config.cache_dir is None:
        return source
    config.cache_dir.mkdir(parents=True, exist_ok=True)
    sif = config.cache_dir / f"{_sanitize(image)}.sif"
    if not sif.exists():
        logger.info("Pulling %s -> %s", source, sif)
        process.run(
            [config.engine_binary(), "pull", str(sif), source],
            logger=logger,
            env=_engine_env(config),
            log_prefix="singularity",
        )
    return str(sif)


def wrap_command(
    image: str,
    argv: Sequence[str],
    *,
    config: ContainerConfig,
    cwd: str | os.PathLike[str] | None,
    logger: logging.Logger,
) -> list[str]:
    """Build the engine command that runs ``argv`` inside ``image``."""
    mounts = _default_mounts(config, cwd, argv)
    workdir = str(Path(os.path.abspath(cwd))) if cwd is not None else str(mounts[0])

    if config.backend == DOCKER:
        cmd = [config.engine_binary(), "run", "--rm", "--entrypoint", ""]
        cmd += ["-u", f"{os.getuid()}:{os.getgid()}"]
        if config.platform:
            cmd += ["--platform", config.platform]
        for m in mounts:
            cmd += ["-v", f"{m}:{m}"]
        cmd += ["-w", workdir, image, *argv]
        return cmd

    # singularity / apptainer
    source = _singularity_source(image, config, logger)
    cmd = [config.engine_binary(), "exec"]
    for m in mounts:
        cmd += ["--bind", str(m)]
    cmd += ["--pwd", workdir, source, *argv]
    return cmd


def run_tool(
    caps: ToolCapabilities,
    command: Sequence[str | os.PathLike[str]],
    *,
    logger: logging.Logger,
    cwd: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    stdout_path: str | os.PathLike[str] | None = None,
    log_prefix: str | None = None,
) -> int:
    """Run an adapter's tool command, containerized when a backend is active."""
    config = _CONFIG
    image = resolve_image(caps, config) if config.active else None
    if image is None:
        return process.run(
            command, logger=logger, cwd=cwd, env=env, check=check,
            stdout_path=stdout_path, log_prefix=log_prefix,
        )

    argv = [str(part) for part in command]
    wrapped = wrap_command(image, argv, config=config, cwd=cwd, logger=logger)
    merged_env = {**_engine_env(config), **(dict(env) if env else {})} or None
    return process.run(
        wrapped, logger=logger, cwd=cwd, env=merged_env, check=check,
        stdout_path=stdout_path, log_prefix=log_prefix or caps.name,
    )
