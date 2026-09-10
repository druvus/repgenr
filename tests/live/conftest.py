"""Fixtures for the live suite: real `repgenr` processes on real tools.

The tests here are marked ``live`` and are deselected unless ``-m live`` is
given (``addopts`` in ``pyproject.toml``). ``--live-config`` (registered in
``tests/conftest.py``) points at a TOML file that maps tools to bin
directories, prepended to ``PATH`` before test setup so ``requires_binary``
markers see them, and names a cache directory for network fixtures.

Every ``repgenr`` invocation goes through ``subprocess.run`` on the console
script, not ``CliRunner``: global container state, logging flags and
``REPGENR_*`` variables are exercised as a fresh process would see them.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tomllib
from collections.abc import Callable
from pathlib import Path

import pytest

from benchmarks.genomegen import generate_set
from repgenr.core.contracts import SELECTION_TSV, SelectionRow, list_fasta, write_selection

DEFAULT_CACHE_DIR = Path("~/.cache/repgenr-live").expanduser()


def _load_live_config(config: pytest.Config) -> dict:
    path = config.getoption("--live-config")
    if not path:
        return {}
    with open(Path(path).expanduser(), "rb") as fh:
        return tomllib.load(fh)


def pytest_configure(config: pytest.Config) -> None:
    """Prepend configured bin dirs to PATH before any ``requires_binary`` check."""
    live = _load_live_config(config)
    config.stash[_LIVE_KEY] = live
    bin_dirs = [str(Path(p).expanduser()) for p in live.get("bin_dirs", {}).values()]
    if bin_dirs:
        os.environ["PATH"] = os.pathsep.join([*bin_dirs, os.environ.get("PATH", "")])


_LIVE_KEY = pytest.StashKey[dict]()


@pytest.fixture(scope="session")
def live_config(pytestconfig: pytest.Config) -> dict:
    return pytestconfig.stash.get(_LIVE_KEY, {})


@pytest.fixture(scope="session")
def live_cache_dir(live_config: dict) -> Path:
    """Where network fixtures keep their one-time downloads."""
    cache = Path(live_config.get("cache_dir", DEFAULT_CACHE_DIR)).expanduser()
    cache.mkdir(parents=True, exist_ok=True)
    return cache


def _stage_done(workdir: Path, stage: str) -> bool:
    from repgenr.core.config import CONFIG_FILENAME, Config

    if not (workdir / CONFIG_FILENAME).is_file():
        return False
    record = Config.load(workdir).stages.get(stage)
    return bool(record is not None and record.completed)


@pytest.fixture(scope="session")
def cached_workdir(live_cache_dir: Path, repgenr_cmd: list[str]):
    """Build a workdir once under the cache directory and reuse it later.

    ``steps`` are argument lists (without ``-wd``); the workdir is rebuilt
    from scratch unless ``done_stage`` already shows as completed. Network
    fixtures use this so GTDB and NCBI are hit once per machine; tests copy
    the result with :func:`copy_workdir` before touching it.
    """

    def _get(name: str, steps: list[list[str]], done_stage: str) -> Path:
        wd = live_cache_dir / name
        if _stage_done(wd, done_stage):
            return wd
        shutil.rmtree(wd, ignore_errors=True)
        for args in steps:
            argv = [*repgenr_cmd, args[0], "-wd", str(wd), *args[1:]]
            proc = subprocess.run(argv, capture_output=True, text=True, check=False)
            if proc.returncode != 0:
                shutil.rmtree(wd, ignore_errors=True)
                pytest.fail(
                    f"cache build failed ({proc.returncode}): {' '.join(argv)}\n{proc.stderr}"
                )
        return wd

    return _get


def copy_workdir(src: Path, dst: Path) -> Path:
    """Copy a cached workdir (symlinks kept, mtimes preserved so resume holds)."""
    shutil.copytree(src, dst, symlinks=True, ignore=shutil.ignore_patterns("._*", ".DS_Store"))
    return dst


@pytest.fixture(scope="session")
def repgenr_cmd() -> list[str]:
    """Command prefix for the console script (falls back to ``python -m``)."""
    exe = shutil.which("repgenr")
    if exe:
        return [exe]
    return [sys.executable, "-m", "repgenr.cli.main"]


RunRepgenr = Callable[..., subprocess.CompletedProcess[str]]


@pytest.fixture
def run_repgenr(repgenr_cmd: list[str]) -> RunRepgenr:
    """Run ``repgenr <args>`` as a subprocess and return the completed process.

    ``check=True`` (default) fails the test with the captured output when the
    exit status is non-zero. ``env`` entries are added on top of the current
    environment; ``timeout`` is in seconds.
    """

    def _run(
        *args: str | Path,
        check: bool = True,
        env: dict[str, str] | None = None,
        timeout: float = 1800,
    ) -> subprocess.CompletedProcess[str]:
        argv = [*repgenr_cmd, *(str(a) for a in args)]
        merged = {**os.environ, **(env or {})}
        proc = subprocess.run(
            argv, capture_output=True, text=True, env=merged, timeout=timeout, check=False
        )
        if check and proc.returncode != 0:
            pytest.fail(
                f"repgenr exited {proc.returncode}: {' '.join(argv)}\n"
                f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
            )
        return proc

    return _run


@pytest.fixture
def synthetic_set(tmp_path: Path) -> Callable[..., Path]:
    """Generate a seeded synthetic genome set; returns its directory.

    Wraps ``benchmarks.genomegen.generate_set``: ``scenario`` is ``balanced``,
    ``clonal`` or ``mixed``; the truth partition is written to ``truth.json``
    in the same directory.
    """

    def _make(
        scenario: str,
        n: int,
        length: int = 50_000,
        *,
        seed: int = 1,
        order: str = "clustered",
        clone_fraction: float = 0.5,
    ) -> Path:
        out = tmp_path / f"set_{scenario}_{n}"
        generate_set(
            out,
            scenario=scenario,
            n=n,
            seed=seed,
            genome_length=length,
            clone_fraction=clone_fraction,
            order=order,
        )
        return out

    return _make


@pytest.fixture
def selection_for() -> Callable[..., Path]:
    """Write a selection.tsv for a genome directory, optionally with taxonomy and quality.

    ``species`` maps filename to a species token; ``quality`` maps filename to
    ``(completeness, contamination)``. Unlisted genomes keep the values parsed
    from their canonical filename and no quality.
    """

    def _write(
        genomes_dir: Path,
        *,
        species: dict[str, str] | None = None,
        quality: dict[str, tuple[float, float]] | None = None,
        outgroup: str | None = None,
        path: Path | None = None,
    ) -> Path:
        from repgenr.core.contracts import parse_genome_filename

        rows = []
        for f in list_fasta(genomes_dir):
            family, genus, sp, accession = parse_genome_filename(f.name)
            comp, cont = (quality or {}).get(f.name, (None, None))
            rows.append(
                SelectionRow(
                    accession,
                    family,
                    genus,
                    (species or {}).get(f.name, sp),
                    f.name == outgroup,
                    f.name,
                    comp,
                    cont,
                )
            )
        target = path or genomes_dir / SELECTION_TSV
        write_selection(target, rows)
        return target

    return _write


@pytest.fixture
def ingested_workdir(tmp_path: Path, run_repgenr: RunRepgenr) -> Callable[..., Path]:
    """A working directory populated by ``repgenr ingest`` from a genome directory."""

    def _ingest(
        genomes_dir: Path,
        *,
        selection: Path | None = None,
        outgroup: str | None = None,
        copy: bool = False,
        name: str = "wd",
    ) -> Path:
        wd = tmp_path / name
        args: list[str | Path] = ["ingest", "-wd", wd, "--genomes-dir", genomes_dir]
        if selection is not None:
            args += ["--selection", selection]
        if outgroup is not None:
            args += ["--outgroup", outgroup]
        if copy:
            args.append("--copy")
        run_repgenr(*args)
        return wd

    return _ingest
