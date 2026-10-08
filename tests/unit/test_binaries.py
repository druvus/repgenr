"""Unit tests for binary preflight + min_version enforcement."""

from __future__ import annotations

import pytest

from repgenr.core import binaries
from repgenr.core.binaries import BinarySpec, check_binaries
from repgenr.core.errors import MissingBinaryError


def _fake_env(monkeypatch, present: set[str], versions: dict[str, str]) -> None:
    monkeypatch.setattr(binaries.shutil, "which", lambda n, path=None: n if n in present else None)
    monkeypatch.setattr(
        binaries, "_query_version", lambda name, args, timeout=None, path=None: versions.get(name)
    )


def test_missing_binary_raises(monkeypatch) -> None:
    _fake_env(monkeypatch, present=set(), versions={})
    with pytest.raises(MissingBinaryError, match="not found on PATH"):
        check_binaries((BinarySpec("skder"),))


def test_too_old_version_rejected(monkeypatch) -> None:
    # the exact failure mode: an ancient samtools 0.1.19 below the 1.10 floor
    _fake_env(monkeypatch, present={"samtools"}, versions={"samtools": "0.1.19"})
    with pytest.raises(MissingBinaryError, match="0.1.19 < required 1.10"):
        check_binaries((BinarySpec("samtools", min_version="1.10"),))


def test_new_enough_version_passes(monkeypatch) -> None:
    # _query_version normalizes to a dotted version; mock that contract.
    _fake_env(monkeypatch, present={"samtools"}, versions={"samtools": "1.23"})
    out = check_binaries((BinarySpec("samtools", min_version="1.10"),))
    assert out["samtools"] == "1.23"


def test_unparseable_version_is_lenient_by_default(monkeypatch) -> None:
    # if the version string can't be parsed, a non-strict min_version is not
    # enforced (no false reject for tools with unusual version output)
    _fake_env(monkeypatch, present={"tool"}, versions={"tool": "weird-build-xyz"})
    out = check_binaries((BinarySpec("tool", min_version="2.0"),))
    assert out["tool"] == "weird-build-xyz"


def test_strict_version_rejects_unparseable(monkeypatch) -> None:
    # the real footgun: ancient samtools answers --version with an error (no
    # digits); a strict_version spec must reject it, not silently pass.
    _fake_env(
        monkeypatch,
        present={"samtools"},
        versions={"samtools": "[main] unrecognized command '--version'"},
    )
    with pytest.raises(MissingBinaryError, match="could not read a version"):
        check_binaries((BinarySpec("samtools", min_version="1.10", strict_version=True),))


def _fake_run(monkeypatch, returncode: int, stdout: str) -> None:
    monkeypatch.setattr(binaries, "_ask", lambda argv, timeout: (returncode, stdout))


def test_rejected_version_flag_is_recorded_as_unknown(monkeypatch) -> None:
    # sibeliaz 1.2.7 has no version flag: `sibeliaz -v` prints
    # "illegal option -- v" and exits 1. That error line is not a version and
    # must not reach repgenr.yaml or `list-tools --check`.
    _fake_run(monkeypatch, 1, "/env/bin/sibeliaz: illegal option -- v\n")
    monkeypatch.setattr(binaries.shutil, "which", lambda n, path=None: n)
    assert check_binaries((BinarySpec("sibeliaz", version_args=("-v",)),)) == {
        "sibeliaz": "unknown"
    }


def test_unnumbered_version_line_is_kept_on_success(monkeypatch) -> None:
    _fake_run(monkeypatch, 0, "build abc123\n")
    assert binaries._query_version("tool", ("--version",)) == "build abc123"


def test_strict_version_reports_a_missing_version_without_text(monkeypatch) -> None:
    monkeypatch.setattr(binaries.shutil, "which", lambda n, path=None: n)
    monkeypatch.setattr(binaries, "_query_version", lambda name, args, timeout=None, path=None: None)
    with pytest.raises(MissingBinaryError, match="could not read a version") as exc:
        check_binaries((BinarySpec("samtools", min_version="1.10", strict_version=True),))
    assert "None" not in str(exc.value)


def test_a_crashed_version_query_is_not_read_from_its_traceback(monkeypatch) -> None:
    # cactus-pangenome in a broken environment answers --version with a Python
    # traceback (exit 1). The interpreter path in it ("python3.12") was read as
    # version 3.12.0, which passed the 2.5 floor and reached repgenr.yaml.
    traceback = (
        "Traceback (most recent call last):\n"
        '  File "/env/lib/python3.12/site-packages/cactus/refmap/cactus_graphmap.py", '
        "line 7, in <module>\n"
        "ImportError: cannot import name 'x'\n"
    )
    _fake_run(monkeypatch, 1, traceback)
    assert binaries._query_version("cactus-pangenome", ("--version",)) is None


def test_versions_inside_paths_and_identifiers_are_not_matched() -> None:
    assert binaries._parse_version("/env/lib/python3.12/site-packages/x.py") is None
    assert binaries._parse_version("GLIBC_2.17 not found") is None
    # The forms real tools print still parse.
    assert binaries._parse_version("hal2maf v2.2: Convert hal database") == (2, 2, 0)
    assert binaries._parse_version("dRep v3.4.5 :::") == (3, 4, 5)
    assert binaries._parse_version("RAxML-NG v. 2.0.2 released") == (2, 0, 2)
    assert binaries._parse_version("2.9.6-b1802") == (2, 9, 6)
    assert binaries._parse_version("Version of skDER being used is: 1.3.6") == (1, 3, 6)


def _conda_prefix(tmp_path, records: dict[str, str]):
    """A fake conda prefix with ``bin/sibeliaz`` and the given conda-meta files."""
    (tmp_path / "bin").mkdir()
    exe = tmp_path / "bin" / "sibeliaz"
    exe.write_text("#!/bin/sh\nexit 1\n", encoding="utf-8")
    exe.chmod(0o755)
    (tmp_path / "conda-meta").mkdir()
    for fname, text in records.items():
        (tmp_path / "conda-meta" / fname).write_text(text, encoding="utf-8")
    return exe


def test_version_is_read_from_conda_meta_when_the_tool_reports_none(monkeypatch, tmp_path) -> None:
    # sibeliaz has no version flag; its conda package record names 1.2.7.
    import json

    exe = _conda_prefix(
        tmp_path,
        {
            "sibeliaz-1.2.7-h28ef24b_0.json": json.dumps(
                {"name": "sibeliaz", "version": "1.2.7", "files": ["bin/sibeliaz"]}
            )
        },
    )
    monkeypatch.setattr(binaries.shutil, "which", lambda n: str(exe))
    monkeypatch.setattr(binaries, "_query_version", lambda name, args, timeout=None: None)
    spec = BinarySpec("sibeliaz", version_args=("-v",), min_version="1.2")
    assert check_binaries((spec,)) == {"sibeliaz": "1.2.7"}


def test_conda_meta_record_is_found_by_its_file_list(monkeypatch, tmp_path) -> None:
    import json

    exe = _conda_prefix(
        tmp_path,
        {
            "other-1.0-0.json": json.dumps({"version": "1.0", "files": ["bin/other"]}),
            "sibeliaz-suite-2.0-0.json": json.dumps(
                {"version": "2.0", "files": ["bin/sibeliaz", "bin/sibeliaz-lcb"]}
            ),
        },
    )
    assert binaries._metadata_version("sibeliaz", str(exe)) == "2.0"


def test_tool_reported_version_is_preferred_over_conda_meta(monkeypatch, tmp_path) -> None:
    import json

    exe = _conda_prefix(
        tmp_path,
        {"sibeliaz-1.2.7-0.json": json.dumps({"version": "1.2.7", "files": ["bin/sibeliaz"]})},
    )
    monkeypatch.setattr(binaries.shutil, "which", lambda n: str(exe))
    monkeypatch.setattr(binaries, "_query_version", lambda name, args, timeout=None: "1.3.0")
    assert check_binaries((BinarySpec("sibeliaz"),)) == {"sibeliaz": "1.3.0"}


@pytest.mark.parametrize(
    "records",
    [
        {},  # no record for the binary
        {"sibeliaz-1.2.7-0.json": '{"version": "1.2.7", "files": ["bin/sibeliaz"'},  # corrupt
        {"sibeliaz-1.2.7-0.json": '{"files": ["bin/sibeliaz"]}'},  # no version field
    ],
)
def test_missing_or_corrupt_conda_meta_gives_unknown(monkeypatch, tmp_path, records) -> None:
    exe = _conda_prefix(tmp_path, records)
    monkeypatch.setattr(binaries.shutil, "which", lambda n: str(exe))
    monkeypatch.setattr(binaries, "_query_version", lambda name, args, timeout=None: None)
    assert check_binaries((BinarySpec("sibeliaz"),)) == {"sibeliaz": "unknown"}


def test_no_conda_meta_directory_gives_no_metadata_version(tmp_path) -> None:
    (tmp_path / "bin").mkdir()
    exe = tmp_path / "bin" / "tool"
    exe.write_text("", encoding="utf-8")
    assert binaries._metadata_version("tool", str(exe)) is None


def _hanging_tool(tmp_path):
    """A version query that never answers and leaves a helper holding its output."""
    script = tmp_path / "hangtool"
    pidfile = tmp_path / "helper.pid"
    script.write_text(f"#!/bin/sh\nsleep 300 &\necho $! > {pidfile}\nsleep 300\n", encoding="utf-8")
    script.chmod(0o755)
    return script, pidfile


def _alive(pid: int) -> bool:
    import os

    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def test_a_hanging_version_query_is_stopped_with_its_helpers(tmp_path, caplog) -> None:
    import logging
    import time

    script, pidfile = _hanging_tool(tmp_path)
    start = time.monotonic()
    with caplog.at_level(logging.WARNING, logger="repgenr"):
        assert binaries._query_version(str(script), ("--version",), timeout=1) is None
    # The background helper keeps the output pipe open; without stopping its
    # process group the query waited for the helper (300 s) after the timeout.
    assert time.monotonic() - start < 15
    helper = int(pidfile.read_text(encoding="utf-8"))
    deadline = time.monotonic() + 5
    while _alive(helper) and time.monotonic() < deadline:
        time.sleep(0.1)
    assert not _alive(helper)
    warning = next(r for r in caplog.records if r.levelno == logging.WARNING)
    assert "hangtool --version did not answer within 1 s" in warning.getMessage()
    assert "unless its conda package record names one" in warning.getMessage()


def test_an_interrupted_version_query_is_stopped_with_its_helpers(tmp_path) -> None:
    # The query runs in its own session, so Ctrl-C at the terminal does not
    # reach it; the KeyboardInterrupt in repgenr must stop its group.
    import signal
    import time

    from repgenr.core import process

    script, pidfile = _hanging_tool(tmp_path)
    started: list = []
    real_register = process.register_live

    def register(proc):
        started.append(proc)
        real_register(proc)

    def interrupt(signum, frame):
        raise KeyboardInterrupt

    previous = signal.signal(signal.SIGALRM, interrupt)
    try:
        process.register_live = register  # type: ignore[assignment]
        signal.setitimer(signal.ITIMER_REAL, 0.5)
        with pytest.raises(KeyboardInterrupt):
            binaries._query_version(str(script), ("--version",), timeout=60)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)
        process.register_live = real_register  # type: ignore[assignment]
    (proc,) = started
    assert proc not in process._live
    helper = int(pidfile.read_text(encoding="utf-8"))
    deadline = time.monotonic() + 5
    while _alive(helper) and time.monotonic() < deadline:
        time.sleep(0.1)
    assert proc.poll() is not None
    assert not _alive(helper)


def test_a_running_version_query_is_reached_by_stop_running_tools(tmp_path) -> None:
    # The termination handler signals every tool in process._live; a version
    # query must be among them while it runs.
    import threading
    import time

    from repgenr.core import process

    script, pidfile = _hanging_tool(tmp_path)
    result: list = []
    worker = threading.Thread(
        target=lambda: result.append(
            binaries._query_version(str(script), ("--version",), timeout=60)
        )
    )
    worker.start()
    deadline = time.monotonic() + 5
    while not pidfile.exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    assert process.stop_running_tools() >= 1
    worker.join(timeout=10)
    assert not worker.is_alive()
    assert result == [None]  # killed by the signal: exit != 0, no output
    assert not _alive(int(pidfile.read_text(encoding="utf-8")))


def test_check_binaries_passes_its_timeout_and_reports_unknown(monkeypatch) -> None:
    seen: list[float | None] = []

    def query(name, args, timeout=None):
        seen.append(timeout)
        return None

    monkeypatch.setattr(binaries.shutil, "which", lambda n: n)
    monkeypatch.setattr(binaries, "_query_version", query)
    spec = BinarySpec("tool")
    assert check_binaries((spec,), timeout=8) == {"tool": "unknown"}
    assert check_binaries((spec,)) == {"tool": "unknown"}
    with binaries.version_timeout(3):
        check_binaries((spec,))
    check_binaries((spec,))
    assert seen == [8, binaries.VERSION_TIMEOUT, 3, binaries.VERSION_TIMEOUT]
    # A strict floor fails with a message rather than hanging.
    strict = BinarySpec("tool", min_version="1.0", strict_version=True)
    with pytest.raises(MissingBinaryError, match="no version reported"):
        check_binaries((strict,), timeout=8)


def test_list_tools_check_queries_versions_with_a_short_timeout(monkeypatch) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    seen: set[float | None] = set()

    def query(name, args, timeout=None):
        seen.add(timeout)
        return "1.0"

    monkeypatch.setattr(binaries.shutil, "which", lambda n: f"/usr/bin/{n}")
    monkeypatch.setattr(binaries, "_query_version", query)
    result = CliRunner().invoke(app, ["list-tools", "--check"])
    assert result.exit_code == 0, result.output
    assert seen == {8}
