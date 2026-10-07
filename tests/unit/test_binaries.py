"""Unit tests for binary preflight + min_version enforcement."""

from __future__ import annotations

import pytest

from repgenr.core import binaries
from repgenr.core.binaries import BinarySpec, check_binaries
from repgenr.core.errors import MissingBinaryError


def _fake_env(monkeypatch, present: set[str], versions: dict[str, str]) -> None:
    monkeypatch.setattr(binaries.shutil, "which", lambda n: n if n in present else None)
    monkeypatch.setattr(binaries, "_query_version", lambda name, args: versions.get(name))


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
    import subprocess

    def run(argv, **kwargs):
        return subprocess.CompletedProcess(argv, returncode, stdout=stdout, stderr="")

    monkeypatch.setattr(binaries.subprocess, "run", run)


def test_rejected_version_flag_is_recorded_as_unknown(monkeypatch) -> None:
    # sibeliaz 1.2.7 has no version flag: `sibeliaz -v` prints
    # "illegal option -- v" and exits 1. That error line is not a version and
    # must not reach repgenr.yaml or `list-tools --check`.
    _fake_run(monkeypatch, 1, "/env/bin/sibeliaz: illegal option -- v\n")
    monkeypatch.setattr(binaries.shutil, "which", lambda n: n)
    assert check_binaries((BinarySpec("sibeliaz", version_args=("-v",)),)) == {
        "sibeliaz": "unknown"
    }


def test_unnumbered_version_line_is_kept_on_success(monkeypatch) -> None:
    _fake_run(monkeypatch, 0, "build abc123\n")
    assert binaries._query_version("tool", ("--version",)) == "build abc123"


def test_strict_version_reports_a_missing_version_without_text(monkeypatch) -> None:
    monkeypatch.setattr(binaries.shutil, "which", lambda n: n)
    monkeypatch.setattr(binaries, "_query_version", lambda name, args: None)
    with pytest.raises(MissingBinaryError, match="could not read a version") as exc:
        check_binaries((BinarySpec("samtools", min_version="1.10", strict_version=True),))
    assert "None" not in str(exc.value)
