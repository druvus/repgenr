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


def test_unparseable_version_is_lenient(monkeypatch) -> None:
    # if the version string can't be parsed, min_version is not enforced (no false reject)
    _fake_env(monkeypatch, present={"tool"}, versions={"tool": "weird-build-xyz"})
    out = check_binaries((BinarySpec("tool", min_version="2.0"),))
    assert out["tool"] == "weird-build-xyz"
