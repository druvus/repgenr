"""cgroup CPU quota and usable CPU count, on fake cgroup trees."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core import resources
from repgenr.core.resources import cgroup_cpu_quota


def _v2(tmp_path: Path, own: str, limits: dict[str, str]) -> tuple[Path, Path]:
    proc = tmp_path / "cgroup"
    proc.write_text(f"0::{own}\n")
    root = tmp_path / "fs"
    for rel, value in limits.items():
        (root / rel).mkdir(parents=True, exist_ok=True)
        (root / rel / "cpu.max").write_text(value + "\n")
    (root / own.lstrip("/")).mkdir(parents=True, exist_ok=True)
    return proc, root


def _v1(tmp_path: Path, quota: str, period: str = "100000") -> tuple[Path, Path]:
    proc = tmp_path / "cgroup"
    proc.write_text("5:memory:/slurm/job1\n4:cpu,cpuacct:/slurm/job1\n")
    root = tmp_path / "fs"
    d = root / "cpu,cpuacct" / "slurm" / "job1"
    d.mkdir(parents=True)
    (d / "cpu.cfs_quota_us").write_text(quota + "\n")
    (d / "cpu.cfs_period_us").write_text(period + "\n")
    return proc, root


def test_v2_takes_the_smallest_limit_on_the_path_to_the_root(tmp_path: Path) -> None:
    proc, root = _v2(
        tmp_path,
        "/kubepods/pod1/c1",
        {"": "max 100000", "kubepods/pod1": "400000 100000", "kubepods/pod1/c1": "150000 100000"},
    )
    assert cgroup_cpu_quota(proc, root) == 1.5


def test_v2_a_parent_limit_applies_to_its_child(tmp_path: Path) -> None:
    proc, root = _v2(tmp_path, "/a/b", {"a": "200000 100000", "a/b": "max 100000"})
    assert cgroup_cpu_quota(proc, root) == 2.0


def test_v2_max_everywhere_is_no_quota(tmp_path: Path) -> None:
    proc, root = _v2(tmp_path, "/a", {"": "max 100000", "a": "max 100000"})
    assert cgroup_cpu_quota(proc, root) is None


def test_v2_in_a_container_namespace_reads_the_mount_root(tmp_path: Path) -> None:
    # docker run --cpus 2: /proc/self/cgroup reads "0::/".
    proc, root = _v2(tmp_path, "/", {"": "200000 100000"})
    assert cgroup_cpu_quota(proc, root) == 2.0


def test_v1_reads_the_cpu_controller(tmp_path: Path) -> None:
    proc, root = _v1(tmp_path, "300000")
    assert cgroup_cpu_quota(proc, root) == 3.0


def test_v1_minus_one_is_no_quota(tmp_path: Path) -> None:
    proc, root = _v1(tmp_path, "-1")
    assert cgroup_cpu_quota(proc, root) is None


@pytest.mark.parametrize("value", ["garbage", "", "abc 100000", "100000"])
def test_a_malformed_cpu_max_is_ignored(tmp_path: Path, value: str) -> None:
    proc, root = _v2(tmp_path, "/a", {"a": value})
    assert cgroup_cpu_quota(proc, root) is None


def test_a_malformed_proc_cgroup_line_is_ignored(tmp_path: Path) -> None:
    proc = tmp_path / "cgroup"
    proc.write_text("nonsense\n")
    assert cgroup_cpu_quota(proc, tmp_path) is None


def test_no_proc_cgroup_is_no_quota(tmp_path: Path) -> None:
    # macOS has no /proc.
    assert cgroup_cpu_quota(tmp_path / "missing", tmp_path) is None


def _host(monkeypatch, cpus: int, quota: float | None) -> None:
    monkeypatch.setattr(
        resources.os, "sched_getaffinity", lambda pid: set(range(cpus)), raising=False
    )
    monkeypatch.setattr(resources, "cgroup_cpu_quota", lambda: quota)


@pytest.mark.parametrize(
    ("cpus", "quota", "expected"),
    [(11, 2.0, 2), (11, 1.5, 2), (11, 0.25, 1), (3, 8.0, 3), (12, None, 12)],
)
def test_usable_cpus_is_the_affinity_capped_by_the_quota(
    monkeypatch, cpus: int, quota: float | None, expected: int
) -> None:
    _host(monkeypatch, cpus, quota)
    assert resources.usable_cpus() == expected


def test_usable_cpus_without_affinity_uses_the_cpu_count(monkeypatch) -> None:
    monkeypatch.delattr(resources.os, "sched_getaffinity", raising=False)
    monkeypatch.setattr(resources.os, "cpu_count", lambda: 6)
    monkeypatch.setattr(resources, "cgroup_cpu_quota", lambda: None)
    assert resources.usable_cpus() == 6
