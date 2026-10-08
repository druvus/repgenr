"""Auto default for --num-processes (parallel chunk workers).

Benchmark-derived rule: ~4 threads per worker, capped by core count, never
exceeding the thread budget so threads/worker stays >= 1 (no oversubscription).
"""

from __future__ import annotations

from repgenr.stages import dereplicate
from repgenr.stages.dereplicate import _auto_num_processes


def test_auto_targets_four_threads_per_worker(monkeypatch) -> None:
    monkeypatch.setattr(dereplicate, "usable_cpus", lambda: 11)
    assert _auto_num_processes(16) == 4
    assert _auto_num_processes(8) == 2
    assert _auto_num_processes(4) == 1


def test_auto_floors_at_one(monkeypatch) -> None:
    monkeypatch.setattr(dereplicate, "usable_cpus", lambda: 11)
    assert _auto_num_processes(2) == 1  # threads // 4 == 0 -> 1
    assert _auto_num_processes(1) == 1


def test_auto_capped_by_cores(monkeypatch) -> None:
    monkeypatch.setattr(dereplicate, "usable_cpus", lambda: 11)
    assert _auto_num_processes(64) == 11  # 64 // 4 = 16, capped to 11 cores


def test_auto_follows_the_cpu_limit_of_the_process(monkeypatch) -> None:
    # A cgroup quota of 2 CPUs (docker --cpus 2) on an 11-CPU host.
    from repgenr.core import resources

    monkeypatch.setattr(
        resources.os, "sched_getaffinity", lambda pid: set(range(11)), raising=False
    )
    monkeypatch.setattr(resources, "cgroup_cpu_quota", lambda: 2.0)
    assert _auto_num_processes(64) == 2
