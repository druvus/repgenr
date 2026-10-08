"""CPUs this process may use: the affinity mask and a cgroup CPU quota.

A scheduler can limit a process in two ways. A cpuset (Slurm usually, taskset)
shows in the affinity mask. A CPU quota (``docker run --cpus``, a Kubernetes
CPU limit) does not: the process sees every CPU but receives only a fraction of
their time, so a thread count derived from the CPU count oversubscribes the
quota. :func:`usable_cpus` takes the lower of the two.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

PROC_CGROUP = Path("/proc/self/cgroup")
CGROUP_ROOT = Path("/sys/fs/cgroup")


def _quota_v2(root: Path, relative: str) -> float | None:
    """The smallest ``cpu.max`` quota from the process's cgroup up to the root."""
    best: float | None = None
    path = root / relative.lstrip("/")
    while True:
        try:
            fields = (path / "cpu.max").read_text(encoding="utf-8").split()
            if fields and fields[0] != "max":
                quota, period = int(fields[0]), int(fields[1])
                if quota > 0 and period > 0:
                    value = quota / period
                    best = value if best is None else min(best, value)
        except (OSError, ValueError, IndexError):
            pass
        if path == root or root not in path.parents:
            return best
        path = path.parent


def _quota_v1(root: Path, relative: str) -> float | None:
    """``cpu.cfs_quota_us / cpu.cfs_period_us`` of the process's cgroup.

    Only the process's own cgroup is read (or the mount root inside a
    container), so a quota set only on a parent cgroup is not seen. This is
    rare: Docker and Kubernetes set the quota on the container's own cgroup.
    """
    for controller in ("cpu,cpuacct", "cpu"):
        base = root / controller
        # Inside a container's namespace the recorded path may not exist under
        # the mount; the mount root is then the container's own cgroup.
        for directory in (base / relative.lstrip("/"), base):
            try:
                quota = int((directory / "cpu.cfs_quota_us").read_text(encoding="utf-8"))
                period = int((directory / "cpu.cfs_period_us").read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if quota > 0 and period > 0:
                return quota / period
    return None


def cgroup_cpu_quota(proc_cgroup: Path = PROC_CGROUP, root: Path = CGROUP_ROOT) -> float | None:
    """The CPU quota of this process's cgroup in CPUs (1.5 for 150 %), else None.

    cgroup v2 takes the smallest ``cpu.max`` limit on the path from the
    process's cgroup to the root, since a parent's limit applies to its
    children. cgroup v1 reads the cpu controller. None when there is no quota
    or the files cannot be read (macOS, an unusual mount layout).
    """
    try:
        lines = proc_cgroup.read_text(encoding="utf-8").splitlines()
    except OSError:
        return None
    for line in lines:
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        hierarchy, controllers, relative = parts
        if hierarchy == "0" and controllers == "":
            return _quota_v2(root, relative)
        if "cpu" in controllers.split(","):
            return _quota_v1(root, relative)
    return None


def usable_cpus() -> int:
    """CPUs this process may use: the affinity mask (or the CPU count), capped by
    the cgroup CPU quota rounded up; at least 1."""
    affinity = getattr(os, "sched_getaffinity", None)  # absent on macOS
    try:
        count = (len(affinity(0)) if affinity else 0) or os.cpu_count() or 1
    except OSError:
        count = os.cpu_count() or 1
    quota = cgroup_cpu_quota()
    if quota is not None:
        count = min(count, max(1, math.ceil(quota)))
    return max(1, count)
