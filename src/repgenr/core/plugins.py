"""Shared plugin infrastructure: capability metadata + entry-point registry.

Each tool family (dereplicators, aligners, snptypers, treebuilders) defines an
ABC and instantiates a :class:`Registry` bound to its entry-point group. In-tree
adapters and third-party packages are discovered identically through
``importlib.metadata`` entry points, so core never imports a concrete adapter
and external tools need no core edits.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib.metadata import entry_points

from .binaries import BinarySpec, check_binaries
from .errors import PluginError


@dataclass(frozen=True)
class ToolCapabilities:
    """Declarative metadata for one tool adapter.

    ``recommended_max_genomes`` drives auto-selection and scale warnings.
    ``supports_native_scaling`` lets a dereplicator opt out of the chunking
    wrapper (e.g. skDER scales natively).
    """

    name: str
    required_binaries: tuple[BinarySpec, ...] = ()
    default_params: dict = field(default_factory=dict)
    recommended_max_genomes: int | None = None
    supports_native_scaling: bool = False
    threads_param: str | None = None


class Registry[T]:
    """Lazily-loaded registry of adapter classes for one entry-point group."""

    def __init__(self, group: str):
        self.group = group
        self._classes: dict[str, type[T]] = {}
        self._loaded = False

    def _load(self) -> None:
        if self._loaded:
            return
        for ep in entry_points(group=self.group):
            try:
                self._classes[ep.name] = ep.load()
            except Exception as exc:  # a broken third-party plugin must not kill the run
                # Deferred: surfaced only if the broken name is actually requested.
                self._classes.setdefault(ep.name, _BrokenPlugin(ep.name, exc))  # type: ignore[arg-type]
        self._loaded = True

    def names(self) -> list[str]:
        self._load()
        return sorted(self._classes)

    def get(self, name: str) -> type[T]:
        self._load()
        if name not in self._classes:
            available = ", ".join(self.names()) or "none"
            raise PluginError(
                f"Unknown tool '{name}' for {self.group}. Available: {available}"
            )
        cls = self._classes[name]
        if isinstance(cls, _BrokenPlugin):
            raise PluginError(f"Plugin '{name}' failed to load: {cls.error}") from cls.error
        return cls

    def create(self, name: str) -> T:
        return self.get(name)()  # type: ignore[call-arg]


class _BrokenPlugin:
    """Placeholder for an entry point that failed to import."""

    def __init__(self, name: str, error: Exception):
        self.name = name
        self.error = error


def preflight(capabilities: ToolCapabilities) -> dict[str, str]:
    """Check the adapter's required binaries; return resolved versions."""
    return check_binaries(capabilities.required_binaries)
