"""Shared plugin infrastructure: capability metadata + entry-point registry.

Each tool family (dereplicators, aligners, snptypers, treebuilders) defines an
ABC and instantiates a :class:`Registry` bound to its entry-point group. In-tree
adapters and third-party packages are discovered identically through
``importlib.metadata`` entry points, so core never imports a concrete adapter
and external tools need no core edits.
"""

from __future__ import annotations

import shutil
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
    # Container execution: a pinned image URI (BioContainer or Wave-minted) for
    # single-tool adapters, or a conda spec resolved to an image via Wave when no
    # explicit image is set and Wave is enabled.
    container: str | None = None
    conda: tuple[str, ...] = ()


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
    """Check the adapter's required binaries; return resolved versions.

    When a container backend is active and an image resolves for this tool, the
    tool lives in the image (not on the host): check the engine binary instead
    and record the image reference in place of host tool versions.
    """
    from .containers import get_config, resolve_image  # deferred: avoids import cycle

    config = get_config()
    if config.active:
        image = resolve_image(capabilities, config)
        if image:
            check_binaries((BinarySpec(config.engine_binary(), version_args=("--version",)),))
            return {capabilities.name: image}
    return check_binaries(capabilities.required_binaries)


AUTO = "auto"


def _capabilities_of(registry: Registry, name: str) -> ToolCapabilities | None:
    try:
        return registry.get(name).capabilities  # type: ignore[attr-defined]
    except PluginError:
        return None


def _binaries_available(cap: ToolCapabilities) -> bool:
    return all(shutil.which(spec.name) is not None for spec in cap.required_binaries)


def auto_select(registry: Registry, n_items: int) -> str | None:
    """Pick the best-scaling *available* registered tool for ``n_items`` inputs.

    Preference order: required binaries present, then fits the recommended scale,
    then natively-scaling, then the largest (or unbounded) recommended limit,
    then alphabetical for determinism. Preferring installed tools avoids
    auto-selecting an adapter whose binary is missing.
    """
    best: tuple[tuple[int, int, int, float], str] | None = None
    for name in registry.names():
        cap = _capabilities_of(registry, name)
        if cap is None:
            continue
        limit = cap.recommended_max_genomes
        available = 1 if _binaries_available(cap) else 0
        fits = 1 if (limit is None or limit >= n_items) else 0
        native = 1 if cap.supports_native_scaling else 0
        headroom = float("inf") if limit is None else float(limit)
        score = (available, fits, native, headroom)
        if best is None or score > best[0] or (score == best[0] and name < best[1]):
            best = (score, name)
    return best[1] if best else None


def scale_warning(
    registry: Registry, tool: str, n_items: int
) -> tuple[int, list[str]] | None:
    """If ``tool`` is over its recommended scale, return (limit, alternatives).

    Alternatives are registered tools whose recommended scale accommodates
    ``n_items``. Returns None when the tool is within its recommended scale.
    """
    cap = _capabilities_of(registry, tool)
    if cap is None or cap.recommended_max_genomes is None:
        return None
    if n_items <= cap.recommended_max_genomes:
        return None
    alternatives = []
    for name in registry.names():
        if name == tool:
            continue
        other = _capabilities_of(registry, name)
        if other is None:
            continue
        if other.recommended_max_genomes is None or other.recommended_max_genomes >= n_items:
            alternatives.append(name)
    return cap.recommended_max_genomes, sorted(alternatives)
