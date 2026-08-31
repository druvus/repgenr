"""Shared plugin infrastructure: capability metadata + entry-point registry.

Each tool family (dereplicators, aligners, snptypers, treebuilders) defines an
ABC and instantiates a :class:`Registry` bound to its entry-point group. In-tree
adapters and third-party packages are discovered identically through
``importlib.metadata`` entry points, so core never imports a concrete adapter
and external tools need no core edits.
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Mapping
from dataclasses import dataclass, field
from importlib.metadata import entry_points

from .binaries import BinarySpec, check_binaries
from .errors import PluginError, UserInputError


def parse_extra_int(extra: Mapping[str, object], key: str, default: int) -> int:
    """Read an integer tool override from ``--extra`` with a clear error.

    A non-integer value (e.g. ``--extra ksize=big``) raises
    :class:`UserInputError` naming the key instead of a raw ``ValueError``
    traceback from ``int()``.
    """
    raw = extra.get(key, default)
    try:
        return int(raw)  # type: ignore[call-overload]
    except (TypeError, ValueError) as exc:
        raise UserInputError(
            f"Extra parameter '{key}' must be an integer, got {raw!r}."
        ) from exc


@dataclass(frozen=True)
class ToolCapabilities:
    """Declarative metadata for one tool adapter.

    ``recommended_max_genomes`` drives auto-selection and scale warnings.
    ``supports_native_scaling`` marks a dereplicator that scales to large sets in
    one pass (e.g. skDER), so it is single-pass by default; it no longer gates
    chunking -- ``--process-size`` opts any tool into the two-stage chunked path.
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
    # Extra-dict keys this adapter actually reads; used to warn on (and avoid
    # injecting) tuning that a tool would silently ignore.
    accepted_extras: frozenset[str] = frozenset()


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
                # Deferred: surfaced only if the broken name is actually requested,
                # but log at debug so a broken in-tree adapter is diagnosable.
                logging.getLogger("repgenr").warning(
                    "Plugin %r (group %s) failed to load and is unavailable: %s",
                    ep.name, self.group, exc,
                )
                self._classes.setdefault(ep.name, _BrokenPlugin(ep.name, exc))  # type: ignore[arg-type]
        self._loaded = True

    def names(self) -> list[str]:
        self._load()
        return sorted(self._classes)

    def is_broken(self, name: str) -> bool:
        """True when ``name`` is registered but its adapter failed to import."""
        self._load()
        return isinstance(self._classes.get(name), _BrokenPlugin)

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


def _tool_available(cap: ToolCapabilities) -> bool:
    """Is this adapter runnable in the CURRENT execution environment?

    Under an active container backend the tool lives in an image, not on the
    host, so a declared ``container`` or ``conda`` spec counts as available
    (resolution itself happens at preflight); natively, the required binaries
    must be on PATH.
    """
    from .containers import get_config

    if get_config().active:
        return cap.container is not None or bool(cap.conda)
    return all(shutil.which(spec.name) is not None for spec in cap.required_binaries)


# Tie-break order for auto-selection, matching the documented per-family
# defaults; unlisted tools rank after these, alphabetically.
_PREFERRED_ORDER = ("skder", "iqtree", "progressivemauve", "simple")


def _preference_rank(name: str) -> int:
    try:
        return _PREFERRED_ORDER.index(name)
    except ValueError:
        return len(_PREFERRED_ORDER)


def auto_select(registry: Registry, n_items: int) -> str | None:
    """Pick the best *available* registered tool for ``n_items`` inputs.

    Preference order: runnable in the current environment (container-aware),
    then fits the recommended scale, then the TIGHTEST fitting limit --
    higher-quality tools declare tighter recommended scales, so the loosest
    tool must not win small inputs -- with unbounded tools last among the
    fitting ones. Ties follow the documented defaults, then alphabetical.
    When nothing fits, the largest-capacity tool is chosen.
    """
    inf = float("inf")
    best: tuple[tuple, str] | None = None
    for name in registry.names():
        cap = _capabilities_of(registry, name)
        if cap is None:
            logging.getLogger("repgenr").warning(
                "auto-select skipping '%s' (%s): the plugin failed to load.",
                name, registry.group,
            )
            continue
        limit = cap.recommended_max_genomes
        limit_value = inf if limit is None else float(limit)
        fits = limit is None or limit >= n_items
        # Among fitting tools smaller limits sort first (tightest fit); when
        # nothing fits, larger capacity sorts first.
        key = (
            0 if _tool_available(cap) else 1,
            0 if fits else 1,
            limit_value if fits else -limit_value,
            _preference_rank(name),
            name,
        )
        if best is None or key < best[0]:
            best = (key, name)
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
