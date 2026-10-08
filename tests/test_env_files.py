"""The conda environment files under envs/ agree with the adapters.

Offline checks only: the solve itself runs in .github/workflows/envs.yml.
Every version floor an adapter enforces at preflight (BinarySpec.min_version)
must appear as the same ``>=`` floor in exactly one envs/*.yml file, so an
environment built from these files passes preflight.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import pytest
import yaml

from repgenr.core.plugins import ToolCapabilities

ENVS = Path(__file__).resolve().parents[1] / "envs"
EXPECTED_FILES = {"core", "gubbins", "mashtree", "snippy", "parsnp", "checkm2", "mauve"}
FAMILIES = (
    "aligners",
    "assemblers",
    "classifiers",
    "dereplicators",
    "maskers",
    "polishers",
    "snptypers",
    "treebuilders",
)
# Tools with a floor but no conda package: cactus is distributed as a container.
CONTAINER_ONLY = {"cactus-pangenome"}
# Packages allowed in more than one file. Empty: a package in two files means
# two environments may put different builds of one tool on PATH.
SHARED_ALLOWED: set[str] = set()

_SPEC_RE = re.compile(r"^([A-Za-z0-9_.\-]+)\s*(>=|==|=)?\s*([^\s#]+)?")


def _load(name: str) -> dict:
    return yaml.safe_load((ENVS / f"{name}.yml").read_text())


def _specs(name: str) -> dict[str, tuple[str | None, str | None]]:
    """Package name -> (operator, version) for the conda entries of one file."""
    out: dict[str, tuple[str | None, str | None]] = {}
    for dep in _load(name)["dependencies"]:
        if not isinstance(dep, str):  # the pip: sub-list
            continue
        match = _SPEC_RE.match(dep.strip())
        assert match, f"{name}.yml: cannot parse {dep!r}"
        out[match.group(1).lower()] = (match.group(2), match.group(3))
    return out


def _all_specs() -> dict[str, dict[str, tuple[str | None, str | None]]]:
    return {path.stem: _specs(path.stem) for path in sorted(ENVS.glob("*.yml"))}


def _capabilities() -> list[ToolCapabilities]:
    caps: list[ToolCapabilities] = []
    for family in FAMILIES:
        registry = importlib.import_module(f"repgenr.{family}.base").registry
        for name in registry.names():
            if not registry.is_broken(name):
                caps.append(registry.get(name).capabilities)
    from repgenr.stages.assemble_qc import CHECKM2_CAPS
    from repgenr.stages.genome import DATASETS_CAPS

    caps.extend([CHECKM2_CAPS, DATASETS_CAPS])
    return caps


def _package_for(binary: str, caps: ToolCapabilities) -> str:
    """The conda package that ships ``binary``, from the adapter's conda spec.

    A package whose name contains the binary name or is contained in it
    (gubbins / run_gubbins.py, snippy / snippy-core, ska2 / ska). A binary
    that matches none of the adapter's packages is its own package (skani
    under skder).
    """
    lowered = binary.lower()
    packages = [spec.split("::")[-1].split("=")[0].lower() for spec in caps.conda]
    for package in packages:
        if package == lowered:
            return package
    for package in packages:
        if package in lowered or lowered in package:
            return package
    return lowered


def test_the_expected_files_exist() -> None:
    assert {path.stem for path in ENVS.glob("*.yml")} == EXPECTED_FILES
    assert not (ENVS.parent / "environment.yml").exists()


@pytest.mark.parametrize("name", sorted(EXPECTED_FILES))
def test_each_file_names_its_environment_and_channels(name: str) -> None:
    data = _load(name)
    expected = "repgenr" if name == "core" else f"repgenr-{name}"
    assert data["name"] == expected
    assert data["channels"] == ["conda-forge", "bioconda"]
    first_line = (ENVS / f"{name}.yml").read_text().splitlines()[0]
    assert first_line.startswith("#") and "Platforms:" in first_line


def test_core_pins_python_312_and_installs_the_package() -> None:
    assert _specs("core")["python"] == ("=", "3.12")
    pip = [dep for dep in _load("core")["dependencies"] if isinstance(dep, dict)]
    assert pip == [{"pip": ["-e .."]}]


def test_no_package_is_in_two_files() -> None:
    seen: dict[str, str] = {}
    for env, specs in _all_specs().items():
        for package in specs:
            if package in SHARED_ALLOWED:
                continue
            assert package not in seen, f"{package} is in {seen[package]}.yml and {env}.yml"
            seen[package] = env


def test_every_preflight_floor_is_in_an_env_file() -> None:
    by_package = {pkg: spec for specs in _all_specs().values() for pkg, spec in specs.items()}
    checked = 0
    for caps in _capabilities():
        for spec in caps.required_binaries:
            if spec.min_version is None or spec.name in CONTAINER_ONLY:
                continue
            package = _package_for(spec.name, caps)
            assert package in by_package, f"{caps.name}: {package} is in no envs/*.yml file"
            assert by_package[package] == (">=", spec.min_version), (
                f"{caps.name}: {package} floor {by_package[package]} differs from "
                f"BinarySpec {spec.name} min_version {spec.min_version}"
            )
            checked += 1
    assert checked >= 20  # guards against the registries silently coming back empty


def test_package_mapping_follows_the_conda_spec() -> None:
    caps = ToolCapabilities(name="x", conda=("bioconda::gubbins",))
    assert _package_for("run_gubbins.py", caps) == "gubbins"
    assert _package_for("skani", ToolCapabilities(name="skder", conda=("bioconda::skder",))) == (
        "skani"
    )
    assert _package_for("ska", ToolCapabilities(name="ska2", conda=("bioconda::ska2",))) == "ska2"
