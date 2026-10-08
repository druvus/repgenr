"""The Nextflow SKETCH module agrees with the Python sketch contract.

The module calls sourmash directly, so it declares the image and the sketch
parameters itself; these checks keep them equal to the values in
``repgenr.core.sourmash`` and the sketch contract.
"""

from __future__ import annotations

import re
from pathlib import Path

from repgenr.core.sourmash import SOURMASH_CONDA, SOURMASH_IMAGE

ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "nextflow" / "modules" / "local" / "dataflow" / "sketch.nf"
MODULES_CONFIG = ROOT / "nextflow" / "conf" / "modules.config"

# Literal on purpose: equals repgenr.core.sketches.SKETCH_PARAMS, which the
# sketch contract module defines.
SKETCH_PARAMS = "k=21,k=31,k=51,scaled=1000"


def test_module_container_is_the_shared_sourmash_image() -> None:
    text = MODULE.read_text(encoding="utf-8")
    assert re.findall(r"^\s*container\s+'([^']+)'", text, re.M) == [SOURMASH_IMAGE]
    assert re.findall(r"^\s*conda\s+'([^']+)'", text, re.M) == [" ".join(SOURMASH_CONDA)]


def test_module_sketches_with_the_contract_parameters() -> None:
    text = MODULE.read_text(encoding="utf-8")
    assert f"sourmash sketch dna -p {SKETCH_PARAMS} " in text
    assert ".partial.sig.zip" in text  # written under a temporary name, then renamed


def test_module_is_published_as_sketches() -> None:
    text = MODULES_CONFIG.read_text(encoding="utf-8")
    assert "withName: 'SKETCH'" in text
    assert "path('sketches')" in MODULE.read_text(encoding="utf-8")


def test_nextflow_files_are_ascii() -> None:
    for path in (MODULE, MODULES_CONFIG, ROOT / "nextflow" / "tests" / "sketch_process.nf.test"):
        path.read_text(encoding="utf-8").encode("ascii")
