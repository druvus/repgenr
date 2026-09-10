"""Builder-level validation added by the CLI audit (PR-C).

The CLI matrix test covers the command-line negative paths; these pin the
builder behaviour `run` and the per-stage commands share.
"""

from __future__ import annotations

import pytest

from repgenr.cli.param_builders import (
    metadata_params,
    phylo_params,
    require_mask,
    vgenome_params,
)
from repgenr.core.errors import UserInputError


def test_metadata_rejects_unknown_level_before_the_stage_runs() -> None:
    with pytest.raises(UserInputError, match="--level"):
        metadata_params(dataset="rep", level="kingdom")


def test_metadata_rejects_unknown_dataset() -> None:
    with pytest.raises(UserInputError, match="--dataset"):
        metadata_params(dataset="representative", level="genus")


def test_phylo_mask_under_aligner_source_is_an_error_not_a_silent_drop() -> None:
    with pytest.raises(UserInputError, match="--mask applies only with --msa-source snptype"):
        phylo_params(msa_source="aligner", extra={"mask": "gubbins"})


def test_phylo_mask_under_snptype_source_is_kept() -> None:
    params = phylo_params(msa_source="snptype", extra={"mask": "gubbins"})
    assert params.extra["mask"] == "gubbins"


def test_phylo_mask_must_be_a_registered_masker() -> None:
    with pytest.raises(UserInputError, match="--mask"):
        phylo_params(msa_source="snptype", extra={"mask": "clonalframeml"})


def test_require_mask_accepts_none_and_registered() -> None:
    require_mask("none")
    require_mask("gubbins")


def test_vgenome_signature_is_explicit() -> None:
    with pytest.raises(TypeError):
        vgenome_params(lenght_method="mean")  # a typo no longer reaches the dataclass


def test_vgenome_length_method_choice() -> None:
    assert vgenome_params(length_method="mean").length_method == "mean"
    with pytest.raises(UserInputError, match="--length-method"):
        vgenome_params(length_method="average")


def test_vgenome_outgroup_treebuilder_needs_distance_matrix_support() -> None:
    assert vgenome_params(outgroup_treebuilder="mashtree").outgroup_treebuilder == "mashtree"
    with pytest.raises(UserInputError, match="--outgroup-treebuilder"):
        vgenome_params(outgroup_treebuilder="iqtree")
