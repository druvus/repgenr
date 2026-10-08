"""stage_plain_inputs: decompressed copies for tools that cannot read gzip."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.core.plugins import ToolCapabilities
from repgenr.core.process import copy_plain_fasta, is_gzip, stage_plain_inputs

_LOG = logging.getLogger("test")
_PLAIN = ToolCapabilities(name="plain_only")
_GZ = ToolCapabilities(name="reads_gz", reads_gzip=True)
_CONTENT = b">c1\nACGTACGTAC\n>c2\nGGGTTT\n"


def _gz(path: Path, data: bytes = _CONTENT) -> Path:
    with gzip.open(path, "wb") as fh:
        fh.write(data)
    return path


def test_gzipped_input_yields_a_plain_copy_named_by_record_name(tmp_path: Path) -> None:
    src = _gz(tmp_path / "Fam_Gen_sp_GCF_1.1.fasta.gz")
    staged = stage_plain_inputs([src], _PLAIN, tmp_path / "staged", _LOG)
    dest = staged[src]
    assert dest == tmp_path / "staged" / "Fam_Gen_sp_GCF_1.1.fasta"
    assert dest.read_bytes() == _CONTENT
    assert not list((tmp_path / "staged").glob("*.tmp"))


def test_plain_input_passes_through(tmp_path: Path) -> None:
    src = tmp_path / "x.fna"
    src.write_bytes(_CONTENT)
    staged = stage_plain_inputs([src], _PLAIN, tmp_path / "staged", _LOG)
    assert staged == {src: src}
    assert not (tmp_path / "staged").exists(), "nothing to decompress: no directory"


def test_gzip_is_detected_by_magic_bytes_not_by_name(tmp_path: Path) -> None:
    src = _gz(tmp_path / "x.fasta")  # gzipped content under a plain name
    plain_gz_name = tmp_path / "y.fasta.gz"
    plain_gz_name.write_bytes(_CONTENT)  # plain content under a .gz name
    assert is_gzip(src) and not is_gzip(plain_gz_name)
    staged = stage_plain_inputs([src, plain_gz_name], _PLAIN, tmp_path / "staged", _LOG)
    assert staged[src] == tmp_path / "staged" / "x.fasta"
    assert staged[src].read_bytes() == _CONTENT
    assert staged[plain_gz_name] == plain_gz_name


def test_tool_that_reads_gzip_gets_the_inputs_unchanged(tmp_path: Path) -> None:
    src = _gz(tmp_path / "x.fasta.gz")
    staged = stage_plain_inputs([src], _GZ, tmp_path / "staged", _LOG)
    assert staged == {src: src}
    assert not (tmp_path / "staged").exists()


@pytest.mark.parametrize("caps", [_PLAIN, _GZ], ids=["no_gzip", "reads_gzip"])
@pytest.mark.parametrize("other", ["x.fasta", "x.fna"])
def test_two_inputs_with_one_record_name_are_refused(
    tmp_path: Path, caps: ToolCapabilities, other: str
) -> None:
    """record_name gives x for x.fasta and x.fasta.gz alike; their records and
    leaves could not be told apart, whichever tool runs."""
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    one = _gz(tmp_path / "a" / "x.fasta.gz")
    two = tmp_path / "b" / other
    two.write_bytes(_CONTENT)
    with pytest.raises(UserInputError, match="share the record name 'x'"):
        stage_plain_inputs([one, two], caps, tmp_path / "staged", _LOG)
    assert not (tmp_path / "staged").exists()


def test_one_path_given_twice_is_staged_once(tmp_path: Path) -> None:
    """The reference is usually one of the genomes as well."""
    src = _gz(tmp_path / "x.fasta.gz")
    staged = stage_plain_inputs([src, src], _PLAIN, tmp_path / "staged", _LOG)
    assert list(staged) == [src]


def test_copy_plain_fasta_decompresses_and_copies(tmp_path: Path) -> None:
    _gz(tmp_path / "in.fasta.gz")
    copy_plain_fasta(tmp_path / "in.fasta.gz", tmp_path / "out1.fasta")
    (tmp_path / "plain.fasta").write_bytes(_CONTENT)
    copy_plain_fasta(tmp_path / "plain.fasta", tmp_path / "out2.fasta")
    assert (tmp_path / "out1.fasta").read_bytes() == _CONTENT
    assert (tmp_path / "out2.fasta").read_bytes() == _CONTENT


def test_reads_gzip_flags_of_the_builtin_adapters() -> None:
    """Measured on the pinned versions (verification.md); False is the safe default."""
    from repgenr.aligners.base import registry as aligners
    from repgenr.snptypers.base import registry as typers
    from repgenr.treebuilders.base import registry as builders

    flags = {
        **{n: aligners.get(n).capabilities.reads_gzip for n in aligners.names()},
        **{f"snp:{n}": typers.get(n).capabilities.reads_gzip for n in typers.names()},
        **{f"tree:{n}": builders.get(n).capabilities.reads_gzip for n in builders.names()},
    }
    expected_true = {"cactus", "snp:ska2", "snp:simple", "tree:sourmash", "tree:mashtree"}
    expected_false = {"progressivemauve", "sibeliaz", "snp:parsnp", "snp:snippy"}
    assert {k for k in expected_true if flags.get(k)} == {k for k in expected_true if k in flags}
    assert not any(flags.get(k) for k in expected_false)
