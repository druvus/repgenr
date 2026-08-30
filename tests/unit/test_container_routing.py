"""Every external tool invocation must go through the container backend.

glance (dRep compare) and the two viral outgroup paths (mashtree) used the raw
process runner, so under ``--container`` those tools silently ran on the host.
These tests pin that the calls route through ``containers.run_tool`` with the
owning adapter's capabilities.
"""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

from repgenr.viral import _outgroup, bvbrc, selection
from repgenr.viral.bvbrc import _Record

_LOG = logging.getLogger("test")

_MATRIX = (
    "\tS_S1\tO_O1\n"
    "S_S1\t0\t0.5\n"
    "O_O1\t0.5\t0\n"
)


def _fake_run_tool(calls: list):
    def run_tool(caps, cmd, **kwargs):
        parts = [str(c) for c in cmd]
        calls.append((caps.name, parts))
        matrix = Path(parts[parts.index("--outmatrix") + 1])
        matrix.write_text(_MATRIX, encoding="utf-8")
        return 0

    return run_tool


def _seq(name: str) -> SimpleNamespace:
    return SimpleNamespace(description=name, seq="ACGT")


def test_selection_outgroup_mashtree_routed_through_run_tool(tmp_path, monkeypatch) -> None:
    from repgenr.viral.ncbi_virus import VirusRecord

    calls: list = []
    # both back-ends share the mashtree scaffolding in viral._outgroup
    monkeypatch.setattr(_outgroup, "run_tool", _fake_run_tool(calls))
    monkeypatch.setattr(_outgroup, "preflight", lambda caps: {"mashtree": "9.9"})

    def rec(acc: str, species: str) -> VirusRecord:
        return VirusRecord(
            accession=acc, taxid="1", organism=species, family="Fam", genus="Gen",
            species=species, length=100, completeness="complete", segment="", isolate="",
        )

    kept = [rec("S1", "SpeciesA")]
    records = kept + [rec("O1", "SpeciesB")]
    seqs = {"S1": _seq("S1"), "O1": _seq("O1")}
    ctx = SimpleNamespace(workdir=tmp_path, outgroup_dir=tmp_path / "outgroup")
    params = SimpleNamespace(outgroup_candidates_taxid_min_genomes=1, keep_files=False)

    og, versions = selection._determine_outgroup_records(
        ctx, records, kept, (90, 110), params, seqs, _LOG
    )
    assert og is not None and og.accession == "O1"
    assert versions == {"mashtree": "9.9"}
    assert [name for name, _ in calls] == ["mashtree"]
    assert calls[0][1][0] == "mashtree"


def test_bvbrc_outgroup_mashtree_routed_through_run_tool(tmp_path, monkeypatch) -> None:
    calls: list = []
    monkeypatch.setattr(_outgroup, "run_tool", _fake_run_tool(calls))
    monkeypatch.setattr(_outgroup, "preflight", lambda caps: {"mashtree": "9.9"})

    records = [
        _Record(name="S1", bvbrc_id="1.1", taxid="1", description="S1", length=100),
        _Record(name="O1", bvbrc_id="2.1", taxid="2", description="O1", length=100),
    ]
    sequences = {"S1": _seq("S1"), "O1": _seq("O1")}
    base = {
        "1": {"num": 5, "seq_med": 100},
        "2": {"num": 5, "seq_med": 100},
    }
    kept = {"1": {"1.1": 100}}
    ctx = SimpleNamespace(workdir=tmp_path, outgroup_dir=tmp_path / "outgroup")
    params = SimpleNamespace(outgroup_candidates_taxid_min_genomes=1, keep_files=False)

    versions = bvbrc._determine_outgroup(
        ctx, records, sequences, base, kept, (90, 110), params, _LOG
    )
    assert versions == {"mashtree": "9.9"}
    assert [name for name, _ in calls] == ["mashtree"]
    assert (tmp_path / "outgroup_accession.txt").read_text().strip() == "O1"


def test_divergent_length_tolerances_stay_divergent() -> None:
    """The two back-ends' tolerances differ on purpose (see viral/_outgroup.py);
    unifying them silently changes which genomes existing users get."""
    assert _outgroup.RECORDS_LENGTH_TOLERANCE == 0.15
    assert _outgroup.BVBRC_LENGTH_TOLERANCE == 0.10
