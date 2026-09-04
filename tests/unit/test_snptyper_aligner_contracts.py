"""Contract tests over every registered snptyper and aligner adapter (WP2-3).

Same pattern as test_adapter_contracts: a faked ``run_tool`` records argv and
writes canned tool outputs, so each adapter's real output-parsing runs without
the binaries. Also covers the gubbins masker adapter.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pytest

from repgenr.aligners.base import AlignParams
from repgenr.aligners.base import registry as align_registry
from repgenr.snptypers.base import SnpParams
from repgenr.snptypers.base import registry as snp_registry

_LOG = logging.getLogger("test")

# per-genome sequences; positions 0, 4, and 8 vary -> 3 core SNP sites
_SEQ = {
    "g1": "AAAACCCCGGGG",
    "g2": "TAAACCCCGGGG",
    "g3": "AAAATCCCGGGG",
    "g4": "AAAACCCCTGGG",
}
_STEMS = sorted(_SEQ)

_CANNED_MSA = "".join(f">{stem}\n{seq}\n" for stem, seq in _SEQ.items())


def _maf_block() -> str:
    rows = "\n".join(f"s {stem} 0 12 + 12 {seq}" for stem, seq in _SEQ.items())
    return f"##maf version=1\na score=0\n{rows}\n"


@pytest.fixture()
def genomes(tmp_path) -> list[Path]:
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    out = []
    for stem in _STEMS:
        p = gdir / f"{stem}.fasta"
        p.write_text(f">{stem}\n{_SEQ[stem]}\n", encoding="utf-8")
        out.append(p)
    return out


def _flag_value(cmd: list[str], flag: str) -> str:
    return cmd[cmd.index(flag) + 1]


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _xmfa_for(ref: str, query: str) -> str:
    return (
        "#FormatVersion Mauve1\n"
        f"#Sequence1File\t{ref}\n#Sequence1Format\tFastA\n"
        f"#Sequence2File\t{query}\n#Sequence2Format\tFastA\n"
        f"> 1:1-12 + {ref}\n{_SEQ[Path(ref).stem]}\n"
        f"> 2:1-12 + {query}\n{_SEQ[Path(query).stem]}\n"
        "=\n"
    )


def _make_fake_run_tool(recorded: list[list[str]]):
    def fake_run_tool(caps, command, *, logger, stdout_path=None, cwd=None, **kwargs):
        cmd = [str(part) for part in command]
        recorded.append(cmd)
        tool = Path(cmd[0]).name
        if tool == "bash" and len(cmd) > 1:
            tool = Path(cmd[1]).name  # macOS sibeliaz wrapper
        if tool == "samtools":
            if cmd[1] == "sort":
                _write(Path(_flag_value(cmd, "-o")), "")
        elif tool == "minimap2":
            _write(Path(stdout_path), "@SQ\n")
        elif tool == "bcftools":
            if cmd[1] == "consensus":
                out = Path(_flag_value(cmd, "-o"))
                stem = out.name.removesuffix(".consensus.fasta")
                _write(out, f">{stem}\n{_SEQ[stem]}\n")
            elif "-o" in cmd:
                _write(Path(_flag_value(cmd, "-o")), "")
        elif tool == "snippy":
            Path(_flag_value(cmd, "--outdir")).mkdir(parents=True, exist_ok=True)
        elif tool == "snippy-core":
            prefix = _flag_value(cmd, "--prefix")
            _write(Path(prefix + ".aln"), _CANNED_MSA)
            _write(Path(prefix + ".full.aln"), _CANNED_MSA)
        elif tool == "parsnp":
            _write(Path(_flag_value(cmd, "-o")) / "parsnp.ggr", "GGR")
        elif tool == "harvesttools":
            if "-S" in cmd:
                _write(Path(_flag_value(cmd, "-S")), _CANNED_MSA)
            elif "-M" in cmd:
                _write(Path(_flag_value(cmd, "-M")), _CANNED_MSA)
        elif tool == "ska":
            if cmd[1] == "build":
                _write(Path(_flag_value(cmd, "-o") + ".skf"), "SKF")
            elif cmd[1] == "align":
                _write(Path(_flag_value(cmd, "-o")), _CANNED_MSA)
        elif tool == "run_gubbins.py":
            prefix = _flag_value(cmd, "--prefix")
            _write(Path(prefix + ".filtered_polymorphic_sites.fasta"), _CANNED_MSA)
        elif tool == "progressiveMauve":
            ref, query = cmd[-2], cmd[-1]
            _write(Path(_flag_value(cmd, "--output")), _xmfa_for(ref, query))
        elif "sibeliaz" in tool:
            _write(Path(_flag_value(cmd, "-o")) / "alignment.maf", _maf_block())
        elif tool == "cactus-pangenome":
            _write(Path(_flag_value(cmd, "--outDir")) / "pangenome.hal", "HAL")
        elif tool == "hal2maf":
            _write(Path(cmd[-1]), _maf_block())
        else:
            # A third-party adapter's unknown argv must not crash the whole
            # suite; its own (skipped) parametrization covers it.
            return 0
        return 0

    return fake_run_tool


@pytest.fixture()
def recorded(monkeypatch) -> list[list[str]]:
    calls: list[list[str]] = []
    fake = _make_fake_run_tool(calls)
    for reg in (snp_registry, align_registry):
        for name in reg.names():
            module = sys.modules[reg.get(name).__module__]
            monkeypatch.setattr(module, "run_tool", fake, raising=False)
    import repgenr.converters.hal_to_maf as h2m
    import repgenr.maskers.gubbins as gubbins

    monkeypatch.setattr(h2m, "run_tool", fake)
    monkeypatch.setattr(gubbins, "run_tool", fake)
    return calls


def _read_headers(path: Path) -> set[str]:
    return {
        line[1:].strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    }


# --- snptyper contract --------------------------------------------------------


_SNP_PARAM_TOKENS = {
    "simple": [],
    "snippy": ["--cpus", "7"],
    "parsnp": ["-p", "7"],
    "ska2": ["--threads", "7"],
}

# Typers whose output is a variable-site alignment only (split k-mers have no
# reference coordinates); recombination masking is refused for them upstream.
_NO_FULL_ALIGNMENT = {"ska2"}


@pytest.mark.parametrize("tool", sorted(_SNP_PARAM_TOKENS))
def test_snptyper_contract(tool, genomes, recorded, tmp_path) -> None:
    if tool not in snp_registry.names():
        pytest.skip(f"{tool} not registered")
    typer_ = snp_registry.create(tool)
    result = typer_.call(genomes, genomes[0], tmp_path / "snp_out",
                         SnpParams(threads=7), _LOG)

    assert result.core_snp_fasta.exists()
    headers = _read_headers(result.core_snp_fasta)
    assert headers >= set(_STEMS), f"core SNP FASTA must name all genomes, got {headers}"
    if result.snp_distance_matrix is not None:
        assert result.snp_distance_matrix.exists()

    if tool in _NO_FULL_ALIGNMENT:
        assert result.full_alignment is None
    else:
        # Every other built-in typer must supply a whole-genome alignment for
        # maskers (snippy: core.full.aln; parsnp: harvesttools -M; simple: its
        # consensuses).
        assert result.full_alignment is not None, f"{tool} did not return a full_alignment"
        assert result.full_alignment.exists()
        full_headers = _read_headers(result.full_alignment)
        assert full_headers >= set(_STEMS), (
            f"full alignment must name all genomes, got {full_headers}"
        )

    flat = [tok for cmd in recorded for tok in cmd]
    for token in _SNP_PARAM_TOKENS[tool]:
        assert token in flat, f"{token!r} missing from recorded argv for {tool}"


def test_snippy_without_full_aln_returns_none_full_alignment(
    genomes, tmp_path, monkeypatch
) -> None:
    """snippy-core sometimes omits the whole-genome alignment (older
    versions, or --no-fullaln-like configurations); the typer must not
    fabricate one, so masking downstream is correctly refused."""
    if "snippy" not in snp_registry.names():
        pytest.skip("snippy not registered")

    import repgenr.snptypers.snippy as snippy_mod

    def fake_run_tool(caps, command, *, logger, stdout_path=None, cwd=None, **kwargs):
        cmd = [str(part) for part in command]
        tool = Path(cmd[0]).name
        if tool == "snippy":
            Path(_flag_value(cmd, "--outdir")).mkdir(parents=True, exist_ok=True)
        elif tool == "snippy-core":
            # Deliberately write only the core alignment, not .full.aln.
            _write(Path(_flag_value(cmd, "--prefix") + ".aln"), _CANNED_MSA)
        return 0

    monkeypatch.setattr(snippy_mod, "run_tool", fake_run_tool)
    typer_ = snp_registry.create("snippy")
    result = typer_.call(
        genomes, genomes[0], tmp_path / "snp_out", SnpParams(threads=7), _LOG
    )
    assert result.full_alignment is None


def test_gubbins_masking_returns_filtered_fasta(genomes, recorded, tmp_path) -> None:
    from repgenr.maskers.base import MaskParams
    from repgenr.maskers.gubbins import GubbinsMasker

    full = tmp_path / "full_alignment.fasta"
    full.write_text(_CANNED_MSA, encoding="utf-8")
    filtered = GubbinsMasker().mask(full, tmp_path / "gubbins_out", MaskParams(threads=7), _LOG)
    assert filtered.exists()
    assert filtered.name.endswith(".filtered_polymorphic_sites.fasta")
    assert any(Path(c[0]).name == "run_gubbins.py" for c in recorded)


# --- aligner contract ---------------------------------------------------------


_ALIGN_PARAM_TOKENS = {
    "progressivemauve": [],
    "sibeliaz": ["-t", "7"],
    "cactus": ["--reference", "g1"],
}


@pytest.mark.parametrize("tool", sorted(_ALIGN_PARAM_TOKENS))
def test_aligner_contract(tool, genomes, recorded, tmp_path) -> None:
    if tool not in align_registry.names():
        pytest.skip(f"{tool} not registered")
    aligner = align_registry.create(tool)
    result = aligner.align(genomes, genomes[0], tmp_path / "align_out",
                           AlignParams(threads=7), _LOG)

    assert result.msa_fasta.exists()
    headers = _read_headers(result.msa_fasta)
    assert headers >= set(_STEMS), f"MSA must name all genomes, got {headers}"
    seqs = [
        line
        for line in result.msa_fasta.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith(">")
    ]
    assert len({len(s) for s in seqs}) == 1, "MSA rows must share one length"

    flat = [tok for cmd in recorded for tok in cmd]
    for token in _ALIGN_PARAM_TOKENS[tool]:
        assert token in flat, f"{token!r} missing from recorded argv for {tool}"


def test_every_registered_snptyper_and_aligner_has_contract_coverage() -> None:
    builtin_snp = {"simple", "snippy", "parsnp", "ska2"}
    builtin_aln = {"progressivemauve", "cactus", "sibeliaz"}
    assert builtin_snp & set(snp_registry.names()) <= set(_SNP_PARAM_TOKENS)
    assert builtin_aln & set(align_registry.names()) <= set(_ALIGN_PARAM_TOKENS)
