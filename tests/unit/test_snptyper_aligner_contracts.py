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
from repgenr.core.contracts import record_name
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
        f"> 1:1-12 + {ref}\n{_SEQ[record_name(ref)]}\n"
        f"> 2:1-12 + {query}\n{_SEQ[record_name(query)]}\n"
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
            # One full-length alignment against the reference (g1), so the
            # simple typer masks nothing; every genome here is 12 bases.
            query = record_name(cmd[-1])
            _write(
                Path(_flag_value(cmd, "-o") if "-o" in cmd else stdout_path),
                f"@SQ\tSN:g1\tLN:12\n{query}\t0\tg1\t1\t60\t12M\t*\t0\t0\t*\t*\n",
            )
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


def _make_fake_run_chain(fake_run_tool):
    """Adapters that run a chain of commands as one unit (the simple typer)."""

    def fake_run_chain(caps, steps, *, logger, **kwargs):
        for prefix, command in steps:
            fake_run_tool(caps, command, logger=logger, log_prefix=prefix)

    return fake_run_chain


@pytest.fixture()
def recorded(monkeypatch) -> list[list[str]]:
    calls: list[list[str]] = []
    fake = _make_fake_run_tool(calls)
    fake_chain = _make_fake_run_chain(fake)
    for reg in (snp_registry, align_registry):
        for name in reg.names():
            module = sys.modules[reg.get(name).__module__]
            monkeypatch.setattr(module, "run_tool", fake, raising=False)
            monkeypatch.setattr(module, "run_chain", fake_chain, raising=False)
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
    result = typer_.call(genomes, genomes[0], tmp_path / "snp_out", SnpParams(threads=7), _LOG)

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
    result = typer_.call(genomes, genomes[0], tmp_path / "snp_out", SnpParams(threads=7), _LOG)
    assert result.full_alignment is None


def test_snippy_names_the_reference_by_its_genome(genomes, tmp_path, monkeypatch) -> None:
    """snippy-core labels the reference 'Reference'; the typer renames it to the
    reference genome's stem so tree leaves match the genome names."""
    if "snippy" not in snp_registry.names():
        pytest.skip("snippy not registered")

    import repgenr.snptypers.snippy as snippy_mod

    ref_stem = genomes[0].stem
    msa = _CANNED_MSA.replace(f">{ref_stem}\n", ">Reference\n")
    assert ">Reference\n" in msa

    def fake_run_tool(caps, command, *, logger, stdout_path=None, cwd=None, **kwargs):
        cmd = [str(part) for part in command]
        tool = Path(cmd[0]).name
        if tool == "snippy":
            Path(_flag_value(cmd, "--outdir")).mkdir(parents=True, exist_ok=True)
        elif tool == "snippy-core":
            _write(Path(_flag_value(cmd, "--prefix") + ".aln"), msa)
            _write(Path(_flag_value(cmd, "--prefix") + ".full.aln"), msa)
        return 0

    monkeypatch.setattr(snippy_mod, "run_tool", fake_run_tool)
    typer_ = snp_registry.create("snippy")
    result = typer_.call(genomes, genomes[0], tmp_path / "snp_out", SnpParams(threads=7), _LOG)
    assert _read_headers(result.core_snp_fasta) == set(_STEMS)
    assert result.full_alignment is not None
    assert _read_headers(result.full_alignment) == set(_STEMS)


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
    "cactus": ["--reference", "g1", "--maxCores", "7"],
}


@pytest.mark.parametrize("tool", sorted(_ALIGN_PARAM_TOKENS))
def test_aligner_contract(tool, genomes, recorded, tmp_path) -> None:
    if tool not in align_registry.names():
        pytest.skip(f"{tool} not registered")
    aligner = align_registry.create(tool)
    result = aligner.align(
        genomes, genomes[0], tmp_path / "align_out", AlignParams(threads=7), _LOG
    )

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


@pytest.mark.parametrize("tool", sorted(_ALIGN_PARAM_TOKENS))
def test_aligner_names_gzipped_genomes_without_their_suffix(tool, recorded, tmp_path) -> None:
    """A genome x.fasta.gz is the record 'x' in the MSA, as in the typers and the tree.

    Path.stem gave 'x.fasta'. Only the file name matters for the record name,
    so the content is plain text here; aligners that cannot read gzip receive
    decompressed copies from the stage.
    """
    if tool not in align_registry.names():
        pytest.skip(f"{tool} not registered")
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    genomes = []
    for stem in _STEMS:
        path = gdir / f"{stem}.fasta.gz"
        path.write_text(f">{stem}\n{_SEQ[stem]}\n", encoding="utf-8")
        genomes.append(path)
    result = align_registry.create(tool).align(
        genomes, genomes[0], tmp_path / "align_out", AlignParams(threads=2), _LOG
    )
    assert _read_headers(result.msa_fasta) == set(_STEMS)


def test_every_registered_snptyper_and_aligner_has_contract_coverage() -> None:
    builtin_snp = {"simple", "snippy", "parsnp", "ska2"}
    builtin_aln = {"progressivemauve", "cactus", "sibeliaz"}
    assert builtin_snp & set(snp_registry.names()) <= set(_SNP_PARAM_TOKENS)
    assert builtin_aln & set(align_registry.names()) <= set(_ALIGN_PARAM_TOKENS)


@pytest.mark.parametrize("suffix", [".fasta", ".fasta.gz"])
def test_parsnp_names_records_by_genome_stem(suffix, tmp_path, monkeypatch) -> None:
    """harvesttools names records by file name, the reference with '.ref'; the
    typer renames them to genome stems, so tree leaves match the input genomes
    and tree2tax finds a versioned outgroup such as 'x_GCF_9.1'."""
    if "parsnp" not in snp_registry.names():
        pytest.skip("parsnp not registered")

    import repgenr.snptypers.parsnp as parsnp_mod

    gdir = tmp_path / "genomes"
    gdir.mkdir()
    stems = ["Fam_Gen_sp_GCF_1.1", "Fam_Gen_sp_GCF_2.1", "Fam_Gen_sp_GCF_9.1"]
    genomes = []
    for stem in stems:
        path = gdir / f"{stem}{suffix}"
        path.write_text(f">{stem}\nACGT\n", encoding="utf-8")
        genomes.append(path)
    # As observed with parsnp 2 and harvesttools 1.3 on the 50-genome set.
    harvest = (
        f">{stems[0]}{suffix}.ref\nACGT\n>{stems[1]}{suffix}\nACGA\n>{stems[2]}{suffix}\nACTT\n"
    )

    def fake_run_tool(caps, command, *, logger, stdout_path=None, cwd=None, **kwargs):
        cmd = [str(part) for part in command]
        if Path(cmd[0]).name == "parsnp":
            _write(Path(_flag_value(cmd, "-o")) / "parsnp.ggr", "GGR")
        elif "-S" in cmd:
            _write(Path(_flag_value(cmd, "-S")), harvest)
        elif "-M" in cmd:
            _write(Path(_flag_value(cmd, "-M")), harvest)
        return 0

    monkeypatch.setattr(parsnp_mod, "run_tool", fake_run_tool)
    result = snp_registry.create("parsnp").call(
        genomes, genomes[0], tmp_path / "snp_out", SnpParams(threads=2), _LOG
    )
    assert _read_headers(result.core_snp_fasta) == set(stems)
    assert result.full_alignment is not None
    assert _read_headers(result.full_alignment) == set(stems)
    assert sorted(p.name for p in (tmp_path / "snp_out").glob("harvest_*")) == []


def test_parsnp_hardlinks_its_query_genomes(genomes, recorded, tmp_path) -> None:
    """The query directory ParSNP reads holds links to the genomes, not copies."""
    if "parsnp" not in snp_registry.names():
        pytest.skip("parsnp not registered")
    snp_registry.create("parsnp").call(
        genomes, genomes[0], tmp_path / "snp_out", SnpParams(threads=2), _LOG
    )
    staged = tmp_path / "snp_out" / "input_genomes"
    assert sorted(p.name for p in staged.iterdir()) == [g.name for g in genomes[1:]]
    for genome in genomes[1:]:
        assert (staged / genome.name).stat().st_ino == genome.stat().st_ino


@pytest.mark.parametrize("suffix", [".fasta", ".fasta.gz"])
def test_cactus_names_msa_records_by_genome_stem(suffix, tmp_path, monkeypatch) -> None:
    """cactus sample names replace '.' with '_'; the MSA records are renamed back
    to genome stems so a versioned outgroup (x_GCF_9.1) stays findable by
    IQ-TREE's -o and by tree2tax."""
    if "cactus" not in align_registry.names():
        pytest.skip("cactus not registered")

    import repgenr.aligners.cactus as cactus_mod
    import repgenr.converters.hal_to_maf as h2m

    gdir = tmp_path / "genomes"
    gdir.mkdir()
    seqs = {"Fam_Gen_sp_GCF_1.1": "AAAACCCC", "Fam_Gen_sp_GCF_2.1": "AAATCCCC"}
    seqs["Fam_Gen_sp_GCF_9.1"] = "TAAACCCA"
    genomes = []
    for stem, seq in seqs.items():
        path = gdir / f"{stem}{suffix}"
        path.write_text(f">contig1\n{seq}\n", encoding="utf-8")
        genomes.append(path)
    maf_rows = "\n".join(
        f"s {stem.replace('.', '_')}.contig1 0 8 + 8 {seq}" for stem, seq in seqs.items()
    )

    def fake_run_tool(caps, command, *, logger, stdout_path=None, cwd=None, **kwargs):
        cmd = [str(part) for part in command]
        tool = Path(cmd[0]).name
        if tool == "cactus-pangenome":
            _write(Path(_flag_value(cmd, "--outDir")) / "pangenome.full.hal", "HAL")
        elif tool == "hal2maf":
            _write(Path(cmd[-1]), f"##maf version=1\na score=0\n{maf_rows}\n")
        return 0

    monkeypatch.setattr(cactus_mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(h2m, "run_tool", fake_run_tool)
    result = align_registry.create("cactus").align(
        genomes, genomes[0], tmp_path / "align", AlignParams(threads=2), _LOG
    )
    assert _read_headers(result.msa_fasta) == set(seqs)
    assert not (tmp_path / "align" / "cactus_samples.fasta").exists()
