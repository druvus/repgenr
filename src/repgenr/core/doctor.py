"""Workdir health checks: verify outputs against records (``repgenr doctor``).

``repgenr status`` reports what ``repgenr.yaml`` claims; ``doctor`` verifies the
claims against the filesystem and the manifest -- interrupted stages, missing or
corrupt genomes, manifest drift, representative/cluster mismatches, truncated
deliverables, unresolvable outgroups, leftover temp files, and stages whose
recorded input digests no longer match reality (they will re-run).

Strictly read-only: no file, log, or manifest is created in the workdir.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from .config import CONFIG_FILENAME, Config
from .contracts import (
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    GENOMES_MAP_TSV,
    SELECTION_TSV,
    TREE2TAX_TSV,
    TREE_NWK,
    accession_from_filename,
    list_fasta,
    read_clusters,
    read_selection,
)
from .inputs import inputs_digest, manifest_digest
from .integrity import (
    check_genome_completeness,
    check_representatives_consistency,
    looks_like_fasta,
)
from .manifest import MANIFEST_FILENAME, Manifest

_LOG = logging.getLogger(__name__)
_MAX_LISTED = 3  # examples shown per finding


@dataclass
class Finding:
    level: str  # "ok" | "warn" | "fail"
    area: str
    message: str


def diagnose(workdir: Path) -> list[Finding]:
    """Run every health check; never raises for a broken workdir."""
    workdir = Path(workdir)
    if not (workdir / CONFIG_FILENAME).exists():
        return [Finding("warn", "config", f"No RepGenR run found at {workdir}.")]

    findings: list[Finding] = []
    config = Config.load(workdir)
    checks = (
        _check_stage_records,
        _check_genomes,
        _check_manifest_drift,
        _check_outgroup,
        _check_representatives,
        _check_tree,
        _check_tree2tax_pair,
        _check_stale_inputs,
        _check_leftovers,
    )
    for check in checks:
        try:
            findings.extend(check(workdir, config))
        except Exception as exc:  # a broken artifact must not abort the report
            findings.append(
                Finding("fail", check.__name__.removeprefix("_check_"),
                        f"Check could not complete: {exc}")
            )
    if not any(f.level == "fail" for f in findings):
        findings.append(Finding("ok", "summary", "No failures detected."))
    return findings


def _examples(names: list[str]) -> str:
    shown = ", ".join(names[:_MAX_LISTED])
    more = f" (+{len(names) - _MAX_LISTED} more)" if len(names) > _MAX_LISTED else ""
    return shown + more


def _check_stage_records(workdir: Path, config: Config) -> list[Finding]:
    out: list[Finding] = []
    for name, record in config.stages.items():
        if record.completed:
            out.append(Finding("ok", name, f"completed {record.completed}"))
        elif record.params or record.tool:
            out.append(Finding(
                "fail", name,
                "started a run and never finished (crashed mid-run); its outputs "
                "may be partial -- re-run the stage.",
            ))
    return out


def _check_genomes(workdir: Path, config: Config) -> list[Finding]:
    genomes_dir = workdir / "genomes"
    out: list[Finding] = []
    shortfall = check_genome_completeness(
        genomes_dir, workdir, logger=_LOG, allow_incomplete=True
    )
    if shortfall:
        out.append(Finding(
            "fail", "genomes",
            f"{len(shortfall)} selected genome(s) missing from {genomes_dir} "
            f"(e.g. {_examples(shortfall)}); re-run the genome stage.",
        ))
    bad = [p.name for p in list_fasta(genomes_dir) if not looks_like_fasta(p)]
    if bad:
        out.append(Finding(
            "fail", "genomes",
            f"{len(bad)} file(s) under {genomes_dir} are not FASTA "
            f"(e.g. {_examples(bad)}); delete them and re-run the genome stage.",
        ))
    if not shortfall and not bad and genomes_dir.exists():
        out.append(Finding("ok", "genomes",
                           f"{len(list_fasta(genomes_dir))} genome file(s) look sound"))
    return out


def _check_manifest_drift(workdir: Path, config: Config) -> list[Finding]:
    selection = workdir / SELECTION_TSV
    manifest_path = workdir / MANIFEST_FILENAME
    if not selection.exists() or not manifest_path.exists():
        return []
    selected = {row.accession for row in read_selection(selection)}
    manifest = Manifest.open_readonly(manifest_path)
    try:
        recorded = {g.accession for g in manifest.all_genomes(include_outgroup=True)}
    finally:
        manifest.close()
    extra = sorted(recorded - selected)
    missing = sorted(selected - recorded)
    out: list[Finding] = []
    if extra:
        out.append(Finding(
            "fail", "manifest",
            f"manifest holds {len(extra)} genome(s) not in selection.tsv "
            f"(e.g. {_examples(extra)}); re-run the metadata stage to reconcile.",
        ))
    if missing:
        out.append(Finding(
            "fail", "manifest",
            f"manifest is missing {len(missing)} selected genome(s) "
            f"(e.g. {_examples(missing)}); re-run the metadata stage.",
        ))
    if not extra and not missing:
        out.append(Finding("ok", "manifest", "manifest matches selection.tsv"))
    return out


def _check_outgroup(workdir: Path, config: Config) -> list[Finding]:
    acc_file = workdir / "outgroup_accession.txt"
    if not acc_file.exists():
        return []
    accession = acc_file.read_text(encoding="utf-8").strip()
    if not accession:
        return []
    outgroup_dir = workdir / "outgroup"
    candidates = (
        [p for p in sorted(outgroup_dir.iterdir()) if not p.name.startswith(".")]
        if outgroup_dir.exists() else []
    )
    resolved = any(
        accession_from_filename(f.name) == accession or accession in f.name
        for f in candidates
    )
    if resolved:
        return [Finding("ok", "outgroup", f"outgroup {accession} resolves")]
    return [Finding(
        "warn", "outgroup",
        f"recorded outgroup {accession} has no matching file under {outgroup_dir}; "
        "phylo/tree2tax will proceed unrooted.",
    )]


def _check_representatives(workdir: Path, config: Config) -> list[Finding]:
    derep_dir = workdir / "derep"
    clusters = derep_dir / CLUSTERS_TSV
    if not clusters.exists():
        return []
    out: list[Finding] = []
    shortfall = check_representatives_consistency(
        derep_dir / "representatives", clusters, logger=_LOG, allow_incomplete=True
    )
    if shortfall:
        out.append(Finding(
            "fail", "dereplicate",
            f"{len(shortfall)} representative(s) listed in {CLUSTERS_TSV} are "
            f"absent on disk (e.g. {_examples(shortfall)}); re-run dereplicate.",
        ))
    extras = sorted(
        {p.name for p in list_fasta(derep_dir / "representatives")}
        - set(read_clusters(clusters))
    )
    if extras:
        out.append(Finding(
            "warn", "dereplicate",
            f"{len(extras)} file(s) under representatives/ are not in "
            f"{CLUSTERS_TSV} (e.g. {_examples(extras)}).",
        ))
    if not (derep_dir / GENOME_STATUS_TSV).exists():
        out.append(Finding(
            "fail", "dereplicate",
            f"{CLUSTERS_TSV} present but {GENOME_STATUS_TSV} missing; the "
            "dereplicate stage likely crashed mid-write -- re-run it.",
        ))
    if not out:
        out.append(Finding("ok", "dereplicate", "representatives match clusters.tsv"))
    return out


def _check_tree(workdir: Path, config: Config) -> list[Finding]:
    tree = workdir / "tree" / TREE_NWK
    if not tree.exists():
        return []
    content = tree.read_text(encoding="utf-8").strip()
    if not content or not content.endswith(";"):
        return [Finding(
            "fail", "phylo",
            f"{tree} is empty or truncated (no terminating ';'); re-run phylo.",
        )]
    return [Finding("ok", "phylo", f"{TREE_NWK} looks like a complete tree")]


def _check_tree2tax_pair(workdir: Path, config: Config) -> list[Finding]:
    t2t = workdir / TREE2TAX_TSV
    gmap = workdir / GENOMES_MAP_TSV
    if t2t.exists() == gmap.exists():
        return (
            [Finding("ok", "tree2tax", "deliverable pair present")]
            if t2t.exists() else []
        )
    missing = GENOMES_MAP_TSV if t2t.exists() else TREE2TAX_TSV
    return [Finding(
        "fail", "tree2tax",
        f"deliverables are mismatched: {missing} is missing while its partner "
        "exists; the tree2tax stage likely crashed mid-write -- re-run it.",
    )]


def _check_stale_inputs(workdir: Path, config: Config) -> list[Finding]:
    """Completed stages whose recorded input digests no longer match reality."""
    from ..cli.base import _MANIFEST_INPUT_STAGES, STAGE_INPUTS  # deferred: core<-cli

    ctx = SimpleNamespace(
        workdir=workdir,
        genomes_dir=workdir / "genomes",
        outgroup_dir=workdir / "outgroup",
        derep_dir=workdir / "derep",
        representatives_dir=workdir / "derep" / "representatives",
        snp_dir=workdir / "snp",
        tree_dir=workdir / "tree",
    )
    out: list[Finding] = []
    for name, record in config.stages.items():
        if not record.completed or not record.inputs:
            continue
        spec = STAGE_INPUTS.get(name)
        if spec is None:
            continue
        params = SimpleNamespace(**record.params)
        digests = inputs_digest(workdir, spec(ctx, params))
        if name in _MANIFEST_INPUT_STAGES:
            manifest_path = workdir / MANIFEST_FILENAME
            if manifest_path.exists():
                manifest = Manifest.open_readonly(manifest_path)
                try:
                    digests["manifest"] = manifest_digest(manifest)
                finally:
                    manifest.close()
        changed = sorted(
            key for key in {*record.inputs, *digests}
            if record.inputs.get(key) != digests.get(key)
        )
        if changed:
            out.append(Finding(
                "warn", name,
                f"input(s) changed since completion ({_examples(changed)}); "
                "the stage will re-run on its next invocation.",
            ))
    return out


def _check_leftovers(workdir: Path, config: Config) -> list[Finding]:
    leftovers = sorted(
        str(p.relative_to(workdir))
        for pattern in ("*.tmp", "*.part")
        for p in workdir.rglob(pattern)
        if "scratch" not in p.parts
    )
    if not leftovers:
        return []
    return [Finding(
        "warn", "leftovers",
        f"{len(leftovers)} temp file(s) from an interrupted write "
        f"(e.g. {_examples(leftovers)}); safe to delete.",
    )]
