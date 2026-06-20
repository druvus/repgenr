"""RepGenR command-line interface.

Replaces the old ``repgenr.py`` string-rewriting dispatcher with a real Typer
app. Each subcommand parses arguments, builds a :class:`WorkdirContext`, and
calls the matching ``stages.<name>.run(ctx, params)``. Errors of type
:class:`RepGenRError` are logged cleanly and turned into a non-zero exit code.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import logging
import os
import sys
from pathlib import Path
from typing import Any

import typer

from .. import __version__
from ..core.context import WorkdirContext
from ..core.errors import RepGenRError, UserInputError
from ..core.logging import configure_logging

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="RepGenR: modular genome dereplication, alignment, SNP typing and phylogenetics.",
)

# Top-level run options shared by every subcommand (set in the callback).
_RUN_STATE: dict[str, Any] = {"force": False, "log_level": logging.INFO}


def _require_choice(value: str, choices: set[str], label: str) -> None:
    if value not in choices:
        raise UserInputError(
            f"Invalid {label} {value!r}. Choose from: {', '.join(sorted(choices))}."
        )


def _require_unit_interval(value: float | None, label: str) -> None:
    if value is not None and not (0.0 < value <= 1.0):
        raise UserInputError(f"{label} must be in (0, 1], got {value}.")


def _stage_fingerprint(stage_name: str, params: object) -> str:
    """Stable hash of a stage invocation, used to skip already-completed work.

    Built from the stage name plus the parameter object (a dataclass), so the
    same invocation produces the same fingerprint across runs. Paths and other
    non-JSON values are stringified.
    """
    if dataclasses.is_dataclass(params) and not isinstance(params, type):
        payload: object = dataclasses.asdict(params)
    else:
        payload = vars(params)
    blob = json.dumps({"stage": stage_name, "params": payload}, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"repgenr {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False, "--version", callback=_version_callback, is_eager=True,
        help="Show version and exit.",
    ),
    container: str = typer.Option(
        "none", "--container", envvar="REPGENR_CONTAINER",
        help="Run external tools in containers: none, docker, or singularity.",
    ),
    container_engine: str | None = typer.Option(
        None, "--container-engine", envvar="REPGENR_CONTAINER_ENGINE",
        help="Engine binary override (e.g. apptainer, podman).",
    ),
    container_cache: str | None = typer.Option(
        None, "--container-cache", envvar="REPGENR_CONTAINER_CACHE",
        help="Directory for Singularity .sif images / Wave cache (large; can be external).",
    ),
    platform: str | None = typer.Option(
        None, "--platform", envvar="REPGENR_CONTAINER_PLATFORM",
        help="Container platform, e.g. linux/amd64 for emulated BioContainers on arm64.",
    ),
    wave: bool = typer.Option(
        False, "--wave/--no-wave", envvar="REPGENR_WAVE",
        help="Resolve images for multi-tool adapters via the Seqera Wave CLI.",
    ),
    force: bool = typer.Option(
        False, "--force/--no-force", "-f", envvar="REPGENR_FORCE",
        help="Re-run a stage even if it already completed with the same parameters.",
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose (DEBUG) logging."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Only warnings and errors."),
) -> None:
    """RepGenR top-level entry point."""
    from ..core.containers import configure_container

    _RUN_STATE["force"] = force
    if quiet:
        level = logging.WARNING
    elif verbose:
        level = logging.DEBUG
    else:
        env = os.environ.get("REPGENR_LOG_LEVEL")
        level = getattr(logging, env.upper(), logging.INFO) if env else logging.INFO
    _RUN_STATE["log_level"] = level
    configure_container(
        backend=container, engine=container_engine, platform=platform,
        cache_dir=container_cache, wave_enabled=wave,
    )


def _run(stage_name: str, workdir: Path, build_params, *, create: bool = False) -> None:
    """Common harness: context, dispatch, clean error handling.

    Resume: a stage that already completed with the same parameters is skipped
    (fingerprint match), unless ``--force`` is set. A stage that crashed before
    recording completion has no ``completed`` stamp and so always re-runs.
    """
    logger = configure_logging(
        workdir if (create or workdir.exists()) else None, level=_RUN_STATE["log_level"]
    )
    try:
        ctx = WorkdirContext(workdir, logger=logger, create=create)
        params = build_params()
        fingerprint = _stage_fingerprint(stage_name, params)
        prior = ctx.config.stages.get(stage_name)
        if (
            not _RUN_STATE["force"]
            and prior is not None
            and prior.completed
            and prior.fingerprint == fingerprint
        ):
            logger.info(
                "Stage '%s' already completed with the same parameters; skipping "
                "(use --force to re-run).", stage_name,
            )
            return
        module = __import__(f"repgenr.stages.{stage_name}", fromlist=["run"])
        module.run(ctx, params)
        # Stamp the fingerprint on the record the stage just wrote, so the next
        # invocation with identical params can skip.
        record = ctx.config.stages.get(stage_name)
        if record is not None:
            record.fingerprint = fingerprint
            ctx.save_config()
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


@app.command()
def metadata(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    dataset: str = typer.Option(..., "-d", "--dataset", help="all or rep."),
    level: str = typer.Option(..., "-l", "--level", help="family, genus or species."),
    source: str = typer.Option(
        "tsv", "--source", help="tsv (download full table) or api (GTDB API, target only)."
    ),
    release: str | None = typer.Option(None, "-r", "--release", help="GTDB release (tsv source)."),
    version: str | None = typer.Option(None, "-v", "--version", help="bac120/ar53 (tsv source)."),
    target_family: str | None = typer.Option(None, "-tf", "--target-family"),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus"),
    target_species: str | None = typer.Option(None, "-ts", "--target-species"),
    outgroup_accession: str | None = typer.Option(None, "--outgroup-accession"),
    metadata_path: str | None = typer.Option(None, "--metadata-path"),
    nodownload: bool = typer.Option(False, "--nodownload"),
    limit: int | None = typer.Option(None, "--limit"),
) -> None:
    """Select a taxon's genomes from GTDB (full table or the GTDB API)."""
    from ..stages.metadata import MetadataParams

    def build() -> MetadataParams:
        return MetadataParams(
            dataset=dataset, level=level, source=source,
            release=release, version=version,
            target_family=target_family, target_genus=target_genus,
            target_species=target_species, outgroup_accession=outgroup_accession,
            metadata_path=metadata_path, nodownload=nodownload, limit=limit,
        )

    _run("metadata", workdir, build, create=True)


@app.command()
def vmetadata(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory (created)."),
    target: str | None = typer.Option(None, "-t", "--target", help="Virus taxon/group/family."),
    source: str = typer.Option(
        "ncbi_virus", "--source", help="ncbi_virus (NCBI Virus via datasets) or bvbrc."
    ),
    filter: str = typer.Option("complete genome", "-f", "--filter", help="BV-BRC header tag."),
    host: str | None = typer.Option(None, "--host", help="ncbi_virus: restrict to a host species."),
    complete_only: bool = typer.Option(
        False, "--complete-only", help="ncbi_virus: only COMPLETE sequences."
    ),
    released_after: str | None = typer.Option(
        None, "--released-after", help="ncbi_virus: MM/DD/YYYY."
    ),
    list_targets: bool = typer.Option(False, "-l", "--list", help="List BV-BRC targets and exit."),
) -> None:
    """Retrieve viral metadata from NCBI Virus (default) or BV-BRC."""
    from ..stages.vmetadata import VmetadataParams

    def build() -> VmetadataParams:
        _require_choice(source, {"ncbi_virus", "bvbrc"}, "--source")
        return VmetadataParams(
            target=target, filter=filter, list_targets=list_targets,
            source=source, host=host, complete_only=complete_only,
            released_after=released_after,
        )

    _run("vmetadata", workdir, build, create=True)


@app.command()
def vgenome(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus"),
    target_species: str | None = typer.Option(None, "-ts", "--target-species"),
    target_serotype: str | None = typer.Option(None, "-tse", "--target-serotype"),
    target_custom: str | None = typer.Option(None, "-tc", "--target-custom", help="key:value."),
    length_all: bool = typer.Option(False, "--length-all"),
    length_deviation: int = typer.Option(10, "--length-deviation"),
    length_method: str = typer.Option("median_of_medians", "--length-method"),
    length_range: str | None = typer.Option(None, "--length-range", help="e.g. 25000-35000."),
    discard: str | None = typer.Option(None, "--discard", help="Comma-separated header tags."),
    no_outgroup: bool = typer.Option(False, "--no-outgroup"),
    min_outgroup_genomes: int = typer.Option(5, "--outgroup-candidates-taxid-min-genomes"),
    glance: bool = typer.Option(False, "--glance", help="Print selection and stop."),
    print_fasta_headers: bool = typer.Option(False, "--print-fasta-headers"),
    ignore_duplicates: bool = typer.Option(False, "--ignore-duplicates"),
    keep_files: bool = typer.Option(False, "--keep-files"),
) -> None:
    """Select and organize viral genomes (virus equivalent of genome)."""
    from ..stages.vgenome import VgenomeParams

    def build() -> VgenomeParams:
        return VgenomeParams(
            target_genus=target_genus, target_species=target_species,
            target_serotype=target_serotype, target_custom=target_custom,
            length_all=length_all, length_deviation=length_deviation,
            length_method=length_method, length_range=length_range, discard=discard,
            no_outgroup=no_outgroup,
            outgroup_candidates_taxid_min_genomes=min_outgroup_genomes,
            glance=glance, print_fasta_headers=print_fasta_headers,
            ignore_duplicates=ignore_duplicates, keep_files=keep_files,
        )

    _run("vgenome", workdir, build)


@app.command()
def genome(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    accession_list_only: bool = typer.Option(False, "--accession-list-only"),
    keep_files: bool = typer.Option(False, "--keep-files"),
) -> None:
    """Download and organize genomes selected by the metadata stage."""
    from ..stages.genome import GenomeParams

    def build() -> GenomeParams:
        return GenomeParams(accession_list_only=accession_list_only, keep_files=keep_files)

    _run("genome", workdir, build)


@app.command()
def dereplicate(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    tool: str = typer.Option("skder", "--tool", help="auto/skder/drep/galah/sourmash."),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "-af", "--aligned-fraction"),
    threads: int = typer.Option(16, "-t", "--threads"),
    process_size: int | None = typer.Option(
        None, "-s", "--process-size",
        help="Chunk size; when set and exceeded, two-stage chunking runs for any tool.",
    ),
    num_processes: int = typer.Option(
        1, "-p", "--num-processes", help="Parallel stage-1 chunk workers (threads split across)."
    ),
    pre_primary_ani: float | None = typer.Option(
        None, "--pre-primary-ani",
        help="Stage-1 (intra-chunk) primary ANI; defaults to --primary-ani.",
    ),
    pre_secondary_ani: float | None = typer.Option(
        None, "--pre-secondary-ani",
        help="Stage-1 (intra-chunk) secondary ANI; defaults to --secondary-ani.",
    ),
    reduce: str = typer.Option(
        "none", "--reduce",
        help="Taxonomy-aware reduction after ANI: none, species, or genus "
        "(one representative per taxon).",
    ),
    target_reps: int = typer.Option(
        0, "--target-reps",
        help="Target representative count: search --secondary-ani to land near it "
        "(0 = off; re-runs dereplication per search step).",
    ),
    virus: bool = typer.Option(False, "--virus", help="Pass virus-tuned parameters to the tool."),
) -> None:
    """Cluster genomes by ANI and select representatives."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.dereplicate import DereplicateParams

    def build() -> DereplicateParams:
        _require_choice(tool, {"auto", *_derep_registry.names()}, "--tool")
        _require_choice(reduce, {"none", "species", "genus"}, "--reduce")
        if target_reps < 0:
            raise UserInputError(f"--target-reps must be >= 0, got {target_reps}.")
        _require_unit_interval(primary_ani, "--primary-ani")
        _require_unit_interval(secondary_ani, "--secondary-ani")
        _require_unit_interval(aligned_fraction, "--aligned-fraction")
        _require_unit_interval(pre_primary_ani, "--pre-primary-ani")
        _require_unit_interval(pre_secondary_ani, "--pre-secondary-ani")
        return DereplicateParams(
            tool=tool,
            primary_ani=primary_ani,
            secondary_ani=secondary_ani,
            aligned_fraction=aligned_fraction,
            threads=threads,
            process_size=process_size,
            num_processes=num_processes,
            pre_primary_ani=pre_primary_ani,
            pre_secondary_ani=pre_secondary_ani,
            reduce=reduce,
            target_reps=target_reps,
            extra={"virus": virus} if virus else {},
        )

    _run("dereplicate", workdir, build)


@app.command()
def snptype(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    tool: str = typer.Option("simple", "--tool", help="SNP typer: simple/snippy/parsnp."),
    reference: str | None = typer.Option(None, "--reference", help="Reference genome filename."),
    all_genomes: bool = typer.Option(False, "--all-genomes", help="Use all genomes, not reps."),
    mask: str = typer.Option("none", "--mask", help="Recombination masking: none or gubbins."),
    threads: int = typer.Option(16, "-t", "--threads"),
) -> None:
    """Call SNPs and build a core-SNP alignment."""
    from ..snptypers.base import registry as _snp_registry
    from ..stages.snptype import SnptypeParams

    def build() -> SnptypeParams:
        _require_choice(tool, set(_snp_registry.names()), "--tool")
        _require_choice(mask, {"none", "gubbins"}, "--mask")
        return SnptypeParams(
            tool=tool,
            threads=threads,
            reference=reference,
            all_genomes=all_genomes,
            mask=mask,
        )

    _run("snptype", workdir, build)


@app.command()
def phylo(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    treebuilder: str = typer.Option(
        "iqtree", "--treebuilder", help="auto/iqtree/fasttree/raxmlng/mashtree/sourmash."
    ),
    msa_source: str = typer.Option("aligner", "--msa-source", help="aligner or snptype."),
    aligner: str = typer.Option(
        "progressivemauve", "--aligner", help="progressivemauve, cactus, sibeliaz."
    ),
    snptyper: str = typer.Option("simple", "--snptyper", help="SNP typer for snptype source."),
    all_genomes: bool = typer.Option(False, "--all-genomes", help="Use all genomes, not reps."),
    no_outgroup: bool = typer.Option(False, "--no-outgroup", help="Do not root with an outgroup."),
    bootstrap: int = typer.Option(0, "-B", "--bootstrap", help="Bootstrap replicates (>=1000)."),
    reference: str | None = typer.Option(None, "--reference", help="Reference genome filename."),
    aligner_arg: list[str] = typer.Option(
        [], "--aligner-arg",
        help="Aligner tuning as key=value (repeatable), e.g. kmer=15 (sibeliaz) "
        "or seed_weight=11 (progressivemauve).",
    ),
    threads: int = typer.Option(16, "-t", "--threads"),
) -> None:
    """Build a phylogenetic tree from an alignment, SNP alignment, or directly."""
    from ..aligners.base import registry as _aln_registry
    from ..snptypers.base import registry as _snp_registry
    from ..stages.phylo import PhyloParams
    from ..treebuilders.base import registry as _tb_registry

    def build() -> PhyloParams:
        _require_choice(treebuilder, {"auto", *_tb_registry.names()}, "--treebuilder")
        _require_choice(msa_source, {"aligner", "snptype"}, "--msa-source")
        if msa_source == "aligner":
            _require_choice(aligner, set(_aln_registry.names()), "--aligner")
        else:
            _require_choice(snptyper, set(_snp_registry.names()), "--snptyper")
        return PhyloParams(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            all_genomes=all_genomes,
            no_outgroup=no_outgroup,
            bootstrap=bootstrap,
            reference=reference,
            threads=threads,
            extra=_parse_key_values(aligner_arg, "--aligner-arg"),
        )

    _run("phylo", workdir, build)


@app.command()
def tree2tax(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    node_basename: str | None = typer.Option(None, "--node-basename", help="Prefix for nodes."),
    root_name: str = typer.Option("root", "-r", "--root-name", help="Name for the root node."),
    remove_outgroup: bool = typer.Option(False, "--remove-outgroup", help="Drop outgroup."),
    include_dereplicated: bool = typer.Option(
        False, "--include-dereplicated", help="List redundant genomes under their representative."
    ),
) -> None:
    """Emit FlexTaxD-compatible taxonomy relations from the tree."""
    from ..stages.tree2tax import Tree2taxParams

    def build() -> Tree2taxParams:
        return Tree2taxParams(
            node_basename=node_basename,
            root_name=root_name,
            remove_outgroup=remove_outgroup,
            include_dereplicated=include_dereplicated,
        )

    _run("tree2tax", workdir, build)


@app.command()
def glance(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    threads: int = typer.Option(24, "-t", "--threads"),
    plot_max: float = typer.Option(1.0, "--plot-max"),
    plot_min: float = typer.Option(0.0, "--plot-min"),
    keep_files: bool = typer.Option(False, "--keep-files"),
) -> None:
    """Quick all-vs-all ANI overview (dRep compare dendrogram + plots)."""
    from ..stages.glance import GlanceParams

    def build() -> GlanceParams:
        return GlanceParams(
            threads=threads, plot_max=plot_max, plot_min=plot_min, keep_files=keep_files
        )

    _run("glance", workdir, build)


@app.command(name="derep-unpack")
def derep_unpack(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    no_representant: bool = typer.Option(False, "--no-representant"),
) -> None:
    """Explode clusters into one directory per representative."""
    from ..stages.derep_unpack import DerepUnpackParams

    def build() -> DerepUnpackParams:
        return DerepUnpackParams(no_representant=no_representant)

    _run("derep_unpack", workdir, build)


@app.command(name="derep-stock")
def derep_stock(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help="Working directory."),
    action: str = typer.Option(..., "--action", help="list, pack, unpack or delete."),
    name: str | None = typer.Option(None, "--name", help="Run name for pack/unpack/delete."),
) -> None:
    """Store, load, list or delete named dereplication runs."""
    from ..stages.derep_stock import DerepStockParams

    def build() -> DerepStockParams:
        return DerepStockParams(action=action, name=name)

    _run("derep_stock", workdir, build)


def _parse_key_values(items: list[str], label: str) -> dict[str, str]:
    """Parse repeated ``key=value`` options into a dict (used for tool extras)."""
    out: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise UserInputError(f"{label} must be key=value, got '{item}'.")
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise UserInputError(f"{label} has an empty key in '{item}'.")
        out[key] = value.strip()
    return out


def _read_path_fofn(path: Path) -> list[Path]:
    """Read a file-of-filenames (one path per line; blank lines ignored)."""
    if not path.exists():
        raise UserInputError(f"File not found: {path}")
    return [Path(line.strip()) for line in path.read_text().splitlines() if line.strip()]


@app.command(name="genome-fetch")
def genome_fetch_cmd(
    selection: Path = typer.Option(
        ..., "--selection", help="selection.tsv from the metadata stage."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir for downloaded genomes."),
    keep_files: bool = typer.Option(False, "--keep-files", help="Keep download intermediates."),
) -> None:
    """Download genomes listed in a selection.tsv (stateless data-channel step)."""
    from ..stages.genome_steps import GenomeFetchParams, genome_fetch

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    try:
        genome_fetch(
            GenomeFetchParams(selection_tsv=selection, out_dir=out_dir, keep_files=keep_files),
            logger,
        )
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


@app.command(name="dereplicate-chunk")
def dereplicate_chunk_cmd(
    genomes_fofn: Path = typer.Option(
        ..., "--genomes-fofn", help="File of genome FASTA paths, one per line."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output directory for the chunk result."),
    tool: str = typer.Option("skder", "--tool", help="skder/drep/galah/sourmash."),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "-af", "--aligned-fraction"),
    threads: int = typer.Option(16, "-t", "--threads"),
    virus: bool = typer.Option(False, "--virus", help="Pass virus-tuned parameters to the tool."),
) -> None:
    """Dereplicate one chunk of genomes (scatter step; writes a chunk result dir)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.derep_steps import ChunkParams, dereplicate_chunk

    # A concrete tool (not 'auto') so every scattered chunk and the merge agree.
    _require_choice(tool, set(_derep_registry.names()), "--tool")
    _require_unit_interval(primary_ani, "--primary-ani")
    _require_unit_interval(secondary_ani, "--secondary-ani")
    _require_unit_interval(aligned_fraction, "--aligned-fraction")
    genomes = _read_path_fofn(genomes_fofn)

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    try:
        dereplicate_chunk(
            ChunkParams(
                tool=tool, genomes=genomes, out_dir=out_dir,
                primary_ani=primary_ani, secondary_ani=secondary_ani,
                aligned_fraction=aligned_fraction, threads=threads,
                extra={"virus": virus} if virus else {},
            ),
            logger,
        )
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


@app.command(name="dereplicate-merge")
def dereplicate_merge_cmd(
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir for the merged result."),
    chunk_dir: list[Path] = typer.Option(
        [], "--chunk-dir", help="A chunk result directory (repeatable)."
    ),
    chunk_fofn: Path | None = typer.Option(
        None, "--chunk-fofn", help="File listing chunk result directories, one per line."
    ),
    tool: str = typer.Option("skder", "--tool", help="skder/drep/galah/sourmash."),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani"),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani"),
    aligned_fraction: float = typer.Option(0.50, "-af", "--aligned-fraction"),
    threads: int = typer.Option(16, "-t", "--threads"),
    virus: bool = typer.Option(False, "--virus", help="Pass virus-tuned parameters to the tool."),
) -> None:
    """Dereplicate the union of chunk representatives (gather step)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.derep_steps import MergeParams, dereplicate_merge

    _require_choice(tool, set(_derep_registry.names()), "--tool")
    _require_unit_interval(primary_ani, "--primary-ani")
    _require_unit_interval(secondary_ani, "--secondary-ani")
    _require_unit_interval(aligned_fraction, "--aligned-fraction")
    chunk_dirs = list(chunk_dir)
    if chunk_fofn is not None:
        chunk_dirs += _read_path_fofn(chunk_fofn)
    if not chunk_dirs:
        raise UserInputError("Provide at least one --chunk-dir or a --chunk-fofn.")

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    try:
        dereplicate_merge(
            MergeParams(
                tool=tool, chunk_dirs=chunk_dirs, out_dir=out_dir,
                primary_ani=primary_ani, secondary_ani=secondary_ani,
                aligned_fraction=aligned_fraction, threads=threads,
                extra={"virus": virus} if virus else {},
            ),
            logger,
        )
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=1) from exc


@app.command(name="list-tools")
def list_tools() -> None:
    """List the available pluggable tools in each family."""
    from ..aligners.base import registry as aligners
    from ..dereplicators.base import registry as dereplicators
    from ..snptypers.base import registry as snptypers
    from ..treebuilders.base import registry as treebuilders

    for label, reg in (
        ("dereplicators", dereplicators),
        ("aligners", aligners),
        ("snptypers", snptypers),
        ("treebuilders", treebuilders),
    ):
        typer.echo(f"{label}: {', '.join(reg.names()) or '(none)'}")


if __name__ == "__main__":  # pragma: no cover
    sys.exit(app())
