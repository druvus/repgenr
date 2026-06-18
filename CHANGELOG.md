# Changelog

All notable changes to RepGenR are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **Nextflow nf-core rewrite (Phase 4, in progress)**: parameter schema
  (`nextflow_schema.json`) with nf-schema validation, execution reports, and
  nf-core template files (4a); discrete `repgenr dereplicate-chunk` /
  `dereplicate-merge` CLI steps (4b); a data-channel scatter-gather dereplication
  subworkflow with `stub:` blocks (4c); and stub-based nf-test plus a CI job
  running them on Nextflow 26.04 (4d). The per-stage conversion of the remaining
  stages to data channels is ongoing.
- **Sparse sourmash dereplication back-end**: when the optional
  `sourmash_plugin_branchwater` plugin is installed, the sourmash dereplicator
  uses `manysketch` + `pairwise` to compute only above-threshold edges instead of
  the dense N x N `compare` matrix, keeping memory roughly linear in the number of
  close pairs (relevant at 10k+ genomes). Selected automatically; falls back to
  the dense `compare` path when the plugin is absent. Both paths yield the same
  cluster partition and representative count for a given threshold. Install via
  the `sparse` extra or the `sourmash_plugin_branchwater` conda package.

## [2.0.0] - 2026-06-18

First stable release of the v2 rewrite: a modular `repgenr` Python package
(Typer CLI, entry-point plugin registries for dereplicators / aligners / SNP
typers / tree builders), a SQLite genome manifest, `repgenr.yaml` provenance, a
no-shell subprocess layer, and a Nextflow orchestration layer.

### Added
- **Container execution backend**: run any external tool in a pinned container
  (`--container docker|singularity`, `--wave`), with BioContainers and Seqera
  Wave image resolution; per-call `HOME` and `extra_mounts`; container-cache
  control. All tool families verified in containers.
- **Dereplication scaling**: `--process-size` two-stage chunking for any tool;
  stage-1 ANI thresholds (`--pre-primary-ani`/`--pre-secondary-ani`);
  `--num-processes` parallel chunk workers.
- **Resume/idempotency**: stages that already completed with the same parameters
  are skipped; `--force` to re-run.
- **Validation & logging**: enum/range validation of CLI options; `--verbose`/
  `--quiet`/`REPGENR_LOG_LEVEL`.
- **Tool version floors**: `min_version` preflight enforcement for the tools with
  reliable version strings; lower-bounded `environment.yml`.
- Nextflow: CPU-matched threads, dynamic resources + retry, first-class chunking
  params, alignment-free default tree builder.

### Changed
- Dereplication post-processing is O(n) (two-stage compose) and the sourmash
  clustering uses a numpy matrix; the manifest uses batched transactions + WAL;
  genome staging uses hardlinks (copy fallback). These remove the per-genome
  Python/I-O bottlenecks for 1000-10000 genome sets.

### Fixed
- `maf_to_fasta` produced an empty MSA on versioned accessions; cactus picked a
  per-chromosome HAL and leaked the `_MINIGRAPH_` pseudo-genome; the manifest
  swallowed write errors and had no schema versioning/concurrency timeout.

### Notes
- Reproducibility: tool versions are recorded in `repgenr.yaml`; generate a
  pinned per-platform conda lock for exact reproducibility (see `environment.yml`).
