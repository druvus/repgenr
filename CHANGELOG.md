# Changelog

All notable changes to RepGenR are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- `metadata --limit N` no longer keeps the first N genomes in GTDB file (or
  API) order. It round-robins over species, taking the best CheckM-scored
  genome of every species first, then each species' next best, until N;
  within a species unscored genomes rank last, then the GTDB representative
  flag, then accession. The same `--limit` therefore returns a different,
  better set than before. On the API path the per-genome quality cards are
  fetched for every candidate before the cut.

### Fixed
- Six small defects noted in the 2026-09-01 audit's self-review: `phylo`
  publishes `tree/tree.nwk` through the atomic copy used by every other
  deliverable; the Nextflow retry window no longer includes exit 130 and 131
  (a cancelled task is not resubmitted); an unknown `--container` value is a
  user-input error with the valid choices in the message; IQ-TREE refuses
  `--bootstrap` values below its floor of 1000 before running; the Wave image
  cache is keyed by platform as well as conda spec; and a failure inside the
  tool-output read loop kills the child process before the error propagates.

## [3.0.0] - 2026-09-04

### Added
- `ska2` SNP typer: reference-free split k-mer calling (`repgenr snptype
  --tool ska2`, `phylo --msa-source snptype --snptyper ska2`). No genome is
  privileged as the reference, so reference-private errors do not bias the
  SNP distances. It emits a variable-site alignment only, so `--mask` is
  refused for it. Tuning via `--tool-arg ksize=` and `min_freq=`.
- `repgenr doctor`: read-only workdir health check that verifies outputs
  against the records in `repgenr.yaml` -- interrupted stages, missing or
  corrupt genomes, manifest drift, representative/cluster mismatches,
  truncated deliverables, unresolvable outgroups, leftover temp files, and
  stages whose inputs changed since completion. Exits 1 on failures.
- `repgenr run` accepts `--msa-source` and `--snptyper`, so the SNP-typing
  phylogeny path (previously manual-only via `repgenr phylo`) is reachable from
  the one-shot orchestrator with identical resume fingerprints.
- **Nextflow nf-core rewrite (Phase 4)**: the pipeline is now a typed
  data-channel workflow with no shared working directory. Parameter schema
  (`nextflow_schema.json`) with nf-schema validation, execution reports and
  nf-core template files; stateless `repgenr dereplicate-chunk` /
  `dereplicate-merge` / `genome-fetch` CLI steps and a portable `selection.tsv`
  hand-off; a scatter-gather dereplication subworkflow; data-channel
  `BACTERIAL_DATAFLOW` (metadata -> genome -> derep -> phylo -> tree2tax) and
  `VIRAL_DATAFLOW` pipelines selected by `--mode`; results published under
  `--outdir`. Stub-based nf-test and a CI job run them on Nextflow 26.04. The
  legacy shared-workdir orchestrator and done-signal modules are removed.
- **Sparse sourmash dereplication back-end**: when the optional
  `sourmash_plugin_branchwater` plugin is installed, the sourmash dereplicator
  uses `manysketch` + `pairwise` to compute only above-threshold edges instead of
  the dense N x N `compare` matrix, keeping memory roughly linear in the number of
  close pairs (relevant at 10k+ genomes). Selected automatically; falls back to
  the dense `compare` path when the plugin is absent. Both paths yield the same
  cluster partition and representative count for a given threshold. Install via
  the `sparse` extra or the `sourmash_plugin_branchwater` conda package.
- **Quality-aware representative selection.** `--keeper quality` (default)
  re-picks each cluster's representative by CheckM completeness minus five
  times contamination after any dereplicator runs, at the in-process stage,
  the Nextflow chunk level and taxonomy reduction; `--keeper tool` keeps the
  adapter's own pick. GTDB CheckM values are carried into the manifest
  (schema v2, two quality columns) and `selection.tsv` (two optional
  columns). Provenance records `keeper`, `keeper_effective` and the number
  of representatives changed.
- `snp/full_alignment.fasta`: SNP typers emit the whole-genome alignment
  (snippy `core.full.aln`, ParSNP `harvesttools -M`, the simple typer's
  consensuses) as a deliverable and as the input for recombination maskers.
- Nextflow publishes `pipeline_info/software_versions.yml`, collected from
  every process.

### Changed
- **Crash/restart hardening (stage audit)**: a stage that crashes while
  re-running is no longer silently skipped on restart (the resume record is
  dirtied before the stage body runs; `repgenr status` shows `[interrupted]`);
  deliverable writes are atomic (a failing tree builder can no longer truncate
  a previous `tree/tree.nwk`); downloads are validated before landing in
  `genomes/` and stages refuse incomplete input sets (`--allow-incomplete` to
  override); the manifest is reconciled on re-selection instead of accumulating
  stale rows; stale outgroups are pruned and resolved by exact accession.
- **Pluggability**: tool lists in `--help` are generated from the registries;
  `--tool-arg key=value` provides tool tuning on dereplicate/snptype and the
  data-channel steps; `auto` selection is container-aware and prefers the
  tightest-fitting tool (its picks change); the Nextflow `derep_tool` enum is
  dropped so pip-installed dereplicators work from the pipeline; recombination
  masking is a plugin family (`repgenr.maskers`, `phylo` gains a real
  `--mask`); `repgenr glance` gains `--tool` via a `Dereplicator.compare`
  hook and the viral outgroup builder is selectable via
  `--outgroup-treebuilder` (a `TreeBuilder.distance_matrix` hook); public
  `Registry.register/unregister`; the dead `threads_param` and
  `Aligner.output_kind` capability fields are removed.
- **Breaking: resume fingerprints are input-aware.** A completed stage is now
  skipped only when its parameters, the digests of its inputs (upstream stage
  outputs, digested from file metadata for genome directories and content for
  small contract files), and the container identity (backend/platform/wave) are
  all unchanged. Re-running an upstream stage automatically re-runs downstream
  stages; the previous timestamp-based "may be stale" warning is removed.
  Workdirs created by older versions re-run each stage once.
- **Breaking: `repgenr run` and the manual commands now build identical stage
  parameters** through shared builders, so the two entry points share resume
  fingerprints. `tree2tax --include-dereplicated` now defaults to on for both
  (previously only `repgenr run` enabled it); pass
  `--no-include-dereplicated` for the old manual behavior. The bacterial
  `repgenr run` now requires `-l/--level` instead of silently passing an empty
  level.
- **Breaking for third-party adapters:** `Masker.mask(full_alignment,
  out_dir, params, logger)` receives the whole-genome alignment and a
  `MaskParams`; `SnpParams.mask` is removed; `xmfa_to_fasta` loses its unused
  `flank` argument.
- Stage-level flags (`--virus`, `--mask`) enter a tool's extras only when that
  tool declares them, so resume fingerprints change only when behaviour
  does; every stage warns by name about extras no selected tool reads.
- Repository hygiene: `ruff format` is enforced in CI (the mechanical
  reformat is listed in `.git-blame-ignore-revs`), CI caches pip and the
  micromamba environment and installs nf-test from a pinned release tarball,
  the wiki binaries are no longer tracked (figures live in `docs/images/`),
  and a test pins `pyproject.toml`, the Nextflow manifest and this changelog
  to one version. The Nextflow floor is `>=23.10.0`, the minimum nf-schema
  2.3.0 supports.

### Fixed
- `tree2tax` roots on the outgroup edge instead of at the outgroup's parent
  node. For the unrooted trees that mashtree, fasttree and the sourmash tree
  builder emit, the old rooting left the root with three children (the
  outgroup, its nearest ingroup neighbour and everything else), so the
  taxonomy had no node for the whole ingroup and one taxon sat beside the
  outgroup. The phylo docstring and architecture note now say where rooting
  happens.
- `metadata --source api` reads CheckM quality from each genome's GTDB card
  (`metadata_gene.checkm2_*`, `checkm_*` fallback); the genomes-detail rows
  carry none, so API selections had empty quality columns and the quality
  keeper silently kept the adapter's picks. `dereplicate` now warns when the
  manifest has no quality and records `keeper_effective` in `repgenr.yaml`.
- `repgenr run` preflights every external tool (dereplicator, tree builder,
  and the aligner or SNP typer when the builder needs an MSA) before the
  first stage, instead of discovering a missing tree builder after download
  and dereplication.
- Alignment-free `phylo` runs no longer record an aligner in provenance, so
  `--aligner` cannot invalidate a mashtree or sourmash tree; mashtree builds
  and distance matrices share one argument set and come from one call;
  streamed tool output goes to DEBUG (the file log keeps it, the console
  shows the pipeline's own messages), with carriage-return redraws logged
  once in their final state.
- Dereplication refuses incomplete results (a genome without a status, an
  unmarked representative, a contained genome with no or several
  representatives) instead of silently dropping genomes; genomes that fail a
  tool's QC filter are carried through chunked composition, taxonomy
  reduction and the Nextflow merge.
- `doctor` is strictly read-only: it opens the manifest in read-only mode and
  never creates, migrates or journals it.
- Nextflow help text names NCBI Virus; `tree2tax` publishes only its TSVs to
  the outdir root.

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
