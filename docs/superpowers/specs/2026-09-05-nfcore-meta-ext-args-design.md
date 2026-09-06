# Nextflow layer: nf-core meta maps and ext.args

Date: 2026-09-05. Status: approved design, awaiting implementation plan.

## Goal

Bring the Nextflow layer to the shape the current nf-core template and
`nf-core pipelines lint` expect: every process exchanges `tuple val(meta),
path(...)`, tool flags reach processes through `task.ext.args` from
`conf/modules.config` instead of `params.*` reads inside module scripts,
publishing lives in configuration rather than in modules, and the pipeline
declares Nextflow 26.04 as its floor. Command lines and the parameter schema
are unchanged for users.

## Decisions taken during design

| Question | Decision |
|---|---|
| What the meta map represents | One run-level meta built from the parameters: `[id: <target slug>, mode: params.mode]`. One target per run, as today. A multi-target samplesheet is a later, additive change. |
| Relation of `ext.args` to the `*_args` and `derep_*` params | The params stay as the user interface and the schema is unchanged; `conf/modules.config` maps each to its process with a closure. Modules read only `task.ext.*`. |
| Version reporting | Unchanged: each process still writes a `versions.yml` fragment and `PUBLISH_VERSIONS` collects them. The `eval()`/topic migration is a follow-up that the new floor makes free. |
| Nextflow floor | `!>=26.04.0` (the `!` makes an older Nextflow a hard failure). nf-schema moves to 2.6.1, verified to load and validate under 26.04.6. |
| Configuration layout | Full nf-core shape: `conf/base.config` and `conf/modules.config` included from `nextflow.config`. |

## Section 1: configuration and meta

`nextflow.config` keeps `params`, `plugins`, `manifest`, the execution
reports and the profiles, and gains two includes.

`conf/base.config`:

- `errorStrategy` with the retry window `(132..145) + [104, 247]` and
  `maxRetries = 2`, moved from `nextflow.config` unchanged.
- The three resource labels `process_low`, `process_medium`, `process_high`
  with the current CPU, memory and time values (memory and time scale with
  `task.attempt` as today).
- `resourceLimits` set per profile: the `test` profile caps at
  `[cpus: 2, memory: '4.GB', time: '1.h']` and drops its per-label
  overrides; `standard`, `slurm` and `cloud` set no cap.

`conf/modules.config` has one `withName` block per process. Each block
carries the process's `ext.args` closure and its `publishDir`:

| Process | `ext.args` | `publishDir` |
|---|---|---|
| METADATA | `{ params.metadata_args }` | `${params.outdir}/metadata`, mode copy |
| GENOME | none | none (raw genomes are intermediates) |
| VACQUIRE | `{ params.vmetadata_args }`; `ext.args2 = { params.vgenome_args }` | none |
| DEREP_CHUNK | composed closure, see below | none |
| DEREP_MERGE | same composed closure | `${params.outdir}/dereplicate`, mode copy |
| PHYLO | `{ params.phylo_args }` | `${params.outdir}/phylo`, mode copy |
| TREE2TAX | `{ params.tree2tax_args }` | `${params.outdir}`, mode copy, pattern `*.tsv` |

The dereplication closure composes five `derep_*` params into one flag
string: `--tool`, `--primary-ani`, `--secondary-ani`, `--aligned-fraction`,
`--keeper`, plus `--virus` when `params.mode == 'viral'`.
`derep_process_size` stays a param read by the scatter subworkflow, since it
shapes channels rather than a command line.

A process-wide `ext.repgenr_opts = { params.repgenr_opts }` in the
`process {}` block replaces the direct read of `params.repgenr_opts` in every
module. The container profiles keep setting `params.repgenr_opts`.

The meta is built once in `main.nf` by a function `run_meta(params)`:

- `id`: a slug of the target. Bacterial: the value after `-ts`, else `-tg`,
  else `-tf` in `params.metadata_args` (most specific rank first); viral: the value after `-tg` in
  `params.vgenome_args`, else after `-t` in `params.vmetadata_args`. Lowercase,
  non-alphanumerics replaced by `_`. Fallback: `params.mode`.
- `mode`: `params.mode`.

Chunk metas are `[id: "${meta.id}.chunk_${i}", mode: meta.mode]` and the
merge meta is `[id: "${meta.id}.merged", mode: meta.mode]`, derived in
`DEREPLICATE_SCATTER` from the incoming meta.

The parameter schema (`nextflow_schema.json`) is unchanged. Its description
strings for the `*_args` params gain the sentence "Reaches the process as
`task.ext.args`; a site config can override it with `withName`."

## Section 2: modules and subworkflows

Every process:

- takes `val(meta)` (METADATA, VACQUIRE, which have no other inputs) or
  `tuple val(meta), path(...)` for each data input, and emits
  `tuple val(meta), path(...)` for each data output. `versions.yml` stays a
  bare `path` output, matching nf-core.
- has `tag "${meta.id}"`, exactly one resource label, and
  `when: task.ext.when == null || task.ext.when`.
- opens `script:` with `def args = task.ext.args ?: ''` (and `args2` where
  used) and `def opts = task.ext.repgenr_opts ?: ''`; the command uses
  `${opts}` and `${args}`.
- contains no `params.*` read and no `publishDir`.
- ends its script with one call `repgenr_versions_fragment <workdir-or-file>`
  to a shell function in `nextflow/bin/repgenr_versions_fragment`, which
  writes `versions.yml` for `${task.process}` from `repgenr --version` and the
  stage's `tool_versions.yml`. `nextflow/bin/` is on the task PATH by Nextflow
  convention. The stub blocks keep `touch versions.yml`.
- keeps its stub block producing every declared output with the new tuple
  shape; stub content is otherwise unchanged.

Subworkflows:

- `ACQUIRE(ch_meta)` emits `genomes` as `tuple(meta, [fasta...])`,
  `outgroup` as `tuple(meta, [fasta...])` (empty list when none),
  `selection` and `outgroup_accession` as `tuple(meta, path)`.
- `VACQUIRE` emits the same four shapes.
- `DEREPLICATE_SCATTER(ch_genomes, ch_selection)` takes `tuple(meta,
  [fasta...])` and `tuple(meta, selection_or_[])`, derives chunk metas, and
  emits `reps` as `tuple(meta, dir)` with the run meta (not the merge meta),
  so downstream joins on the run meta work.
- `PHYLO(tuple(meta, reps_dir, [outgroup...], outgroup_accession))` and
  `TREE2TAX(tuple(meta, tree, reps_dir, [outgroup...], outgroup_accession))`:
  one joined tuple per process, built in the dataflow subworkflows with
  `.join(by: 0)` on the meta, which is the nf-core shape for a process with
  several per-run inputs. The dereplication processes keep `selection.tsv` as a
  bare auxiliary input (reference-file style), and DEREP_MERGE's output
  directory is `task.ext.prefix` (set to `merged` in `conf/modules.config`) so
  published paths do not change.
- Both dataflow subworkflows take `ch_meta` from `main.nf`.
- `PUBLISH_VERSIONS` is unchanged.

The two harnesses under `nextflow/tests/` (`local_dataflow.nf`,
`acquire_scatter.nf`, `bacterial_dataflow.nf`, `dereplicate_scatter.nf`,
`viral_dataflow.nf`) build a meta the same way and pass tuples.

Strict-syntax cleanup, required because 26.04 is now the floor: lowercase
`channel.` factories, every closure parameter named (unused ones prefixed with
`_`), `def` on workflow-scope channel variables. `nextflow lint nextflow/`
must report no errors and no warnings.

## Section 3: tests, documentation, compatibility

nf-test:

- All 15 test files move to tuple inputs. The six `.snap` files are
  regenerated and each is inspected before commit; the content assertions
  that sit beside `snapshot().match()` stay.
- New test `modules_config.nf.test`: runs PHYLO in stub mode with a test
  config that sets `withName: PHYLO { ext.args = '--treebuilder faketree' }`
  and asserts the string reaches the command line (`process.out` plus the
  task's `.command.sh`). This is the behaviour the change exists to enable.
- New test in `bacterial_dataflow.nf.test`: the meta id built from
  `metadata_args` (`francisella` under the test profile) is the `meta.id` on
  `TREE2TAX.out.tree2tax`.
- The e2e test keeps running the real tools through `local_dataflow.nf`.

Python: no code change. `tests/unit/test_version_consistency.py` already reads
the manifest version and keeps passing. `tests/unit/test_nextflow_retry_window.py`
must find the `errorStrategy` line in `conf/base.config`; it is updated to read
that file.

Documentation:

- `docs/usage.md`: new section "Configuring processes" explaining where
  `ext.args` lives, how a site overrides one process with a `-c` config and
  `withName`, the `ext.repgenr_opts` hook, and the `!>=26.04.0` floor.
- `docs/architecture.md`: one paragraph on the run-level meta and the tuple
  contract between processes.
- `README.md` Nextflow section: the floor.
- `CHANGELOG.md` under Unreleased, Changed: floor raised to Nextflow 26.04
  and nf-schema 2.6.1; configuration split into `conf/base.config` and
  `conf/modules.config` with per-process `ext.args` and `publishDir`; every
  channel carries a meta map; command lines and parameters unchanged;
  breaking only for code that included a module or subworkflow directly. A
  Notes line records that the versions-topic migration is now a free
  follow-up.

Compatibility: no change for command-line users of `nextflow run
nextflow/main.nf`. Anyone who included a RepGenR module or subworkflow from
their own workflow must adapt to the tuple shapes. CI already runs Nextflow
26.04.3 and nf-test 0.9.5; no workflow change is needed.

## Out of scope

- Multi-target samplesheet input.
- Version reporting through `eval()` outputs and the `versions` topic.
- The workflow `output {}` block; `publishDir` in `conf/modules.config` is the
  mechanism used throughout, deliberately not half-migrated.
- `meta.yml` files for local modules; `nf-core pipelines lint` does not
  require them for `modules/local`.

## Verification before the pull request

- `nextflow lint nextflow/`: no errors, no warnings.
- `nf-test test --tag stub`: all stub tests pass, snapshots reviewed.
- `nextflow run nextflow/main.nf -profile test -stub --outdir /tmp/x`: completes
  and publishes `tree2tax.tsv`, `genomes_map.tsv`, `phylo/tree/tree.nwk`,
  `pipeline_info/software_versions.yml`.
- A `-c` override test by hand: `withName: PHYLO { ext.args = '--treebuilder
  sourmash' }` changes the PHYLO command line in `.command.sh`.
- Python gate unchanged: pytest, ruff check and format, mypy.
