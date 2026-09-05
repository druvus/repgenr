# nf-core meta maps and ext.args Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring RepGenR's Nextflow layer to the current nf-core shape: a meta map on every channel, tool flags via `task.ext.args` from `conf/modules.config`, publishing in configuration, Nextflow 26.04 as the floor, with unchanged command lines and parameters.

**Architecture:** `nextflow.config` keeps params, plugins, manifest and profiles and includes `conf/base.config` (resources, retries) and `conf/modules.config` (per-process `ext.args` closures over the existing params, and `publishDir`). A `run_meta(params)` function builds one run-level meta in `main.nf`; every process takes and emits `tuple val(meta), path(...)` and reads only `task.ext.*`. The work proceeds in layers that each leave the stub suite green: floor and lint, config and ext, versions helper, meta at the acquisition front, meta through dereplication, meta through phylo and tree2tax, docs.

**Tech Stack:** Nextflow 26.04.6 (strict parser), nf-schema 2.6.1, nf-test 0.9.5 (conda env `repgenr_nf`), the `repgenr` Python CLI (conda env `repgenr_dev`), pytest for the one Python test touched.

**Spec:** `docs/superpowers/specs/2026-09-05-nfcore-meta-ext-args-design.md`

## Global Constraints

- Nextflow floor becomes `nextflowVersion = '!>=26.04.0'`; nf-schema becomes `id 'nf-schema@2.6.1'` (verified to load and validate under 26.04.6).
- Strict parser rules everywhere: lowercase `channel.` factories, every closure parameter named (unused ones prefixed `_`), `def` on workflow-scope variables, no slashy strings, no `for` loops, no spread.
- `nextflow lint nextflow/` must finish with no errors and no warnings after Task 1 and stay that way.
- No Unicode in any file under `nextflow/` (comments included).
- Modest scientific language in docs, comments and log messages.
- Every module reads only `task.ext.*`; no `params.*` and no `publishDir` inside a module after Task 2.
- Version reporting stays file-style `versions.yml` fragments collected by `PUBLISH_VERSIONS`.
- Command lines and `nextflow_schema.json` parameters are unchanged for users.
- Run Nextflow tooling as `conda run -n repgenr_nf --no-capture-output <cmd>` from the repo root; run Python gates with `/Users/andreassjodin/miniforge3/envs/repgenr_dev/bin/python -m pytest -q`, `ruff check src/ tests/`, `ruff format --check src/ tests/ benchmarks/`, `mypy src/repgenr`.
- Commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_01H41vhDaRfyPXrcfFc1MXNG`.
- Work on branch `feat/nfcore-modules` (already created; the spec is its first commit).

## File map

| File | Responsibility after this plan |
|---|---|
| `nextflow/nextflow.config` | params, plugins, includes, profiles, manifest, reports |
| `nextflow/conf/base.config` (new) | retry window, resource labels |
| `nextflow/conf/modules.config` (new) | per-process `ext.args`, `ext.prefix`, `publishDir`; process-wide `ext.repgenr_opts` |
| `nextflow/bin/repgenr_versions_fragment` (new) | shell helper writing a process's `versions.yml` |
| `nextflow/subworkflows/local/run_meta.nf` (new) | `run_meta(params)` function |
| `nextflow/main.nf` | help, validation, build meta, dispatch by mode |
| `nextflow/modules/local/dataflow/{metadata,genome,vacquire,phylo,tree2tax}.nf` | tuple-shaped processes reading `task.ext.*` |
| `nextflow/modules/local/{derep_chunk,derep_merge}.nf` | same |
| `nextflow/subworkflows/local/{acquire,dereplicate_scatter,bacterial_dataflow,viral_dataflow}.nf` | tuple plumbing, joins on meta |
| `nextflow/tests/*.nf`, `nextflow/tests/*.nf.test`, `*.snap` | harnesses and tests on the new shapes |
| `tests/unit/test_nextflow_retry_window.py` | reads `conf/base.config` |
| `docs/usage.md`, `docs/architecture.md`, `README.md`, `CHANGELOG.md`, `nextflow_schema.json` descriptions | documentation |

---

### Task 1: Raise the floor, bump nf-schema, clear every lint warning

**Files:**
- Modify: `nextflow/nextflow.config:3-6` (plugins) and `:146-148` (manifest)
- Modify: `nextflow/main.nf`, `nextflow/subworkflows/local/*.nf`, `nextflow/tests/*.nf` (lint warnings only; no behaviour change)

**Interfaces:**
- Produces: a layer that lints clean under 26.04 so later tasks can require "no warnings" as their exit criterion.

- [ ] **Step 1: Record the current lint baseline (the failing check)**

Run: `cd /Users/andreassjodin/Code/repgenr && conda run -n repgenr_nf --no-capture-output nextflow lint nextflow/ 2>&1 | tail -3`
Expected: `9 files had 21 warnings`, `19 files had no errors`. The task is done when this reports no warnings.

- [ ] **Step 2: Bump the plugin and the floor**

In `nextflow/nextflow.config` change:

```groovy
plugins {
    // Parameter validation and schema-driven help/summary (nf-core convention).
    id 'nf-schema@2.6.1'
}
```

and in the `manifest` block:

```groovy
    nextflowVersion = '!>=26.04.0'
    version         = '3.0.0'
```

- [ ] **Step 3: Fix the warnings**

Run the linter per file and apply these three mechanical fixes everywhere it points:

1. `Channel.empty()`, `Channel.value(...)`, `Channel.fromPath(...)` become `channel.empty()`, `channel.value(...)`, `channel.fromPath(...)`.
2. Closures using implicit `it` get a named parameter: `.filter { !it.name.startsWith('._') }` becomes `.filter { f -> !f.name.startsWith('._') }`; `.map { it.text }` becomes `.map { f -> f.text }`.
3. Closure parameters that are not used get an underscore prefix: `.map { meta, dir -> dir }` becomes `.map { _meta, dir -> dir }`.

Also add `def` to workflow-scope channel assignments in `main.nf`, the subworkflows and the harnesses (`def ch_versions = channel.empty()`), keeping `emit:` blocks as plain `name = value`.

- [ ] **Step 4: Verify lint and the stub suite**

Run: `conda run -n repgenr_nf --no-capture-output nextflow lint nextflow/ 2>&1 | tail -3`
Expected: `28 files had no errors` (or the current file count) and no warnings line.

Run: `conda run -n repgenr_nf --no-capture-output nf-test test --tag stub 2>&1 | grep -E "Executed|FAIL"`
Expected: `Executed 17 tests ... (0 failed)`.

- [ ] **Step 5: Commit**

```bash
git add nextflow
git commit -m "chore(nextflow): floor at 26.04, nf-schema 2.6.1, lint-clean strict syntax

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01H41vhDaRfyPXrcfFc1MXNG"
```

---

### Task 2: Split the config; modules read task.ext, not params

**Files:**
- Create: `nextflow/conf/base.config`, `nextflow/conf/modules.config`
- Create: `nextflow/tests/modules_config.nf.test`, `nextflow/tests/modules_config_override.config`
- Modify: `nextflow/nextflow.config` (move the `process {}` block out; add includes; test profile uses `resourceLimits`)
- Modify: all seven modules under `nextflow/modules/local/` (read `task.ext.*`, drop `publishDir`, echo args in stubs)
- Modify: `tests/unit/test_nextflow_retry_window.py:10-11`

**Interfaces:**
- Produces: `task.ext.args` (string) on every process; `task.ext.args2` on VACQUIRE; `task.ext.repgenr_opts` (string) on every process; `task.ext.prefix` honoured by DEREP_MERGE (default `"${meta.id}"`).
- Later tasks rely on the exact `withName` blocks below and on the stub line `echo "ext.args: ${args}"`.

- [ ] **Step 1: Write the failing nf-test for an ext.args override**

Create `nextflow/tests/modules_config_override.config`:

```groovy
// Test-only: prove a site config can retune one process through ext.args.
process {
    withName: 'PHYLO' {
        ext.args = '--treebuilder faketree-from-config'
    }
}
```

Create `nextflow/tests/modules_config.nf.test`:

```groovy
nextflow_process {

    name "ext.args reaches the process command line"
    script "../modules/local/dataflow/phylo.nf"
    tag "stub"
    process "PHYLO"
    options "-stub"
    config "./modules_config_override.config"

    test("a withName ext.args override is what the stub echoes") {
        when {
            process {
                """
                input[0] = file("${projectDir}/nextflow/tests/data/derep_out")
                input[1] = file("${projectDir}/nextflow/tests/data/outgroup/Fam_Out_grp_GCF_000009.1.fasta")
                input[2] = file("${projectDir}/nextflow/tests/data/outgroup_accession.txt")
                """
            }
        }
        then {
            assert process.success
            assert process.stdout.any { line -> line.contains('ext.args: --treebuilder faketree-from-config') }
        }
    }
}
```

(Task 6 changes `input[0..2]` to one tuple; the test is updated there.)

- [ ] **Step 2: Run it to see it fail**

Run: `conda run -n repgenr_nf --no-capture-output nf-test test nextflow/tests/modules_config.nf.test 2>&1 | grep -E "PASSED|FAILED"`
Expected: FAILED (the stub prints nothing about ext.args yet).

- [ ] **Step 3: Update the Python retry-window test to read the new file**

In `tests/unit/test_nextflow_retry_window.py` change the config path:

```python
def _retry_exit_codes() -> set[int]:
    config = (ROOT / "nextflow" / "conf" / "base.config").read_text(encoding="utf-8")
```

Run: `/Users/andreassjodin/miniforge3/envs/repgenr_dev/bin/python -m pytest -q tests/unit/test_nextflow_retry_window.py`
Expected: FAIL with `FileNotFoundError` until Step 4.

- [ ] **Step 4: Create conf/base.config**

```groovy
// Resources and retries for every process (nf-core base.config convention).
//
// Memory and time scale with the retry attempt so a large run that hits an
// out-of-memory or time-limit kill is resubmitted with more headroom. The
// modules export REPGENR_PROPAGATE_TOOL_EXIT=1 so a tool failure inside the
// repgenr CLI surfaces the tool's own exit code (signal kills as 128+signum,
// e.g. OOM kill = 137) and can match the retry window. 130 and 131 (SIGINT,
// SIGQUIT) are a person cancelling the task and are not retried.

process {
    errorStrategy = { task.exitStatus in ((132..145) + [104, 247]) ? 'retry' : 'terminate' }
    maxRetries = 2

    withLabel: process_low {
        cpus = 2
        memory = { 4.GB * task.attempt }
        time = { 2.h * task.attempt }
    }
    withLabel: process_medium {
        cpus = 8
        memory = { 16.GB * task.attempt }
        time = { 8.h * task.attempt }
    }
    withLabel: process_high {
        cpus = 32
        memory = { 128.GB * task.attempt }
        time = { 48.h * task.attempt }
    }
}
```

- [ ] **Step 5: Create conf/modules.config**

```groovy
// Per-process tool arguments and publishing (nf-core modules.config convention).
//
// Modules read only task.ext.*; the user-facing parameters are mapped here, so a
// site can retune one process from its own config (-c site.config) with a
// withName block, without editing the modules or the parameter schema.

process {
    // Top-level repgenr options (container backend) for every repgenr call.
    ext.repgenr_opts = { params.repgenr_opts }

    withName: 'METADATA' {
        ext.args = { params.metadata_args }
        publishDir = [
            path: { "${params.outdir}/metadata" },
            mode: 'copy',
            saveAs: { filename -> filename == 'versions.yml' ? null : filename },
        ]
    }

    withName: 'VACQUIRE' {
        ext.args  = { params.vmetadata_args }
        ext.args2 = { params.vgenome_args }
    }

    withName: 'DEREP_CHUNK|DEREP_MERGE' {
        ext.args = {
            [
                "--tool ${params.derep_tool}",
                "--primary-ani ${params.derep_primary_ani}",
                "--secondary-ani ${params.derep_secondary_ani}",
                "--aligned-fraction ${params.derep_aligned_fraction}",
                "--keeper ${params.derep_keeper}",
                params.mode == 'viral' ? '--virus' : '',
            ].join(' ').trim()
        }
    }

    withName: 'DEREP_MERGE' {
        // Stable published directory name regardless of the run meta id.
        ext.prefix = 'merged'
        publishDir = [
            path: { "${params.outdir}/dereplicate" },
            mode: 'copy',
            saveAs: { filename -> filename == 'versions.yml' ? null : filename },
        ]
    }

    withName: 'PHYLO' {
        ext.args = { params.phylo_args }
        publishDir = [
            path: { "${params.outdir}/phylo" },
            mode: 'copy',
            saveAs: { filename -> filename == 'versions.yml' ? null : filename },
        ]
    }

    withName: 'TREE2TAX' {
        ext.args = { params.tree2tax_args }
        publishDir = [
            path: { "${params.outdir}" },
            mode: 'copy',
            pattern: '*.tsv',
        ]
    }
}
```

- [ ] **Step 6: Trim nextflow.config**

Delete the whole `process { ... }` block (lines 46 to 74 today) and put, directly after the `params { }` block:

```groovy
includeConfig 'conf/base.config'
includeConfig 'conf/modules.config'
```

Replace the `test` profile body with:

```groovy
    test {
        // Minimal resources for a quick local smoke run on a small target.
        process.executor = 'local'
        process.resourceLimits = [ cpus: 2, memory: '4.GB', time: '1.h' ]
        params.metadata_args = '-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella --limit 5'
        params.derep_tool = 'sourmash'
        params.phylo_args = '--treebuilder mashtree'
    }
```

Update the comment above `derep_keeper` in `params` to end with "(mapped to the dereplication processes' ext.args in conf/modules.config)".

- [ ] **Step 7: Make every module read task.ext**

Apply to all seven modules. Shown in full for `phylo.nf`; the others follow the same pattern with their own command.

`nextflow/modules/local/dataflow/phylo.nf`:

```groovy
// Build the phylogeny (data-channel form).
//
// Calls the stateless `repgenr phylo-build` step directly on the staged channel
// files: the merged representatives directory provides the genome set, the
// outgroup files are staged into an outgroup/ directory, and tree/tree.nwk is
// emitted as a channel output. Tool flags arrive as task.ext.args from
// conf/modules.config; publishing is configured there too.

process PHYLO {
    tag "phylo"
    label 'process_high'

    input:
    path reps_dir
    path outgroup, stageAs: 'outgroup/*'
    path outgroup_accession

    output:
    path "tree/tree.nwk", emit: tree
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} phylo-build \\
        --genomes-dir ${reps_dir}/representatives \\
        --outgroup-dir outgroup \\
        --outgroup-accession ${outgroup_accession} \\
        -o . -t ${task.cpus} ${args} \\
        --versions-out tool_versions.yml

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    repgenr: \$(repgenr --version | sed 's/repgenr //')
END_VERSIONS
    cat tool_versions.yml >> versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    mkdir -p tree
    names=\$(ls ${reps_dir}/representatives | sed 's/\\.[^.]*\$//' | paste -sd, -)
    echo "(\${names});" > tree/tree.nwk
    touch versions.yml
    """
}
```

The same edits for the other modules:

- `metadata.nf`: remove `publishDir`; add `when:`; `def args = task.ext.args ?: ''` and `def opts = task.ext.repgenr_opts ?: ''`; command `repgenr ${opts} metadata -wd metadata_wd ${args}`; stub starts with `echo "ext.args: ${args}"` (declare `def args` in the stub too).
- `genome.nf`: add `when:`; `def opts = task.ext.repgenr_opts ?: ''`; command `repgenr ${opts} genome-fetch ...`; no args (GENOME has no flags); stub unchanged except a leading `echo "ext.args:"`.
- `vacquire.nf`: add `when:`; `def args = task.ext.args ?: ''`, `def args2 = task.ext.args2 ?: ''`, `def opts = task.ext.repgenr_opts ?: ''`; commands `repgenr ${opts} vmetadata -wd wd ${args}` and `repgenr ${opts} vgenome -wd wd ${args2}`; stub echoes both.
- `derep_chunk.nf`: delete the `def virus_flag` line and the six `--tool/--primary-ani/.../--keeper ${params...}` lines and `${virus_flag}`; use `def args = task.ext.args ?: ''` and `def opts = task.ext.repgenr_opts ?: ''`; command becomes

  ```
      repgenr ${opts} dereplicate-chunk \\
          --genomes-fofn genomes.fofn \\
          --out ${meta.id} \\
          ${args} \\
          \$sel \\
          --threads ${task.cpus} \\
          --versions-out tool_versions.yml
  ```

  plus `when:` and the stub echo.
- `derep_merge.nf`: same as chunk, and the output directory becomes the prefix: `def prefix = task.ext.prefix ?: "${meta.id}"` in both `script:` and `stub:`, `--out ${prefix}`, output `tuple val(meta), path("${prefix}"), emit: reps`, stub `mkdir -p ${prefix}/representatives` etc.; remove `publishDir`.
- `tree2tax.nf`: remove `publishDir`; add `when:`; `${opts}` and `${args}` in place of `${params.repgenr_opts}` and `${params.tree2tax_args}`; stub echo.

Note: the DEREP_MERGE output path is now `${prefix}`; the nf-test snapshot for `derep_merge` will change in Step 8 because the directory is named `merged` via `ext.prefix` (the test meta id is also `merged`, so the snapshot content may be identical; regenerate anyway).

- [ ] **Step 8: Run everything**

Run: `conda run -n repgenr_nf --no-capture-output nextflow lint nextflow/ 2>&1 | tail -3`
Expected: no warnings, no errors.

Run: `conda run -n repgenr_nf --no-capture-output nf-test test --tag stub --update-snapshot 2>&1 | grep -E "Executed|FAIL"`
Expected: `Executed 18 tests` (17 plus the new override test), 0 failed. Open each changed `.snap` and confirm only md5 lines for `versions.yml`/directories moved; commit the snapshots.

Run: `/Users/andreassjodin/miniforge3/envs/repgenr_dev/bin/python -m pytest -q tests/unit/test_nextflow_retry_window.py`
Expected: 2 passed.

Run: `conda run -n repgenr_nf --no-capture-output nextflow run nextflow/main.nf -profile test -stub --outdir /tmp/repgenr-stub 2>&1 | tail -5 && ls /tmp/repgenr-stub /tmp/repgenr-stub/phylo/tree`
Expected: completes; `tree2tax.tsv`, `genomes_map.tsv`, `metadata/`, `dereplicate/merged/`, `phylo/tree/tree.nwk`, `pipeline_info/software_versions.yml` present.

- [ ] **Step 9: Commit**

```bash
git add nextflow tests/unit/test_nextflow_retry_window.py
git commit -m "refactor(nextflow): conf/base.config and conf/modules.config; modules read task.ext only

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01H41vhDaRfyPXrcfFc1MXNG"
```

---

### Task 3: One shared versions helper instead of seven heredocs

**Files:**
- Create: `nextflow/bin/repgenr_versions_fragment` (executable)
- Modify: all seven modules (replace the heredoc block with one line)

**Interfaces:**
- Produces: `repgenr_versions_fragment <process-name> [tool_versions.yml]` writes `versions.yml` in the task directory. Nextflow puts `nextflow/bin/` (the directory next to `main.nf`) on the task PATH.

- [ ] **Step 1: Write the failing check**

The e2e test (`nextflow/tests/local_dataflow.nf.test`) already asserts `software_versions.yml` contains `repgenr`. Add to `nextflow/tests/bacterial_dataflow.nf.test`, inside `then`, so the stub path also checks the fragment shape produced by the helper (stubs `touch versions.yml`, so this asserts the file exists and the pipeline still collects it):

```groovy
            def versions = path("$outputDir/pipeline_info/software_versions.yml")
            assert versions.exists()
```

(This step keeps the suite honest; the real assertion on helper output is the e2e job in CI, which runs real tools.)

- [ ] **Step 2: Create the helper**

`nextflow/bin/repgenr_versions_fragment`:

```bash
#!/usr/bin/env bash
# Write versions.yml for the calling Nextflow process: repgenr's own version
# plus the external tool versions the stage recorded (a tool_versions.yml
# fragment written by --versions-out or `repgenr versions`).
#
# Usage: repgenr_versions_fragment "<task.process>" [tool_versions.yml]
set -eu

process="$1"
tools="${2:-tool_versions.yml}"

{
    printf '"%s":\n' "$process"
    printf '    repgenr: %s\n' "$(repgenr --version | sed 's/repgenr //')"
    if [ -f "$tools" ]; then
        cat "$tools"
    fi
} > versions.yml
```

Run: `chmod +x nextflow/bin/repgenr_versions_fragment && git add nextflow/bin/repgenr_versions_fragment && git ls-files -s nextflow/bin | cut -c1-6`
Expected: mode `100755`.

- [ ] **Step 3: Replace the heredocs**

In every module `script:` block replace

```
    cat > versions.yml <<END_VERSIONS
"${task.process}":
    repgenr: \$(repgenr --version | sed 's/repgenr //')
END_VERSIONS
    cat tool_versions.yml >> versions.yml
```

with

```
    repgenr_versions_fragment "${task.process}" tool_versions.yml
```

METADATA and VACQUIRE first run `repgenr versions -wd <wd> --versions-out tool_versions.yml` as they do today, then the helper line. Stubs keep `touch versions.yml`.

- [ ] **Step 4: Verify with real tools locally**

The helper only runs in non-stub mode. Run the e2e harness with the local tools (the `repgenr_dev` env has sourmash and mashtree; the `repgenr_nf` env has Nextflow):

```bash
cd /Users/andreassjodin/Code/repgenr
PATH=/Users/andreassjodin/miniforge3/envs/repgenr_dev/bin:$PATH conda run -n repgenr_nf --no-capture-output \
  nextflow run nextflow/tests/local_dataflow.nf -profile test \
  --genomes_dir nextflow/tests/data/e2e_genomes --derep_tool sourmash --derep_secondary_ani 0.95 \
  --derep_process_size 4 --outdir /tmp/repgenr-e2e 2>&1 | tail -3
cat /tmp/repgenr-e2e/pipeline_info/software_versions.yml
```

Expected: the file lists `repgenr: 3.0.0` under each process name and `sourmash`/`mashtree` versions under the processes that ran them.

Run: `conda run -n repgenr_nf --no-capture-output nf-test test --tag stub 2>&1 | grep -E "Executed|FAIL"`
Expected: 18 passed.

- [ ] **Step 5: Commit**

```bash
git add nextflow
git commit -m "refactor(nextflow): shared versions.yml helper in nextflow/bin

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01H41vhDaRfyPXrcfFc1MXNG"
```

---

### Task 4: run_meta and the meta on the acquisition front

**Files:**
- Create: `nextflow/subworkflows/local/run_meta.nf`, `nextflow/tests/run_meta.nf.test`
- Modify: `nextflow/main.nf`, `nextflow/modules/local/dataflow/{metadata,genome,vacquire}.nf`, `nextflow/subworkflows/local/{acquire,bacterial_dataflow,viral_dataflow}.nf`
- Modify: `nextflow/tests/{metadata_process,genome_process,vacquire_process,acquire}.nf.test` and the three snapshots; harnesses `nextflow/tests/{acquire_scatter,bacterial_dataflow,viral_dataflow}.nf`

**Interfaces:**
- Produces: `run_meta(Map params) -> Map` with keys `id` (String) and `mode` (String).
- Produces channel shapes: `METADATA(val(meta))` emits `selection: tuple(meta, path)`, `outgroup_accession: tuple(meta, path)`; `GENOME(tuple(meta, selection))` emits `genomes: tuple(meta, List<Path>)`, `outgroup: tuple(meta, List<Path>)` (optional); `VACQUIRE(val(meta))` emits `genomes`, `outgroup`, `outgroup_accession` in the same shapes; `ACQUIRE(ch_meta)` emits `genomes: tuple(meta, List<Path>)`, `outgroup: tuple(meta, List<Path>)` (empty list when none), `selection: tuple(meta, path)`, `outgroup_accession: tuple(meta, path)`.
- The dataflow subworkflows take `ch_meta` and, until Task 5 and Task 6, map these tuples down to the old bare shapes for the unchanged downstream processes.

- [ ] **Step 1: Write the failing function test**

`nextflow/tests/run_meta.nf.test`:

```groovy
nextflow_function {

    name "run_meta"
    script "../subworkflows/local/run_meta.nf"
    function "run_meta"
    tag "stub"

    test("bacterial genus target becomes the id") {
        when {
            function {
                """
                input[0] = [mode: 'bacterial',
                            metadata_args: '-r 232.0 --gtdb-version bac120 -d rep -l genus -tg Francisella',
                            vmetadata_args: '', vgenome_args: '']
                """
            }
        }
        then {
            assert function.success
            assert function.result == [id: 'francisella', mode: 'bacterial']
        }
    }

    test("species target wins over family, and is slugged") {
        when {
            function {
                """
                input[0] = [mode: 'bacterial',
                            metadata_args: '-l species -tf Francisellaceae -ts novicida-like',
                            vmetadata_args: '', vgenome_args: '']
                """
            }
        }
        then {
            assert function.success
            assert function.result.id == 'novicida_like'
        }
    }

    test("viral target comes from vgenome_args, then vmetadata_args") {
        when {
            function {
                """
                input[0] = [mode: 'viral', metadata_args: '',
                            vmetadata_args: '-t adenoviridae', vgenome_args: '-tg Mastadenovirus']
                """
            }
        }
        then {
            assert function.success
            assert function.result == [id: 'mastadenovirus', mode: 'viral']
        }
    }

    test("no target falls back to the mode") {
        when {
            function {
                """
                input[0] = [mode: 'bacterial', metadata_args: '', vmetadata_args: '', vgenome_args: '']
                """
            }
        }
        then {
            assert function.success
            assert function.result == [id: 'bacterial', mode: 'bacterial']
        }
    }
}
```

Run: `conda run -n repgenr_nf --no-capture-output nf-test test nextflow/tests/run_meta.nf.test 2>&1 | grep -E "PASSED|FAILED|Executed"`
Expected: fails (script missing).

- [ ] **Step 2: Implement run_meta**

`nextflow/subworkflows/local/run_meta.nf`:

```groovy
// The run-level meta map: one per pipeline run, carried on every channel.
//
// id   -- a slug of the selection target, taken from the argument strings the
//         user passes to the metadata stages, so published work directories
//         and task tags name the taxon. Falls back to the mode.
// mode -- 'bacterial' or 'viral'.

def target_after(String args, String flag) {
    def tokens = args ? args.tokenize() : []
    def i = tokens.indexOf(flag)
    return (i >= 0 && i + 1 < tokens.size()) ? tokens[i + 1] : null
}

def run_meta(Map params) {
    def target = null
    if (params.mode == 'viral') {
        target = target_after(params.vgenome_args ?: '', '-tg')
            ?: target_after(params.vmetadata_args ?: '', '-t')
    }
    else {
        def bact = params.metadata_args ?: ''
        target = target_after(bact, '-ts')
            ?: target_after(bact, '-tg')
            ?: target_after(bact, '-tf')
    }
    def id = (target ?: params.mode).toLowerCase().replaceAll('[^a-z0-9]+', '_')
    return [id: id, mode: params.mode]
}
```

Run the function test again. Expected: 4 passed.

- [ ] **Step 3: Update the process tests to the tuple shapes (failing first)**

`nextflow/tests/metadata_process.nf.test`, replace the test body:

```groovy
    test("emits selection.tsv, the outgroup accession, and versions") {
        when {
            process {
                """
                input[0] = [id: 'test', mode: 'bacterial']
                """
            }
        }
        then {
            assert process.success
            assert snapshot(process.out).match()
            with(process.out.selection.get(0)) {
                assert get(0).id == 'test'
                def lines = path(get(1)).readLines()
                assert lines[0] == 'accession\tfamily\tgenus\tspecies\tis_outgroup\tfilename\tcompleteness\tcontamination'
                assert lines.size() == 4
                assert lines.count { row -> row.split('\t')[4] == '1' } == 1
            }
            assert path(process.out.outgroup_accession.get(0).get(1)).text.trim() == 'GCF_000009.1'
        }
    }
```

`nextflow/tests/genome_process.nf.test`:

```groovy
    test("materializes every selection row as a genome or outgroup FASTA") {
        when {
            process {
                """
                input[0] = tuple([id: 'test', mode: 'bacterial'],
                                 file("${projectDir}/nextflow/tests/data/selection.tsv"))
                """
            }
        }
        then {
            assert process.success
            assert snapshot(process.out).match()
            with(process.out.genomes.get(0)) {
                assert get(0).id == 'test'
                def genomes = [get(1)].flatten().collect { p -> file(p).name }.sort()
                assert genomes == ['Fam_Gen_sp1_GCF_000001.1.fasta', 'Fam_Gen_sp2_GCF_000002.1.fasta']
            }
            with(process.out.outgroup.get(0)) {
                def outgroup = [get(1)].flatten().collect { p -> file(p).name }
                assert outgroup == ['Fam_Out_grp_GCF_000009.1.fasta']
            }
        }
    }
```

`nextflow/tests/vacquire_process.nf.test`:

```groovy
    test("emits viral genomes, the outgroup, and its accession") {
        when {
            process {
                """
                input[0] = [id: 'test', mode: 'viral']
                """
            }
        }
        then {
            assert process.success
            assert snapshot(process.out).match()
            with(process.out.genomes.get(0)) {
                assert get(0).id == 'test'
                def genomes = [get(1)].flatten().collect { p -> file(p).name }.sort()
                assert genomes == ['Vir_gen_sp1_iso1.fasta', 'Vir_gen_sp2_iso2.fasta']
            }
            def outgroup = [process.out.outgroup.get(0).get(1)].flatten().collect { p -> file(p).name }
            assert outgroup == ['Vir_out_grp_iso9.fasta']
            assert path(process.out.outgroup_accession.get(0).get(1)).text.trim() == 'iso9'
        }
    }
```

`nextflow/tests/acquire.nf.test`:

```groovy
    test("metadata selection drives the genome download") {
        when {
            workflow {
                """
                input[0] = channel.value([id: 'test', mode: 'bacterial'])
                """
            }
        }
        then {
            assert workflow.success
            with(workflow.out.genomes.get(0)) {
                assert get(0).id == 'test'
                assert get(1).size() == 2
            }
            with(workflow.out.outgroup.get(0)) {
                assert get(1).size() == 1
            }
            assert workflow.out.selection.get(0).get(0).id == 'test'
        }
    }
```

Run: `conda run -n repgenr_nf --no-capture-output nf-test test nextflow/tests/metadata_process.nf.test nextflow/tests/genome_process.nf.test nextflow/tests/vacquire_process.nf.test nextflow/tests/acquire.nf.test 2>&1 | grep -E "Executed|FAIL"`
Expected: 4 failed.

- [ ] **Step 4: Reshape the three front modules**

`metadata.nf` input/output blocks:

```groovy
    input:
    val meta

    output:
    tuple val(meta), path("selection.tsv")         , emit: selection
    tuple val(meta), path("outgroup_accession.txt"), emit: outgroup_accession
    path "versions.yml"                            , emit: versions
```

and `tag "${meta.id}"`.

`genome.nf`:

```groovy
    input:
    tuple val(meta), path(selection)

    output:
    tuple val(meta), path("out/genomes/*") , emit: genomes
    tuple val(meta), path("out/outgroup/*"), emit: outgroup, optional: true
    path "versions.yml"                    , emit: versions
```

and `tag "${meta.id}"`.

`vacquire.nf`:

```groovy
    input:
    val meta

    output:
    tuple val(meta), path("out/genomes/*")         , emit: genomes
    tuple val(meta), path("out/outgroup/*")        , emit: outgroup, optional: true
    tuple val(meta), path("outgroup_accession.txt"), emit: outgroup_accession
    path "versions.yml"                            , emit: versions
```

and `tag "${meta.id}"`.

- [ ] **Step 5: ACQUIRE with metas**

`nextflow/subworkflows/local/acquire.nf`:

```groovy
// Acquire genomes (data-channel form): metadata -> genome.
//
// METADATA selects accessions and emits selection.tsv; GENOME downloads them and
// emits the genome FASTAs. Every output carries the run meta, and the genome and
// outgroup lists are normalised to lists so DEREPLICATE_SCATTER can chunk them.

include { METADATA } from '../../modules/local/dataflow/metadata'
include { GENOME   } from '../../modules/local/dataflow/genome'

workflow ACQUIRE {
    take:
    ch_meta   // value channel: the run meta map

    main:
    def ch_versions = channel.empty()

    METADATA(ch_meta)
    ch_versions = ch_versions.mix(METADATA.out.versions)

    GENOME(METADATA.out.selection)
    ch_versions = ch_versions.mix(GENOME.out.versions)

    def ch_genomes = GENOME.out.genomes
        .map { meta, files -> tuple(meta, files instanceof List ? files : [files]) }

    // The outgroup output is optional: join with remainder so a run without one
    // still emits tuple(meta, []) and downstream joins on the meta keep working.
    def ch_outgroup = ch_genomes
        .map { meta, _files -> meta }
        .join(GENOME.out.outgroup, by: 0, remainder: true)
        .map { meta, files ->
            def list = files == null ? [] : (files instanceof List ? files : [files])
            tuple(meta, list)
        }

    emit:
    genomes            = ch_genomes
    outgroup           = ch_outgroup
    selection          = METADATA.out.selection
    outgroup_accession = METADATA.out.outgroup_accession
    versions           = ch_versions
}
```

- [ ] **Step 6: Thread ch_meta through main.nf and the dataflows, with temporary glue**

`nextflow/main.nf` workflow body after validation:

```groovy
    def ch_meta = channel.value(run_meta(params))

    def ch_versions = channel.empty()
    if (params.mode == 'viral') {
        VIRAL_DATAFLOW(ch_meta)
        ch_versions = VIRAL_DATAFLOW.out.versions
    }
    else {
        BACTERIAL_DATAFLOW(ch_meta)
        ch_versions = BACTERIAL_DATAFLOW.out.versions
    }

    PUBLISH_VERSIONS(ch_versions)
```

with `include { run_meta } from './subworkflows/local/run_meta'` among the includes.

`bacterial_dataflow.nf` (glue lines marked; Task 5 and Task 6 remove them):

```groovy
workflow BACTERIAL_DATAFLOW {
    take:
    ch_meta

    main:
    def ch_versions = channel.empty()

    ACQUIRE(ch_meta)
    ch_versions = ch_versions.mix(ACQUIRE.out.versions)

    // Glue until Task 5: the scatter still takes bare paths.
    def ch_genome_files = ACQUIRE.out.genomes.flatMap { _meta, files -> files }
    def ch_selection    = ACQUIRE.out.selection.map { _meta, sel -> sel }
    DEREPLICATE_SCATTER(ch_genome_files, ch_selection)
    ch_versions = ch_versions.mix(DEREPLICATE_SCATTER.out.versions)

    // Glue until Task 6: phylo and tree2tax still take bare paths.
    def ch_reps     = DEREPLICATE_SCATTER.out.reps.map { _meta, dir -> dir }
    def ch_outgroup = ACQUIRE.out.outgroup.map { _meta, files -> files }
    def ch_og_acc   = ACQUIRE.out.outgroup_accession.map { _meta, acc -> acc }

    PHYLO(ch_reps, ch_outgroup, ch_og_acc)
    ch_versions = ch_versions.mix(PHYLO.out.versions)

    TREE2TAX(PHYLO.out.tree, ch_reps, ch_outgroup, ch_og_acc)
    ch_versions = ch_versions.mix(TREE2TAX.out.versions)

    emit:
    tree        = PHYLO.out.tree
    tree2tax    = TREE2TAX.out.tree2tax
    genomes_map = TREE2TAX.out.genomes_map
    versions    = ch_versions
}
```

`viral_dataflow.nf` the same way: `take: ch_meta`, `VACQUIRE(ch_meta)`, glue `VACQUIRE.out.genomes.flatMap { _meta, files -> files instanceof List ? files : [files] }` into the scatter with `channel.value([])` as selection, and the outgroup glue `VACQUIRE.out.outgroup.map { _meta, files -> files }.ifEmpty([])`, accession glue `.map { _meta, acc -> acc }`.

Harnesses: `bacterial_dataflow.nf` and `viral_dataflow.nf` under `nextflow/tests/` include `run_meta` and pass `channel.value(run_meta(params))`; `acquire_scatter.nf` does `ACQUIRE(channel.value(run_meta(params)))` and applies the same two glue maps before `DEREPLICATE_SCATTER`.

- [ ] **Step 7: Run the suite, regenerate the three snapshots, lint**

Run: `conda run -n repgenr_nf --no-capture-output nf-test test --tag stub --update-snapshot 2>&1 | grep -E "Executed|FAIL"`
Expected: 22 tests (18 plus 4 function tests), 0 failed. Inspect `metadata_process`, `genome_process`, `vacquire_process` snapshots: each output entry now has a meta map before the path.

Run: `conda run -n repgenr_nf --no-capture-output nextflow lint nextflow/ 2>&1 | tail -2`
Expected: no warnings.

- [ ] **Step 8: Commit**

```bash
git add nextflow
git commit -m "feat(nextflow): run-level meta from params; meta on the acquisition front

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01H41vhDaRfyPXrcfFc1MXNG"
```

---

### Task 5: Meta through the scatter-gather dereplication

**Files:**
- Modify: `nextflow/subworkflows/local/dereplicate_scatter.nf`, `nextflow/modules/local/{derep_chunk,derep_merge}.nf` (tags only), `nextflow/subworkflows/local/{bacterial_dataflow,viral_dataflow}.nf` (remove the Task 5 glue), harnesses `nextflow/tests/{acquire_scatter,dereplicate_scatter,local_dataflow}.nf`
- Modify: `nextflow/tests/{dereplicate_scatter,derep_chunk,derep_merge,acquire_scatter}.nf.test`

**Interfaces:**
- Consumes: `ACQUIRE.out.genomes: tuple(meta, List<Path>)`, `ACQUIRE.out.selection: tuple(meta, path)`.
- Produces: `DEREPLICATE_SCATTER(ch_genomes, ch_selection)` takes `tuple(meta, List<Path>)` and `tuple(meta, path_or_[])`, emits `reps: tuple(meta, dir)` with the run meta. Chunk metas are `[id: "${meta.id}.chunk_${i}", mode: meta.mode, run: meta]`; the merge meta is `[id: "${meta.id}.merged", mode: meta.mode, run: meta]`. DEREP_CHUNK and DEREP_MERGE keep `path selection, stageAs: 'selection.tsv'` as a bare auxiliary input (nf-core reference-file style).

- [ ] **Step 1: Update the scatter test (failing first)**

`nextflow/tests/dereplicate_scatter.nf.test` `when`/`then`:

```groovy
            workflow {
                """
                input[0] = channel.fromPath("${projectDir}/nextflow/tests/data/genomes/*.fasta")
                    .collect()
                    .map { files -> tuple([id: 'test', mode: 'bacterial'], files) }
                input[1] = channel.value(tuple([id: 'test', mode: 'bacterial'], []))
                """
            }
        }
        then {
            assert workflow.success
            // exactly one merged result, carrying the run meta again
            assert workflow.out.reps.size() == 1
            with(workflow.out.reps.get(0)) {
                assert get(0) == [id: 'test', mode: 'bacterial']
                def dir = get(1)
                assert path("${dir}/clusters.tsv").exists()
                assert path("${dir}/representatives").exists()
            }
            assert workflow.out.versions
```

In `derep_chunk.nf.test` and `derep_merge.nf.test` change the metas to `[id: 'test.chunk_0', mode: 'bacterial', run: [id: 'test', mode: 'bacterial']]` and `[id: 'test.merged', mode: 'bacterial', run: [id: 'test', mode: 'bacterial']]`, and the `get(0).id` assertions to match. In `derep_merge.nf.test` the output directory is `merged` (from `ext.prefix` in `conf/modules.config`), so `path("${dir}/...")` assertions still hold.

Run: `conda run -n repgenr_nf --no-capture-output nf-test test nextflow/tests/dereplicate_scatter.nf.test 2>&1 | grep -E "Executed|FAIL"`
Expected: 1 failed.

- [ ] **Step 2: Rewrite the scatter subworkflow**

`nextflow/subworkflows/local/dereplicate_scatter.nf`:

```groovy
// Scatter-gather dereplication (data-channel style).
//
// Genomes arrive as one tuple(meta, [fasta...]) per run. They are grouped into
// chunks of `derep_process_size`, each chunk is dereplicated independently and
// in parallel (scatter, one task per chunk -- on HPC these land on separate
// nodes), and the union of the chunk representatives is dereplicated once more
// (gather) to produce the final representative set. This is the two-stage
// reduce-tree of the in-process dereplicate stage, expressed as typed channels
// so Nextflow owns the fan-out instead of Python threads.
//
// Chunk and merge metas nest the run meta under `run`, so the result can be
// emitted under the run meta again and joined with the other run-level channels.

include { DEREP_CHUNK } from '../../modules/local/derep_chunk'
include { DEREP_MERGE } from '../../modules/local/derep_merge'

workflow DEREPLICATE_SCATTER {
    take:
    ch_genomes    // channel: tuple(meta, [genome FASTA paths])
    ch_selection  // channel: tuple(meta, selection.tsv) or tuple(meta, [])

    main:
    def ch_versions = channel.empty()

    // A null/zero process size means a single chunk. Coerce because
    // command-line params arrive as strings.
    def requested_size = params.derep_process_size ? (params.derep_process_size as Integer) : 0
    def chunk_size = requested_size > 0 ? requested_size : 1000000

    def ch_chunks = ch_genomes.flatMap { meta, files ->
        files.collate(chunk_size).withIndex().collect { chunk, i ->
            tuple([id: "${meta.id}.chunk_${i}", mode: meta.mode, run: meta], chunk)
        }
    }

    // The selection file is an auxiliary input shared by every chunk task, so
    // it is passed as a bare value channel (nf-core reference-file style).
    def ch_selection_file = ch_selection.map { _meta, sel -> sel }.first()

    DEREP_CHUNK(ch_chunks, ch_selection_file)
    ch_versions = ch_versions.mix(DEREP_CHUNK.out.versions.first())

    // Gather every chunk directory under one merge meta per run.
    def ch_merge_in = DEREP_CHUNK.out.chunk
        .map { meta, dir -> tuple([id: "${meta.run.id}.merged", mode: meta.mode, run: meta.run], dir) }
        .groupTuple()

    DEREP_MERGE(ch_merge_in, ch_selection_file)
    ch_versions = ch_versions.mix(DEREP_MERGE.out.versions)

    emit:
    reps     = DEREP_MERGE.out.reps.map { meta, dir -> tuple(meta.run, dir) }
    versions = ch_versions
}
```

DEREP_CHUNK and DEREP_MERGE need no input changes beyond Task 2; confirm both have `tag "${meta.id}"`.

- [ ] **Step 3: Remove the Task 5 glue in the dataflows and harnesses**

In `bacterial_dataflow.nf` replace the two glue lines and the call with:

```groovy
    DEREPLICATE_SCATTER(ACQUIRE.out.genomes, ACQUIRE.out.selection)
```

In `viral_dataflow.nf`:

```groovy
    def ch_genomes = VACQUIRE.out.genomes
        .map { meta, files -> tuple(meta, files instanceof List ? files : [files]) }
    // vmetadata writes no selection.tsv and viral genomes carry no CheckM
    // quality, so the keeper input is an empty list under the run meta.
    def ch_no_selection = ch_meta.map { meta -> tuple(meta, []) }
    DEREPLICATE_SCATTER(ch_genomes, ch_no_selection)
```

`nextflow/tests/acquire_scatter.nf`:

```groovy
workflow {
    def ch_meta = channel.value(run_meta(params))
    ACQUIRE(ch_meta)
    DEREPLICATE_SCATTER(ACQUIRE.out.genomes, ACQUIRE.out.selection)

    DEREPLICATE_SCATTER.out.reps
        .map { _meta, dir -> dir }
        .collectFile(name: 'merged_path.txt', storeDir: params.outdir) { dir -> "${dir}\n" }
}
```

`nextflow/tests/dereplicate_scatter.nf` and `nextflow/tests/local_dataflow.nf`:

```groovy
    def meta = [id: 'local', mode: 'bacterial']
    def ch_genomes = channel
        .fromPath("${params.genomes_dir}/*.{fasta,fa,fna,fas}")
        .filter { f -> !f.name.startsWith('._') }   // skip macOS AppleDouble files
        .collect()
        .map { files -> tuple(meta, files) }

    DEREPLICATE_SCATTER(ch_genomes, channel.value(tuple(meta, [])))
```

(`local_dataflow.nf` keeps its PHYLO/TREE2TAX glue until Task 6: `def ch_reps = DEREPLICATE_SCATTER.out.reps.map { _meta, dir -> dir }`.)

- [ ] **Step 4: Run the suite and lint**

Run: `conda run -n repgenr_nf --no-capture-output nf-test test --tag stub --update-snapshot 2>&1 | grep -E "Executed|FAIL"`
Expected: 22 tests, 0 failed. `acquire_scatter.nf.test` still counts two DEREP_CHUNK tasks.

Run: `conda run -n repgenr_nf --no-capture-output nextflow lint nextflow/ 2>&1 | tail -2`
Expected: no warnings.

- [ ] **Step 5: Commit**

```bash
git add nextflow
git commit -m "feat(nextflow): chunk and merge metas derived from the run meta

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01H41vhDaRfyPXrcfFc1MXNG"
```

---

### Task 6: Meta through PHYLO and TREE2TAX; join on meta; end-to-end id test

**Files:**
- Modify: `nextflow/modules/local/dataflow/{phylo,tree2tax}.nf`, `nextflow/subworkflows/local/{bacterial_dataflow,viral_dataflow}.nf`, `nextflow/tests/local_dataflow.nf`
- Modify: `nextflow/tests/{phylo_process,tree2tax_process,modules_config,bacterial_dataflow}.nf.test` and two snapshots
- Modify: `docs/superpowers/specs/2026-09-05-nfcore-meta-ext-args-design.md` (record the single-tuple refinement)

**Interfaces:**
- Consumes: `DEREPLICATE_SCATTER.out.reps: tuple(meta, dir)`, `ACQUIRE.out.outgroup: tuple(meta, List<Path>)`, `ACQUIRE.out.outgroup_accession: tuple(meta, path)`.
- Produces: `PHYLO(tuple(meta, reps_dir, outgroup_list, outgroup_accession))` emits `tree: tuple(meta, path)`; `TREE2TAX(tuple(meta, tree, reps_dir, outgroup_list, outgroup_accession))` emits `tree2tax: tuple(meta, path)`, `genomes_map: tuple(meta, path)`.

- [ ] **Step 1: Update the process tests (failing first)**

`phylo_process.nf.test` input and assertions:

```groovy
                input[0] = tuple(
                    [id: 'test', mode: 'bacterial'],
                    file("${projectDir}/nextflow/tests/data/derep_out"),
                    [ file("${projectDir}/nextflow/tests/data/outgroup/Fam_Out_grp_GCF_000009.1.fasta") ],
                    file("${projectDir}/nextflow/tests/data/outgroup_accession.txt")
                )
```

```groovy
            with(process.out.tree.get(0)) {
                assert get(0).id == 'test'
                def newick = path(get(1)).text.trim()
                assert newick.endsWith(';')
                assert newick.contains('Fam_Gen_sp1_GCF_000001.1')
                assert newick.contains('Fam_Gen_sp2_GCF_000002.1')
            }
```

`tree2tax_process.nf.test`:

```groovy
                input[0] = tuple(
                    [id: 'test', mode: 'bacterial'],
                    file("${projectDir}/nextflow/tests/data/tree.nwk"),
                    file("${projectDir}/nextflow/tests/data/derep_out"),
                    [ file("${projectDir}/nextflow/tests/data/outgroup/Fam_Out_grp_GCF_000009.1.fasta") ],
                    file("${projectDir}/nextflow/tests/data/outgroup_accession.txt")
                )
```

```groovy
            with(process.out.tree2tax.get(0)) {
                assert get(0).id == 'test'
                def relations = path(get(1)).readLines()
                assert relations[0] == 'child\tparent'
                assert relations.size() == 3
            }
            def genome_map = path(process.out.genomes_map.get(0).get(1)).readLines()
            assert genome_map.size() == 2
            assert genome_map.every { row -> row.contains('\t') }
```

`modules_config.nf.test` input becomes the same single tuple as `phylo_process.nf.test`.

`bacterial_dataflow.nf.test`, add inside `then`:

```groovy
            // The run meta id built from metadata_args under the test profile
            // ('-tg francisella') tags every task through to TREE2TAX.
            def names = workflow.trace.succeeded().collect { t -> t.name }
            assert names.any { n -> n.contains('TREE2TAX (francisella)') }
            assert names.any { n -> n.contains('DEREP_MERGE (francisella.merged)') }
```

Run: `conda run -n repgenr_nf --no-capture-output nf-test test nextflow/tests/phylo_process.nf.test nextflow/tests/tree2tax_process.nf.test 2>&1 | grep -E "Executed|FAIL"`
Expected: 2 failed.

- [ ] **Step 2: Reshape PHYLO and TREE2TAX**

`phylo.nf` input/output:

```groovy
    input:
    tuple val(meta), path(reps_dir), path(outgroup, stageAs: 'outgroup/*'), path(outgroup_accession)

    output:
    tuple val(meta), path("tree/tree.nwk"), emit: tree
    path "versions.yml"                   , emit: versions
```

and `tag "${meta.id}"`. `tree2tax.nf`:

```groovy
    input:
    tuple val(meta), path(tree), path(reps_dir), path(outgroup, stageAs: 'outgroup/*'), path(outgroup_accession)

    output:
    tuple val(meta), path("tree2tax.tsv")   , emit: tree2tax
    tuple val(meta), path("genomes_map.tsv"), emit: genomes_map
    path "versions.yml"                     , emit: versions
```

and `tag "${meta.id}"`. Scripts and stubs are unchanged (they already reference `reps_dir`, `outgroup`, `outgroup_accession`, `tree`).

- [ ] **Step 3: Join on meta in the dataflows and remove the last glue**

`bacterial_dataflow.nf` after the scatter:

```groovy
    def ch_phylo_in = DEREPLICATE_SCATTER.out.reps
        .join(ACQUIRE.out.outgroup, by: 0)
        .join(ACQUIRE.out.outgroup_accession, by: 0)

    PHYLO(ch_phylo_in)
    ch_versions = ch_versions.mix(PHYLO.out.versions)

    // [meta, tree] joined with [meta, reps, outgroup, accession]
    def ch_tree2tax_in = PHYLO.out.tree.join(ch_phylo_in, by: 0)

    TREE2TAX(ch_tree2tax_in)
    ch_versions = ch_versions.mix(TREE2TAX.out.versions)
```

`viral_dataflow.nf` the same, with VACQUIRE's outgroup normalised first:

```groovy
    def ch_outgroup = ch_genomes
        .map { meta, _files -> meta }
        .join(VACQUIRE.out.outgroup, by: 0, remainder: true)
        .map { meta, files ->
            def list = files == null ? [] : (files instanceof List ? files : [files])
            tuple(meta, list)
        }
    def ch_phylo_in = DEREPLICATE_SCATTER.out.reps
        .join(ch_outgroup, by: 0)
        .join(VACQUIRE.out.outgroup_accession, by: 0)
```

`nextflow/tests/local_dataflow.nf`:

```groovy
    def ch_phylo_in = DEREPLICATE_SCATTER.out.reps
        .map { m, dir -> tuple(m, dir, [], file(params.empty_accession)) }
    PHYLO(ch_phylo_in)
    TREE2TAX(PHYLO.out.tree.join(ch_phylo_in, by: 0))
```

- [ ] **Step 4: Amend the spec to match**

In the spec's Section 2, replace the sentence beginning "`PHYLO(tuple(meta, reps_dir), tuple(meta, [outgroup...]), ...`" with:

```
- `PHYLO(tuple(meta, reps_dir, [outgroup...], outgroup_accession))` and
  `TREE2TAX(tuple(meta, tree, reps_dir, [outgroup...], outgroup_accession))`:
  one joined tuple per process, built in the dataflow subworkflows with
  `.join(by: 0)` on the meta, which is the nf-core shape for a process with
  several per-run inputs. The dereplication processes keep `selection.tsv` as a
  bare auxiliary input (reference-file style), and DEREP_MERGE's output
  directory is `task.ext.prefix` (set to `merged` in `conf/modules.config`) so
  published paths do not change.
```

- [ ] **Step 5: Run everything**

Run: `conda run -n repgenr_nf --no-capture-output nf-test test --tag stub --update-snapshot 2>&1 | grep -E "Executed|FAIL"`
Expected: 22 tests, 0 failed. Inspect the `phylo_process` and `tree2tax_process` snapshots.

Run: `conda run -n repgenr_nf --no-capture-output nextflow lint nextflow/ 2>&1 | tail -2`
Expected: no warnings.

Run the real e2e locally (same command as Task 3 Step 4). Expected: 3 representatives, 6 mapped genomes, `software_versions.yml` present.

Run: `conda run -n repgenr_nf --no-capture-output nextflow run nextflow/main.nf -profile test -stub --outdir /tmp/repgenr-stub2 2>&1 | grep -E "TREE2TAX|Completed"`
Expected: task line shows `TREE2TAX (francisella)`; outputs under `/tmp/repgenr-stub2` as in Task 2.

- [ ] **Step 6: Commit**

```bash
git add nextflow docs/superpowers/specs/2026-09-05-nfcore-meta-ext-args-design.md
git commit -m "feat(nextflow): phylo and tree2tax take one meta tuple; dataflows join on meta

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01H41vhDaRfyPXrcfFc1MXNG"
```

---

### Task 7: Documentation, schema descriptions, changelog

**Files:**
- Modify: `docs/usage.md` (new section after "### Scatter-gather dereplication"), `docs/architecture.md:65`, `README.md:142-158`, `nextflow/nextflow_schema.json` (five description strings), `CHANGELOG.md` (Unreleased)

- [ ] **Step 1: usage.md**

Insert before `### Pipeline structure`:

```markdown
### Configuring processes

Each Nextflow process reads its tool flags from `task.ext.args`, which
`nextflow/conf/modules.config` maps from the user-facing parameters:
`--metadata_args`, `--vmetadata_args` and `--vgenome_args`, `--phylo_args`,
`--tree2tax_args`, and the six `--derep_*` parameters (composed into one
string for the two dereplication processes). Publishing directories live in
the same file, and resources and the retry window in
`nextflow/conf/base.config`. A site can retune one process without touching
the pipeline by passing its own config:

```groovy
// site.config
process {
    withName: 'PHYLO' {
        ext.args = '--treebuilder iqtree --bootstrap 1000'
    }
    withName: 'DEREP_CHUNK|DEREP_MERGE' {
        ext.args = '--tool galah --secondary-ani 0.98'
    }
}
```

`nextflow run nextflow/main.nf -c site.config ...`. `ext.repgenr_opts` is
the hook for the top-level repgenr options every process prepends (the
container profiles set it through `--repgenr_opts`).

Every channel carries a meta map built once per run: `id` is a slug of the
selection target (`francisella` for `-tg francisella`), `mode` is
`bacterial` or `viral`. Task tags and the dereplication chunk names use it.
The pipeline requires Nextflow 26.04 or later (`nextflowVersion =
'!>=26.04.0'`).
```

- [ ] **Step 2: architecture.md**

After the line starting "* Nextflow provides the actual parallelism" add:

```markdown
* Every Nextflow channel carries a run-level meta map (`id`, `mode`) built
  from the parameters in `main.nf`; processes exchange
  `tuple val(meta), path(...)`, read tool flags only from `task.ext.args`
  (mapped from the parameters in `nextflow/conf/modules.config`), and are
  published from configuration rather than from the module files.
```

- [ ] **Step 3: README**

Append to the Nextflow section's last paragraph:

```markdown
Nextflow 26.04 or later is required. Per-process tool flags and publishing
are configured in `nextflow/conf/modules.config` (see `docs/usage.md`,
"Configuring processes").
```

- [ ] **Step 4: Schema descriptions**

In `nextflow/nextflow_schema.json` append to each of the five `*_args` descriptions the sentence: `Reaches the process as task.ext.args; a site config can override it with withName.` For example:

```json
"description": "Arguments for the phylogeny stage (aligner or tree builder). Reaches the process as task.ext.args; a site config can override it with withName."
```

Run: `conda run -n repgenr_nf --no-capture-output nf-test test nextflow/tests/params_validation.nf.test nextflow/tests/pipeline_help.nf.test 2>&1 | grep -E "Executed|FAIL"`
Expected: 3 passed (the schema still parses).

- [ ] **Step 5: Changelog**

Under `## [Unreleased]`, add to `### Changed`:

```markdown
- **Nextflow layer in nf-core shape.** Every channel carries a run-level meta
  map (`id` from the selection target, `mode`), processes exchange
  `tuple val(meta), path(...)`, tool flags reach processes as
  `task.ext.args` mapped from the unchanged parameters in
  `nextflow/conf/modules.config` (where publishing now lives too), and
  resources and the retry window sit in `nextflow/conf/base.config`. The
  Nextflow floor is 26.04 (`!>=26.04.0`) with nf-schema 2.6.1, and the layer
  lints clean under the strict parser. Command lines and parameters are
  unchanged; code that included a RepGenR module or subworkflow directly must
  adapt to the tuple shapes.
```

and under a `### Notes` heading (create it after `### Fixed` in Unreleased):

```markdown
- With the 26.04 floor, moving version reporting to `eval()` outputs on the
  `versions` topic costs nothing in compatibility and is the natural next step
  for the Nextflow layer.
```

- [ ] **Step 6: Commit**

```bash
git add docs README.md CHANGELOG.md nextflow/nextflow_schema.json
git commit -m "docs(nextflow): configuring processes, meta map, 26.04 floor

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01H41vhDaRfyPXrcfFc1MXNG"
```

---

### Task 8: Verification and pull request

**Files:** none new.

- [ ] **Step 1: The spec's verification list**

```bash
cd /Users/andreassjodin/Code/repgenr
conda run -n repgenr_nf --no-capture-output nextflow lint nextflow/ 2>&1 | tail -2
conda run -n repgenr_nf --no-capture-output nf-test test --tag stub 2>&1 | grep -E "Executed|FAIL"
conda run -n repgenr_nf --no-capture-output nextflow run nextflow/main.nf -profile test -stub --outdir /tmp/repgenr-final 2>&1 | tail -3
ls /tmp/repgenr-final /tmp/repgenr-final/phylo/tree /tmp/repgenr-final/pipeline_info
printf "process { withName: 'PHYLO' { ext.args = '--treebuilder sourmash' } }\n" > /tmp/override.config
conda run -n repgenr_nf --no-capture-output nextflow run nextflow/main.nf -profile test -stub -c /tmp/override.config --outdir /tmp/repgenr-override 2>&1 | grep -E "Completed|ERROR"
grep -rl "treebuilder sourmash" .nextflow/ 2>/dev/null | head -1 || find . -path '*/work/*/.command.sh' -newer /tmp/override.config -exec grep -l "ext.args: --treebuilder sourmash" {} + | head -1
LC_ALL=C grep -rlP '[^\x00-\x7F]' nextflow/ --include='*.nf' --include='*.config' --include='*.json' --include='*.test' | wc -l
/Users/andreassjodin/miniforge3/envs/repgenr_dev/bin/python -m pytest -q
/Users/andreassjodin/miniforge3/envs/repgenr_dev/bin/ruff check src/ tests/
/Users/andreassjodin/miniforge3/envs/repgenr_dev/bin/ruff format --check src/ tests/ benchmarks/
/Users/andreassjodin/miniforge3/envs/repgenr_dev/bin/mypy src/repgenr
```

Expected: no lint warnings; 22 stub tests pass; the stub run publishes `tree2tax.tsv`, `genomes_map.tsv`, `phylo/tree/tree.nwk`, `pipeline_info/software_versions.yml`; the override run's `.command.sh` (found via the `work/` search) contains `ext.args: --treebuilder sourmash`; ASCII count 0; Python gates green.

- [ ] **Step 2: Push and open the PR**

```bash
git push -u origin feat/nfcore-modules
gh pr create --base main --title "Nextflow layer in nf-core shape: meta maps, ext.args, 26.04 floor" --body-file - <<'EOF'
## Summary

Implements `docs/superpowers/specs/2026-09-05-nfcore-meta-ext-args-design.md`.

- Run-level meta map (`id` from the selection target, `mode`) built once in `main.nf`; every process exchanges `tuple val(meta), path(...)`; chunk and merge metas nest the run meta.
- `conf/modules.config` maps the unchanged `*_args` and `derep_*` parameters to `task.ext.args` per process and holds `publishDir`; `conf/base.config` holds resources and the retry window; modules read only `task.ext.*`.
- One shared `nextflow/bin/repgenr_versions_fragment` replaces seven heredocs; version reporting stays file-style.
- Floor `!>=26.04.0`, nf-schema 2.6.1, `nextflow lint` clean (was 21 warnings).
- Docs: "Configuring processes" in usage.md, meta paragraph in architecture.md, README floor, schema descriptions, changelog.

## Test plan

- [ ] `nextflow lint nextflow/`: 0 errors, 0 warnings
- [ ] nf-test stub suite (22 tests incl. new `modules_config` and `run_meta` tests) green locally and in CI
- [ ] Real e2e (`local_dataflow.nf`) green locally with sourmash + mashtree and in CI
- [ ] `-c` override of `PHYLO` `ext.args` reaches `.command.sh`
- [ ] Python gates unchanged and green

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01H41vhDaRfyPXrcfFc1MXNG
EOF
```

Tick the test-plan boxes from the Step 1 output before requesting review. Merge on green as for the previous PRs in this series.

---

## Self-review

**Spec coverage.** Section 1 (config split, base/modules config, resourceLimits in test, ext.args table, ext.repgenr_opts, meta rule, schema sentence): Tasks 2, 4, 7. Section 2 (tuple shapes, tag, when, task.ext reads, no params/publishDir, shared versions helper, subworkflow joins, harnesses, lint cleanup): Tasks 1 to 6. Section 3 (tests moved to tuples, snapshots regenerated, modules_config test, meta-id test, retry-window test, docs, changelog, compatibility note): Tasks 2, 4, 5, 6, 7. Floor and plugin: Task 1. Verification list: Task 8. Spec refinement (single joined tuple for PHYLO/TREE2TAX, bare selection input, ext.prefix for the merge directory) is recorded in Task 6 Step 4.

**Placeholders.** None; every code step shows the code.

**Type consistency.** `run_meta(Map params)` returns `[id, mode]` in Task 4 and is consumed in Tasks 4 to 6 and the harnesses. Chunk meta `[id, mode, run]` and merge meta `[id, mode, run]` are defined in Task 5 and asserted in the Task 5 tests and the Task 6 `DEREP_MERGE (francisella.merged)` tag. `task.ext.prefix` defaults to `"${meta.id}"` in Task 2 and is set to `merged` in `conf/modules.config` in Task 2, which the Task 5 and e2e assertions on `dereplicate/merged/representatives` rely on. PHYLO and TREE2TAX tuple orders `(meta, reps_dir, outgroup, outgroup_accession)` and `(meta, tree, reps_dir, outgroup, outgroup_accession)` match between the modules, the dataflow joins, the harness and the tests in Task 6.
