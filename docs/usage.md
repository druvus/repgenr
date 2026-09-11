# repgenr: Usage

## Introduction

RepGenR builds a representative-genome repository for a taxonomic group and the
matching FlexTaxD taxonomy. The pipeline runs five stages in order:

1. **metadata** -- select genome accessions (GTDB for bacteria/archaea; the
   viral path uses `vmetadata`/`vgenome`, sourcing from NCBI Virus by default).
2. **genome** -- download and organise the selected genomes with NCBI Datasets.
3. **dereplicate** -- cluster genomes by ANI and pick representatives.
4. **phylo** -- build a phylogeny from the representatives.
5. **tree2tax** -- emit a FlexTaxD-compatible taxonomy from the tree.

The CLI can also start from genomes already on disk: `repgenr ingest` stands
in for the first two stages (see "Starting from local genomes" below).

The Nextflow layer runs these stages as typed data channels: each stage emits its
outputs (the metadata selection, genome FASTAs, per-chunk and merged
representatives, the tree, the taxonomy) as staged files that the next stage
consumes. There is no shared working directory; results are published under
`--outdir`. Nextflow owns the fan-out (scatter-gather dereplication).

## Quick start

```bash
nextflow run nextflow/main.nf \
    --mode bacterial \
    --metadata_args '-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella' \
    --outdir results \
    -profile standard
```

Run `nextflow run nextflow/main.nf --help` for the parameter summary.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Published results and execution reports. |
| `--mode` | `bacterial` | `bacterial` (GTDB) or `viral` (NCBI Virus; BV-BRC via `--vmetadata_args "--source bvbrc"`). |
| `--metadata_args` | see config | Arguments for the bacterial metadata stage. |
| `--vmetadata_args` / `--vgenome_args` | see config | Viral metadata / genome selection arguments. |
| `--derep_tool` | `skder` | Dereplicator for the scatter-gather step. |
| `--derep_process_size` | `null` | Genomes per dereplication chunk (single chunk if unset). |
| `--derep_primary_ani` / `--derep_secondary_ani` / `--derep_aligned_fraction` | `0.90` / `0.99` / `0.50` | ANI / aligned-fraction thresholds. |
| `--phylo_args` | `--treebuilder mashtree` | Aligner or tree builder for the phylogeny. |
| `--phylo_split_msa` | `false` | Run the alignment and the tree as separate tasks. |
| `--tree2tax_args` | (empty) | tree-to-taxonomy (FlexTaxD) arguments; redundant genomes are listed by default (`--no-include-dereplicated` to omit them). |

`--phylo_split_msa` splits the phylogeny into `PHYLO_MSA` and `PHYLO_TREE`.
The alignment then keeps its own cache entry, so trying another tree builder or
another bootstrap re-runs only the tree, and the two halves take their own
resource labels: the alignment is `process_high`, the tree `process_medium`.
It applies to tree builders that consume an alignment; leave it off for
mashtree or sourmash, which build from genomes.

Parameters are validated against `nextflow/nextflow_schema.json` at launch.

### Passing parameters

Either pass them on the command line (`--name value`) or, for a reusable and
typed configuration, supply a params file:

```bash
nextflow run nextflow/main.nf -profile standard -params-file params.json
```

A starting point is in `nextflow/assets/params_example.json`. On the command line, give a string parameter whose value starts with a dash
as `--key=value` (for example `--tree2tax_args=--include-dereplicated`):
Nextflow reads `--key --flag` as the boolean `key = true`. A params file keeps
numeric values typed (e.g. `"derep_process_size": 2000`); on the command line
they arrive as strings, which the schema also accepts for the numeric options.

## Profiles

The `standard` profile (the default, local executor) caps every process's
CPU and memory request to what the machine has; `base.config` asks for up to
32 CPUs and 128 GB for the heavy processes, which the local executor would
otherwise refuse. `slurm` and `cloud` leave the requests as they are.

Combine an executor profile with an optional container profile, e.g.
`-profile slurm,singularity`.

- **Executors**: `standard` (local), `slurm`, `cloud` (AWS Batch).
- **Containers**: `docker`, `singularity`, `wave`. These set RepGenR's own
  adapter-level container backend (`--container ...`), which runs each external
  tool in a pinned image. RepGenR itself must be available to the Nextflow
  process.
- **`test`**: minimal resources and a small target for a quick smoke run.

## Scaling

For large dereplication inputs (10k+ genomes), set a chunk size so the
dereplication scatters across tasks (one per chunk):

```bash
nextflow run nextflow/main.nf --outdir results \
    --derep_tool sourmash --derep_process_size 2000 -profile slurm
```

Resource labels (`process_low/medium/high`) scale memory and time with the retry
attempt, so a task killed for memory or time is resubmitted with more headroom.
Tune the label values per environment in `nextflow.config`.

### Scatter-gather dereplication

For horizontal scaling, the `DEREPLICATE_SCATTER` subworkflow groups genomes into
chunks of `--derep_process_size`, dereplicates each chunk as a separate task (one
per node on HPC), and dereplicates the union of the chunk representatives once
more (the two-stage reduce-tree expressed as typed data channels). It wraps the
`repgenr dereplicate-chunk` / `dereplicate-merge` CLI steps and is driven by
`--derep_tool` and the `--derep_*_ani` / `--derep_aligned_fraction` thresholds.

A standalone harness runs it on a directory of genome FASTAs (no GTDB/NCBI
front-end), useful for testing and for dereplicating a local collection:

```bash
nextflow run nextflow/tests/dereplicate_scatter.nf -c nextflow/nextflow.config \
    --genomes_dir <DIR> --derep_tool sourmash --derep_process_size 2000 \
    --outdir results -profile standard
```

Add `-stub` to exercise the wiring without running the tools.

### Configuring processes

Each Nextflow process reads its tool flags from `task.ext.args`, which
`nextflow/conf/modules.config` maps from the user-facing parameters:
`--metadata_args`, `--vmetadata_args` and `--vgenome_args`, `--phylo_args`,
`--tree2tax_args`, and five dereplication parameters (`--derep_tool`,
`--derep_primary_ani`, `--derep_secondary_ani`, `--derep_aligned_fraction`,
`--derep_keeper`) composed into one string for the two dereplication
processes; `--derep_process_size` is read by the scatter subworkflow to
size the chunks and is not a tool flag. Publishing directories live in
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

An `ext.args` override replaces the whole string for that process; for the
dereplication processes that means every flag not repeated in the override
falls back to the repgenr CLI default rather than to the `--derep_*`
parameter.

Every channel carries a meta map built once per run: `id` is a slug of the
selection target (`francisella` for `-tg francisella`), `mode` is
`bacterial` or `viral`. Task tags and the dereplication chunk names use it.

Nextflow resolves both `bin/` and the default config relative to the
launched script. The harness scripts under `nextflow/tests/` therefore reach
the versions helper through the symlink
`nextflow/tests/bin/repgenr_versions_fragment`, and running one of them
directly needs `-c nextflow/nextflow.config` (nf-test supplies the config
itself). Any future script added to `nextflow/bin/` needs its own link under
`nextflow/tests/bin/`.

The pipeline requires Nextflow 26.04 or later (`nextflowVersion =
'!>=26.04.0'`).

### Pipeline structure

`nextflow/main.nf` dispatches by `--mode` to one of two data-channel subworkflows
that share the dereplication, phylo and tree2tax modules:

```
bacterial: ACQUIRE  (metadata -> genome)  -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
viral:     VACQUIRE (vmetadata -> vgenome) -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
```

`metadata` emits a portable `selection.tsv`; `genome-fetch` downloads the genomes
and emits them as a channel feeding the scatter-gather dereplication; `phylo` and
`tree2tax` run in task-local working directories and emit `tree.nwk`,
`tree2tax.tsv` and `genomes_map.tsv` to `--outdir`. Add `-stub` to any run for a
quick wiring check without external tools.

### Starting from local genomes

`repgenr ingest -wd WD --genomes-dir DIR` populates a working directory from
local FASTA files instead of downloading: files are linked (or copied with
`--copy`) into `genomes/`, `selection.tsv` and the SQLite manifest are written
in the same format `metadata`/`genome` produce, and the stage is recorded so
`status` reports the local chain (`ingest -> dereplicate -> phylo -> tree2tax`)
and a second `ingest` on an unchanged directory skips. Taxonomy and quality
columns come from `--selection selection.tsv` when given, otherwise from the
canonical `Family_genus_species_ACCESSION.fasta` filename. `--outgroup` sets
one genome aside (a name under `--genomes-dir`, or a path to a FASTA file
elsewhere) under `outgroup/` and writes `outgroup_accession.txt`, so `phylo`
roots on it exactly as after a download. This is the CLI counterpart of the
Nextflow harness `nextflow/tests/local_dataflow.nf`, which takes a genome
directory straight into `DEREPLICATE_SCATTER`.

### Representative selection

`repgenr dereplicate` and `repgenr run` (the manual and `run` CLI entry points,
not the Nextflow data-channel path) accept `--keeper quality|tool` (default
`quality`). After the chosen dereplicator clusters the genomes, the keeper step
re-picks each cluster's representative by `completeness - 5 x contamination`
using GTDB CheckM values already present in the manifest
(`src/repgenr/stages/derep_keeper.py`), which corrects the tendency of
connectivity-based tools to keep the most-sequenced (not the best-quality)
genome in a cluster. `--keeper tool` restores the adapter's own pick. Clusters
with no manifest quality data keep the adapter's choice either way. The GTDB
table carries CheckM values directly; `--source api` fetches them from each
genome's card (one request per selected genome). When the manifest has no
quality at all the stage warns, and `repgenr.yaml` records
`keeper_effective: tool` next to the requested `keeper` and the swap count.

### Limiting the selection

`repgenr metadata --limit N` caps the bacterial selection at N genomes. The
cap is not the first N rows of the GTDB table: candidates are grouped by
species and taken round-robin, the best-quality genome of every species first,
then each species' next best, until N. Within a species genomes rank by CheckM
completeness minus five times contamination (the same score the keeper uses),
unscored genomes last, then the GTDB species-representative flag, then
accession, so the result is deterministic. A heavily sequenced species
therefore cannot fill the cap on its own. With `-d rep` there is one genome
per species and the rule reduces to a quality ranking across species. On the
`--source api` path the per-genome quality cards are fetched for every
candidate before the cut, one request per genome, four at a time (about five
genomes a second; the 1540 Wolbachia genomes take some five minutes). The
API refuses sustained bursts now and then; refused cards are retried once
more, slowly, and only a genome refused twice is left unscored. The outgroup is chosen
afterwards from the parent taxon and never counts against the limit.

### Reusing an alignment across tree builders

`phylo` stamps the alignment it builds (`align/msa_source.json` or
`snp/msa_source.json`) with what produced it: the source and its settings, the
genome set, and the alignment's own digest. A later `phylo` run that changes
only the tree builder, the bootstrap or the thread count reuses that alignment
instead of aligning or SNP-calling again, and says so in the log. Anything the
alignment depends on, such as the aligner, the SNP typer, `--reference`,
`--mask` or the genome set, rebuilds it, as does `--force` or a change to the
alignment file itself.

The same split is available to the stateless step: `phylo-build --msa-only`
builds the alignment and writes `msa.fasta` without a tree, and `phylo-build
--msa <file>` builds a tree from an alignment an earlier call produced. The
Nextflow layer uses these to run the alignment and the tree as separate tasks.

### SNP typing and masking

The `repgenr snptype` command (and `phylo-build --msa-source snptype`) call a
SNP typer to produce a core-SNP alignment. Recombination masking (`--mask
gubbins`) runs on the typer's whole-genome alignment and replaces the
core-SNP alignment with Gubbins' filtered polymorphic sites. Typers that only
emit variable sites cannot be masked.

Gubbins expects isolates of one species. The masker estimates how much of the
alignment is variable, warns above 10%, and repeats the figure if Gubbins
fails: on a genus-level set its scan can die with a bus error (see
`docs/verification.md`).

Gubbins builds a tree in every iteration, by default with RAxML, and with more
than one thread it needs a multi-threaded RAxML build (`raxmlHPC-PTHREADS*`).
Some conda builds ship only the single-threaded binary, and Gubbins then exits
before its first iteration. When repgenr runs Gubbins natively and finds no
such build it switches to IQ-TREE (or to one thread when IQ-TREE is missing
too) and says so in the log. `--tool-arg gubbins_tree_builder=raxmlng`,
`--tool-arg gubbins_first_tree_builder=rapidnj` and
`--tool-arg gubbins_args="--min-snps 5"` pass the choice, the first-iteration
builder and any further `run_gubbins.py` arguments through.

The `simple` typer maps each genome independently, so `--threads` buys
concurrent genomes first and threads inside one genome's chain only when there
are more threads than genomes. Its per-genome intermediates are written
compressed and removed as soon as that genome's consensus has been read, so
scratch stays at a few hundred megabytes whatever the genome count. A genome
whose chain fails keeps its intermediates for inspection. Under a container
backend each genome's chain of tools runs in a single container, so a genome
costs one engine start rather than eight.

`--tool ska2` (split k-mer analysis) is reference-free: every genome is an
ordinary sample, so no assembly's private errors bias the SNP distances, and
the alphabetical-first-genome reference default does not apply. It emits a
variable-site alignment only, so it is not compatible with `--mask`. Tune with
`--tool-arg ksize=31` and `--tool-arg min_freq=0.9` (the fraction of samples a
split k-mer must occur in). On the 60 Wolbachia species representatives from
GTDB, a genus-level set, it returned 708 variable sites in about two minutes;
it is designed for clonal and outbreak sets, where the shared k-mer space is
far larger.

### Collapsing weak splits in tree2tax

`repgenr tree2tax` (and the `tree2tax-relations` step the Nextflow module
runs) can merge weakly supported or near-zero-length internal nodes into
their parents before naming them, so a split the data does not support does
not become its own FlexTaxD node. Both thresholds are off by default.
`--collapse-length L` merges a node whose branch is shorter than L, in the
tree's own units (mash distance for mashtree, substitutions per site for the
ML builders). `--collapse-support S` merges a node whose support is below the
fraction S; IQ-TREE and RAxML-NG write percentages and are normalised
automatically, FastTree writes fractions, and mashtree and the sourmash
builder write no supports, in which case the option warns and does nothing.
Either criterion suffices. The root, the outgroup/ingroup split under it and
the leaves never collapse. A collapsed node's branch length is added to its
children's. Provenance records both thresholds and the number of nodes
collapsed, so changing a threshold re-runs the stage. On the 171-leaf
Wolbachia mashtree tree, `--collapse-length 0.0005` removes 26 of 168 internal
nodes.

## Viral length filtering and over-represented species

The viral selection step keeps records whose genome length falls inside a
window around a center value. The default center, `--length-method
median_of_medians`, is a deliberate defense against over-sequenced species:
the median length is computed per species first, and the window center is the
median of those per-species medians, so a species with 900 outbreak records
carries exactly one vote. Switching to `--length-method mean` averages every
record and lets the most-sequenced species set the window -- use it only when
the selection is known to be balanced. `--length-range` overrides the window
entirely and `--length-all` disables the filter.
