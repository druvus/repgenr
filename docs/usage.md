# Usage

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

This page covers the command line first, then the Nextflow layer, running the
tools in containers, and troubleshooting. Every command and option is listed in
[cli-reference.md](cli-reference.md); the files each stage writes are described
in [output.md](output.md).

## Command line

### Bacteria

```bash
WD=./francisella

# Full GTDB table (release-pinned):
repgenr metadata -wd $WD -r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella
# Or query just the target taxon via the GTDB API (no full-table download):
repgenr metadata -wd $WD --source api -d rep -l genus -tg francisella

repgenr genome -wd $WD
repgenr dereplicate -wd $WD --tool skder -t 16
repgenr phylo -wd $WD --aligner progressivemauve --treebuilder iqtree
repgenr tree2tax -wd $WD --include-dereplicated
```

Or run the whole chain in one command (bacterial by default; `--viral` for the
NCBI Virus path), then check progress at any time:

```bash
repgenr run -wd $WD -d rep -l genus -tg francisella --tool skder --treebuilder iqtree
repgenr status -wd $WD     # which stages are done, and what to run next
```

`run` forwards the stage options it shares a name with; for the rest, run
the stage command by hand (`repgenr status` says which comes next). Two
worth knowing: `--with-snptype` adds the standalone `snptype` stage after
dereplication, so the SNP tables under `snp/` are produced even when the tree
is built another way (with `--msa-source snptype`, `phylo` still runs its own
typing pass into the same directory, and reuses it afterwards), and
`--genomes-dir` starts the chain from local genomes.

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

```bash
repgenr ingest -wd $WD --genomes-dir ./my_genomes --outgroup GCF_003574425.1
repgenr dereplicate -wd $WD --tool skder
repgenr phylo -wd $WD --treebuilder mashtree
repgenr tree2tax -wd $WD --include-dereplicated
```

`repgenr run --genomes-dir ./my_genomes ...` runs the same local chain in
one command; `--selection`, `--outgroup` and `--copy` pass through to
`ingest`, and the GTDB selection flags are not needed.

`--outgroup` names a genome under `--genomes-dir` (filename, stem or
accession) or a FASTA file anywhere; it is staged under `outgroup/` and kept
out of the ingroup.

### Starting from sequencing reads

`repgenr reads` selects whole-genome sequencing runs from ENA (which mirrors
SRA) and writes `reads.tsv`; `repgenr assemble` fetches and assembles them
into `genomes/` with the same `selection.tsv` and manifest the other entry
paths write, so the chain is `reads -> assemble -> dereplicate -> phylo ->
tree2tax`. Runs are chosen by taxon
(`--target-family`/`-tf`, `--target-genus`/`-tg` or `--target-species`/`-ts`,
resolved through the ENA taxonomy, synonyms included) or by accession: `--accession` takes a run (SRR/ERR/DRR), a sample
(SAMN.., SRS..) or a study (PRJNA.., SRP..) and repeats; `--accession-file`
lists them one per line. `--platform illumina|ont|pacbio` keeps one
platform, `--min-bases` drops small runs, `--max-bases` drops runs above a
size (an unenriched whole-host library, tens of Gb for a 1 Mb endosymbiont,
would assemble into a host-dominated genome), `--drop-selection` drops runs
by ENA library selection (default `MDA`, whole-genome amplification, which
assembles into chimeric and uneven contigs; repeat the flag for more values,
`--drop-selection none` keeps every run; the value is kept in `reads.tsv` as
`library_selection`), `--one-per-sample` (the default)
keeps the best run of each sample: a long-read run when it carries at least
100 Mb and a tenth of the sample's largest short-read run, else the largest
run (`--all-runs` keeps every run); `--max-runs` caps the selection to the
largest runs. At assembly, a paired run that ENA lists with a third, orphan
FASTQ file is given to skesa as the pair plus the orphan file, and to shovill
as the pair only. Each run is labelled with the family, genus and species of
its NCBI taxid, in the same filename tokens the GTDB path uses.

```bash
repgenr reads -wd $WD -ts "Francisella tularensis" --platform illumina --max-runs 20
repgenr reads -wd $WD --accession PRJNA954307 --accession SRR28800588
repgenr assemble -wd $WD --assembler auto -t 16 --jobs 2
# or the whole chain:
repgenr run -wd $WD --reads -ts "Francisella tularensis" --platform illumina --max-runs 20 \
    --assembler skesa --treebuilder mashtree
```

`assemble` downloads each run's FASTQ files from the locations ENA lists
(over HTTPS, verified against ENA's checksums), assembles them, keeps the
contigs of at least `--min-contig-length` bases (500) renamed
`<run>_contig<n>`, and names the genome `Family_genus_species_RUN.fasta` from
the tokens the reads stage resolved. `--assembler auto` (the default) picks
by platform and layout: `skesa` for Illumina, `shovill` (SPAdes) as the
alternative for paired Illumina, `flye` for Oxford Nanopore and PacBio;
`--tool-arg` passes tuning such as `mode=nano-raw` to Flye. `--jobs` runs
that many assemblies at once with `--threads` split across them; memory,
not CPU, is the limit, so the default is 2, or 1 as soon as a long-read run
is pending. `--memory-gb` is the RAM cap passed to SKESA and shovill (shovill
accepts no less than 8). Each finished run leaves a marker under
`assemblies/<run>/`, so an interrupted stage resumes without refetching;
reads are deleted after a successful assembly unless `--keep-reads`, and the
assembler's scratch unless `--keep-files`. A run without an ENA FASTQ
mirror, one whose download fails its checksum, one no assembler accepts, or
one whose assembly fails is written to `excused_runs.tsv` with the reason and
the rest proceed; the completeness guard of later stages excuses those runs.
`--outgroup FASTA` sets a genome aside for rooting, as `ingest --outgroup`
does. Per-assembly metrics (contigs, total length, N50, coverage from the
sequenced bases) are in `assembly_stats.tsv`.

`assemble --append` adds the assemblies to a working directory that already
holds a selection from `metadata` and `genome` (or `ingest`): the existing
rows, genome files and outgroup stay, a run assembled earlier is replaced,
and the manifest gains the new genomes with source `sra`, so a mixed GTDB
and reads-derived set dereplicates together. Because `metadata` and `ingest`
replace the selection and the genome stage prunes what the manifest no
longer lists, both refuse to re-run while appended genomes are present;
`--drop-foreign` discards them deliberately, and `assemble --append` can put
them back afterwards.

Two optional checks run on the assemblies. With a CheckM2 database
(`--checkm2-db`, or the `CHECKM2DB` variable CheckM2 itself reads; obtain it
with `checkm2 database --download`), every assembly is scored, the
completeness and contamination reach `selection.tsv` and the manifest (so
`--keeper quality` works as it does for GTDB genomes), and an assembly below
`--min-completeness` (50) or above `--max-contamination` (10) is excused with
`qc_failed`. With a GTDB sourmash sketch (`--gtdb-sketch` and
`--gtdb-lineages`, or `REPGENR_GTDB_SKETCH` and `REPGENR_GTDB_LINEAGES`; the
`gtdb-rs226-reps.k31-sc10k.sig.zip` sketch and its `lineages.csv` from
`https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/` serve), each
assembly is classified by `sourmash gather` (`--classifier auto` runs it when
a sketch is configured; `none` never). When the GTDB genus agrees with the
submitted organism, the GTDB family, genus and species name the genome file,
so a reads-derived genome groups with GTDB-downloaded ones; otherwise the
submitted name stays and `assembly_stats.tsv` flags the genome
`classifier_disagrees`. Both lineages are kept in that table, and the sketch
release is recorded in provenance next to the metadata release. Without a
database the checks are skipped and the log says so. `run --reads` forwards
`--accession-file`, `--platform`, `--max-runs`, `--assembler`, `--threads`
and `--outgroup`; the rest is available on the stage commands.

The same work is available as three stateless steps, which the Nextflow reads
mode runs as separate tasks and which also serve a scheduler of your own:

```bash
repgenr reads -wd $WD -tg mycoplasmopsis --platform illumina --max-runs 20
# one task per run: fetch, assemble, filter contigs
repgenr assemble-run --reads-tsv $WD/reads.tsv --run SRR25474756 -o asm/SRR25474756 \
    --assembler auto -t 8 --memory-gb 16 --min-contig-length 500
# one batch: CheckM2 and/or the classifier over every finished run directory
repgenr genome-qc --assemblies asm -o qc --checkm2-db $CHECKM2DB \
    --gtdb-sketch gtdb-rs226-reps.k31-sc10k.sig.zip --gtdb-lineages lineages.csv \
    --classifier auto -t 16
# the genome contract: genomes/, selection.tsv, assembly_stats.tsv, excused_runs.tsv
repgenr reads-gather --reads-tsv $WD/reads.tsv --assemblies asm --qc qc -o out \
    --min-completeness 50 --max-contamination 10
```

`assemble-run` writes `contigs.fasta` and the `assembly.ok` marker into its
`--out` directory, or `excused_runs.tsv` when the run has no FASTQ mirror, an
unsupported platform, or fails to download or assemble (`--keep-reads`,
`--keep-files` and `--tool-arg` as on `assemble`). `genome-qc` reads a
directory of such run directories (`--assemblies`) and writes `quality.tsv`
and `classification.tsv` keyed by run accession; it needs at least one
database. `reads-gather` applies the quality gate and the naming policy above
from the optional `--qc` directory and writes an empty
`outgroup_accession.txt`, since the reads chain has no outgroup at this
point. None of the steps touches a working directory or the manifest.

### Viruses

The viral path selects from NCBI Virus by default (via the `datasets` CLI);
`vmetadata --source bvbrc` uses the legacy BV-BRC FTP path instead.

```bash
WD=./hav
repgenr vmetadata -wd $WD --target hepatovirus            # NCBI Virus (default)
repgenr vgenome   -wd $WD --target-genus Hepatovirus      # add --group-segments for segmented viruses
repgenr dereplicate -wd $WD --tool skder --virus
repgenr phylo -wd $WD --treebuilder mashtree
repgenr tree2tax -wd $WD --include-dereplicated
# or: repgenr run -wd $WD --viral --target hepatovirus -tg Hepatovirus --treebuilder mashtree
```

`vgenome` picks an outgroup by itself: a record of a sister species with
enough genomes, chosen by distance (mashtree by default). `--outgroup-accession`
pins it to a downloaded record instead, an accession on the NCBI Virus path
or a record id on BV-BRC, and `run --viral --outgroup-accession` forwards
the same choice. With `--group-segments` the search still runs, with the
kept records' length span as its window, and the outgroup is one record of
the sister species rather than a grouped isolate. `--no-outgroup` leaves the
tree unrooted. A grouped isolate's genome carries a synthetic `iso-` token as
its accession; `segments.tsv` records the member accessions behind it, and
`tree2tax` lists them under the isolate's leaf in `genomes_map.tsv`.

### Viral length filtering and over-represented species

The viral selection step keeps records whose genome length falls inside a
window around a center value. The default center, `--length-method
median_of_medians`, is a deliberate defense against over-sequenced species:
the median length is computed per species first, and the window center is the
median of those per-species medians, so a species with 900 outbreak records
carries exactly one vote. Switching to `--length-method mean` averages every
record and lets the most-sequenced species set the window -- use it only when
the selection is known to be balanced. `--length-range` overrides the window
entirely and `--length-all` disables the filter.

### Alignment-free and SNP-based phylogenies

Two alternatives to the whole-genome alignment in the bacterial example:

```bash
# Scalable dereplication then an alignment-free tree
repgenr dereplicate -wd $WD --tool skder
repgenr phylo -wd $WD --treebuilder mashtree

# SNP-based phylogeny (core-SNP alignment as the MSA source)
repgenr snptype -wd $WD --tool simple
repgenr phylo -wd $WD --msa-source snptype --treebuilder iqtree --mask gubbins
```

### Resume and `--force`

Each stage records its parameters, the digests of its inputs, and the container
identity in `repgenr.yaml`; re-running a stage is a safe no-op only when all
three are unchanged (it logs that it skipped). Re-running an upstream stage
(e.g. `dereplicate --force`) changes a downstream stage's input digests, so the
downstream stage re-runs automatically the next time it is invoked. Change a
parameter, switch `--container`, or pass `--force` to re-run explicitly. A
stage that crashed mid-run has no completion stamp and so always re-runs. `repgenr
doctor -wd <wd>` verifies a workdir's outputs against its records (missing or
corrupt genomes, manifest drift, truncated deliverables, interrupted stages)
and exits non-zero on failures.

Two limitations, both covered by `--force`: input directories are digested from
file metadata (name, size, mtime), so an in-place edit that preserves size and
mtime is not detected; and upgrading a natively installed tool binary does not
invalidate previous results (switching the container backend or platform does).
Workdirs created by older RepGenR versions re-run each stage once (the
fingerprint format changed).

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

`--reduce species|genus` collapses the ANI representatives to one per taxon
after dereplication, choosing the keeper by quality when scores are known and
by cluster size otherwise; `--target-reps N` searches the secondary ANI to
land near N representatives. Both exist on `dereplicate` and, since the
merge step is where the final set is decided, on `dereplicate-merge`, which
the Nextflow layer drives through `--derep_reduce` and `--derep_target_reps`.
The merge step takes the taxonomy from `selection.tsv` when one reaches it
and from the canonical genome filenames otherwise.

### Inspecting a dereplication

Four commands read a dereplicated working directory without rerunning the
dereplicator. `repgenr cluster-summary` regenerates
`derep/cluster_summary.tsv`, one row per representative (see `output.md`).
`repgenr derep-unpack` lays the clusters out as one directory per
representative with its members inside (`--no-representant` leaves the
representative out). `repgenr glance` runs dRep's comparison over the
representatives and writes its plots (dRep only). `repgenr derep-stock
--action pack --name <run>` stores the current clusters, statuses and
representatives under `derep/stock/<run>`; `--action unpack` restores a
stored run, refreshes the manifest and re-stamps the `dereplicate` record so
the next `dereplicate` recomputes; `--action list` and `--action delete`
manage the store.

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
SNP typer to produce a core-SNP alignment. Four typers are built in: `simple`
(minimap2, samtools and bcftools, the default), `snippy`, `parsnp` and `ska2`;
the first three map every genome to one reference and also write the
whole-genome alignment a masker needs. Recombination masking (`--mask
gubbins`) runs on the typer's whole-genome alignment and replaces the
core-SNP alignment with Gubbins' filtered polymorphic sites. Typers that only
emit variable sites cannot be masked.

Gubbins expects isolates of one species. The masker estimates how much of the
alignment is variable, warns above 10%, and repeats the figure if Gubbins
fails: on a genus-level set its scan can die with a bus error (see
`verification.md`).

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

## Nextflow

```bash
nextflow run nextflow/main.nf \
    --mode bacterial \
    --metadata_args '-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella' \
    --outdir results \
    -profile standard
```

Run `nextflow run nextflow/main.nf --help` for the parameter summary.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Published results and execution reports. |
| `--mode` | `bacterial` | `bacterial` (GTDB), `viral` (NCBI Virus; BV-BRC via `--vmetadata_args "--source bvbrc"`) or `reads` (ENA/SRA sequencing runs, assembled per run). |
| `--metadata_args` | see config | Arguments for the bacterial metadata stage. |
| `--vmetadata_args` / `--vgenome_args` | see config | Viral metadata / genome selection arguments. |
| `--reads_args` | see config | Arguments for the reads stage (ENA/SRA run selection by taxon or accession). |
| `--assembler` | `auto` | Assembler for every run in reads mode (`auto` picks by platform). |
| `--assemble_args` | (empty) | Extra `assemble-run` flags in reads mode (`--min-contig-length`, `--tool-arg`). |
| `--checkm2_db` | `null` | CheckM2 database; switches on quality scoring and the completeness/contamination gate in reads mode. |
| `--gtdb_sketch` / `--gtdb_lineages` | `null` | GTDB sourmash sketch and its lineages CSV; switch on classification and GTDB naming in reads mode. |
| `--derep_tool` | `skder` | Dereplicator for the scatter-gather step. |
| `--derep_process_size` | `null` | Genomes per dereplication chunk (single chunk if unset). |
| `--derep_primary_ani` / `--derep_secondary_ani` / `--derep_aligned_fraction` | `0.90` / `0.99` / `0.50` | ANI / aligned-fraction thresholds. |
| `--derep_keeper` | `quality` | Representative choice at the merge step: `quality` (CheckM scores from `selection.tsv`) or `tool`. |
| `--derep_reduce` | `none` | Collapse the merged representatives to one per `species` or `genus` (`dereplicate-merge --reduce`). |
| `--derep_target_reps` | `0` | Search the merge pass's secondary ANI to land near this many representatives (`dereplicate-merge --target-reps`). |
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

#### Passing parameters

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

### Profiles

The `standard` profile (the default, local executor) caps every process's
CPU and memory request to what the machine has; `base.config` asks for up to
32 CPUs and 128 GB for the heavy processes, which the local executor would
otherwise refuse. `slurm` leaves the requests as they are and sets only the
executor: the queue or account is site-specific, so pass it in a site config
(`-c site.config` with `process.queue = '...'`). The same holds for a cloud
executor such as AWS Batch: the pipeline ships no cloud profile, because the
region, the job queue and an image that provides `repgenr` and the tools are
values only the site knows; a site config that sets them is all that is
needed.

Combine an executor profile with an optional container profile, e.g.
`-profile slurm,singularity`.

- **Executors**: `standard` (local), `slurm`.
- **Containers**: `docker`, `singularity`, `wave`. These set RepGenR's own
  adapter-level container backend (`--container ...`). Under `docker` or
  `singularity` every adapter with a pinned image runs in it (all but the
  `simple` SNP typer, parsnp, skder and SibeliaZ, which run on the host and
  say so in the log); `wave` mints an image from each adapter's conda
  specification instead, which also covers those four. RepGenR itself must
  be available to the Nextflow process.
- **`test`**: minimal resources and a small target for a quick smoke run.

`PHYLO` (and the split `PHYLO_MSA`/`PHYLO_TREE`) publish the alignment they
built (`phylo/align/` for an aligner, `phylo/snp/` for a SNP typer, with the
reuse stamp) and the tree builder's own files under `phylo/tree/`, beside the
tree.

### Scaling

For large dereplication inputs (10k+ genomes), set a chunk size so the
dereplication scatters across tasks (one per chunk):

```bash
nextflow run nextflow/main.nf --outdir results \
    --derep_tool sourmash --derep_process_size 2000 -profile slurm -c site.config
```

Resource labels (`process_low/medium/high`, and `process_assembly` for the
per-run assembly tasks of reads mode, whose memory follows the read platform)
scale memory and time with the retry attempt, so a task killed for memory or
time is resubmitted with more headroom. Tune the label values per environment
in `nextflow/conf/base.config`.

#### Scatter-gather dereplication

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

#### Configuring processes

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

#### Pipeline structure

`nextflow/main.nf` dispatches by `--mode` to one of three data-channel
subworkflows that share the dereplication, phylo and tree2tax modules:

```
bacterial: ACQUIRE       (metadata -> genome)          -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
viral:     VACQUIRE      (vmetadata -> vgenome)        -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
reads:     ACQUIRE_READS (reads -> assemble per run
                          -> [genome-qc] -> gather)    -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
```

`metadata` emits a portable `selection.tsv`; `genome-fetch` downloads the genomes
and emits them as a channel feeding the scatter-gather dereplication; `phylo` and
`tree2tax` run in task-local working directories and emit `tree.nwk`,
`tree2tax.tsv` and `genomes_map.tsv` to `--outdir`. Add `-stub` to any run for a
quick wiring check without external tools.

In reads mode `READS_SELECT` runs `repgenr reads` and emits `reads.tsv`; every
row becomes one `READS_ASSEMBLE` task (`assemble-run`, label
`process_assembly`: 8 CPUs, 16 GB for a short-read run and 32 GB for a
long-read one, 6 h, scaled by attempt), so a hundred runs assemble in
parallel on a cluster. When `--checkm2_db` or `--gtdb_sketch` is set, one
`GENOME_QC` task (`genome-qc`, `process_medium`) scores and classifies the
batch; `READS_GATHER` (`reads-gather`) applies the gate and the naming
policy and emits the genome FASTAs, `selection.tsv` and an empty outgroup
accession file, the same tuple `ACQUIRE` emits. `reads.tsv`,
`selection.tsv`, `assembly_stats.tsv`, `excused_runs.tsv` and the QC tables
are published under `reads/`.

```bash
nextflow run nextflow/main.nf --mode reads --outdir results \
    --reads_args "-tg mycoplasmopsis --platform illumina --max-runs 20" \
    --assembler auto --checkm2_db /db/checkm2 \
    --gtdb_sketch /db/gtdb-rs226-reps.k31-sc10k.sig.zip --gtdb_lineages /db/lineages.csv \
    -profile slurm,singularity
```

## Running tools in containers

RepGenR can run each external tool inside a pinned container instead of relying
on tools installed on `PATH`. This unblocks Linux-only tools on any host (e.g.
macOS) and makes tool versions reproducible. RepGenR itself runs on the host;
only the tool subprocess is containerized, so it works in the plain CLI and
inside Nextflow.

### Quick start

```bash
# Run every tool in a container via Docker:
repgenr --container docker dereplicate -wd $WD --tool drep

# Singularity/Apptainer (HPC), with .sif images on a large disk:
repgenr --container singularity --container-cache /Volumes/LaCie/repgenr_sif \
        phylo -wd $WD --aligner progressivemauve --treebuilder iqtree
```

`--container` is a top-level option (place it before the subcommand). Default is
`none` (native execution, unchanged).

### Options (top-level, or env var)

| Option | Env var | Meaning |
|--------|---------|---------|
| `--container {none,docker,singularity}` | `REPGENR_CONTAINER` | execution backend |
| `--container-engine <bin>` | `REPGENR_CONTAINER_ENGINE` | engine override (apptainer, podman) |
| `--container-cache <dir>` | `REPGENR_CONTAINER_CACHE` | Singularity `.sif` / Wave cache (can be external) |
| `--platform <plat>` | `REPGENR_CONTAINER_PLATFORM` | e.g. `linux/amd64` to emulate BioContainers on Apple Silicon |
| `--wave / --no-wave` | `REPGENR_WAVE` | resolve multi-tool/arm64 images via the Seqera Wave CLI |

### Image sources

Each adapter declares its container metadata in `ToolCapabilities`:
- `container` — a pinned image URI, used as-is when Wave is off. The
  single-package adapters pin a BioContainer (galah, sourmash, dRep, snippy,
  ska2, Gubbins, IQ-TREE, FastTree, RAxML-NG, mashtree), progressiveMauve
  pins a BioContainer with a compatible boost, and cactus pins its project
  image. Four adapters have no pin and run in an image only under `--wave`
  (on the host otherwise, which the log states for each): the `simple` SNP
  typer (minimap2, samtools, bcftools) and parsnp (parsnp, harvesttools)
  span several packages, and the skder and SibeliaZ BioContainers are
  BusyBox-based, where the GNU-only calls in their shell wrappers (`sort
  --parallel` in skder's greedy mode, `mktemp --suffix` in SibeliaZ) fail
  and the run silently yields nothing.
- `conda` — a conda spec (e.g. `bioconda::skder`). With `--wave`, RepGenR mints
  an image for it via the Wave CLI (arm64-native, and the only route for the
  multi-package adapters) and uses it instead of the pin; the pin stays the
  default without Wave.

The pinned tags are listed in each adapter's `capabilities`
(`repgenr list-tools` names the adapters). BioContainers are `linux/amd64`;
on Apple Silicon pass `--platform linux/amd64` (Docker Desktop with Rosetta)
or use `--wave`.

### Storage location

- **Singularity/Apptainer:** `--container-cache` sets `APPTAINER_CACHEDIR` /
  `SINGULARITY_CACHEDIR` (+ `*_TMPDIR`); `docker://` images are pulled once to
  `<cache>/<name>.sif` and reused. Put this on a large/external disk.
- **Docker:** image storage is managed by the Docker daemon (Docker Desktop's
  disk-image location) and is set there, not per-run.
- Large run-time data (Cactus jobstore, scratch, downloads) lives under
  `--workdir` / `TMPDIR`.

### Notes

- Docker runs as the host UID/GID so outputs are owned by you; the workdir and
  `TMPDIR` are bind-mounted at identical paths.
- Symlinked inputs (a `genomes/` directory staged by `repgenr ingest`) are
  followed: the directory each link points to is bound as well, so the tool
  sees the same paths inside the container.
- On Apple Silicon most BioContainers are `linux/amd64` (run via Docker
  emulation, or use `--wave` for arm64-native images).
- The macOS SibeliaZ BSD-wrapper workaround is skipped automatically when running
  in a (Linux) container.
- dRep's CheckM needs its reference DB at run time — mount it via
  `CHECKM_DATA_PATH`, or run dRep with `--ignoreGenomeQuality`.
- The bioconda `mauve` (progressiveMauve) build is broken upstream (boost ABI
  `undefined symbol`); pin a known-good image or run that tool natively on Linux.

### Container profiles in Nextflow

`nextflow run nextflow/main.nf -profile docker` (or `singularity`, `wave`) sets
`params.repgenr_opts` so every stage calls `repgenr --container …`. On HPC,
Singularity is the clean target (RepGenR dispatches each tool to Singularity);
stacking with Nextflow's own Docker engine implies docker-in-docker.

## Troubleshooting

- **`MissingBinaryError` / a tool is not found.** The Python package does not
  install the bioinformatics tools. Use the conda environment
  (`mamba env create -f environment.yml`) or put the tool on `PATH`. Run
  `repgenr list-tools` to see the adapters (`list-tools --check` also runs
  every adapter's preflight and reports which binaries are missing or too
  old) and `--container docker` (or
  `singularity`) to run tools in pinned images instead.
- **Apple Silicon / arm64.** BioContainers are amd64; pass
  `--platform linux/amd64` (and enable Rosetta) so emulated images run.
- **A stage failed; where are the details?** Errors print a concise message; the
  full traceback is in `<workdir>/repgenr.log`. Re-run with `--verbose` to see it
  on the console. `repgenr status -wd <WD>` shows what completed and what is next.
- **GTDB download fails.** Check `--release` (e.g. `232.0`) and `--gtdb-version`
  (`bac120`/`ar53`); transient HTTP errors are retried automatically. The
  `--source api` mode fetches only the target taxon (no full-table download).
- **NCBI Entrez throttling (viral BV-BRC path).** Set `NCBI_API_KEY` (and
  optionally `NCBI_EMAIL`) to raise the request-rate limit.
- **A tool hangs.** Set `REPGENR_SUBPROCESS_TIMEOUT=<seconds>` to cap every
  external tool; on expiry the process group is killed with a clear error.
- **Exit codes.** A script can tell the failure classes apart without
  reading the log:

  | Code | Meaning |
  |------|---------|
  | 0 | Success. |
  | 1 | An unexpected error (traceback in the run log), or `doctor` found failures. |
  | 2 | Invalid or missing user input (also Typer's own usage errors). |
  | 3 | The working directory is missing files or is in a bad state. |
  | 4 | A required external tool is absent or below its version floor. |
  | 5 | A requested tool adapter could not be found or loaded. |
  | 6 | An external tool failed. Under `REPGENR_PROPAGATE_TOOL_EXIT=1` (set by the Nextflow modules) the tool's own status is forwarded instead, a signal kill as 128 plus the signal number. |
