# Changelog

All notable changes to RepGenR are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- `repgenr --help` groups the commands into panels (pipeline, entry points,
  core stages, inspection, environment, Nextflow steps) in pipeline order, and
  its epilog lists the four stage chains (bacterial, viral, local genomes,
  sequencing reads).
- `list-tools` prints each adapter's declared genome limit, for example
  `iqtree (up to 500 genomes)`; a tool without a declared limit shows its name
  only.
- An eighth tool family, polishers (`repgenr.polishers`), corrects long-read
  assemblies with the run's reads inside `assemble` and `assemble-run`:
  medaka for ONT (basecaller model from Dorado read headers, `--tool-arg
  model=...`, or ONT's bacterial methylation model assumed for SRA-mirrored
  reads whose headers name none; recorded in the marker) and racon with minimap2 overlaps for PacBio CLR
  (`--polish-rounds`). `--polisher auto|medaka|racon|none`; HiFi and
  Illumina are never polished. `assembly_stats.tsv` gains a `polisher`
  column and a failed polish excuses the run with `polish_failed`.
- `reads.tsv` records ENA's `library_selection`, and `reads --drop-selection`
  (default `MDA`, repeatable, `none` to keep all) drops amplified libraries,
  which assemble into chimeric and uneven contigs (68 of the 551 Wolbachia
  WGS runs are MDA).
- `reads --max-bases N` drops runs above a size, a guard against unenriched
  whole-host libraries (36 of the 551 Wolbachia runs exceed 5 Gb) that would
  assemble into a host-dominated genome.
- Nextflow `--mode reads`: the reads chain as data-channel tasks. `READS_SELECT`
  runs the reads stage, `READS_ASSEMBLE` assembles one run per task (new
  `process_assembly` label: 8 CPUs, 16 GB for short reads and 32 GB for long
  reads via the run's platform), `GENOME_QC` scores and classifies the batch
  when `--checkm2_db` or `--gtdb_sketch` is set, and `READS_GATHER` writes the
  genome contract that feeds the shared dereplication, phylo and tree2tax
  modules. Parameters `reads_args`, `assembler`, `assemble_args`,
  `checkm2_db`, `gtdb_sketch`, `gtdb_lineages`; tables published under
  `reads/`. Behind the modules, three stateless CLI steps: `assemble-run`,
  `genome-qc` and `reads-gather`. The `assemble` stage now keys its QC and
  classification inputs by run accession rather than by the final genome
  name.

- `assemble --append` adds reads-derived genomes to a working directory that
  already holds a GTDB or local selection, keeping its rows, files and
  outgroup; `metadata --drop-foreign` and `ingest --drop-foreign` are the
  explicit way to discard them, since both stages otherwise refuse to
  overwrite a selection holding appended genomes.

- Assembly quality and organism verification for the reads chain. With a
  CheckM2 database (`--checkm2-db` or `CHECKM2DB`) every assembly is scored,
  the values reach `selection.tsv` and the manifest for the quality keeper,
  and assemblies outside `--min-completeness`/`--max-contamination` are
  excused. A seventh tool family, classifiers (`repgenr.classifiers`), with
  `sourmash` gather against a GTDB sketch (`--gtdb-sketch`,
  `--gtdb-lineages` or the matching environment variables): when the GTDB
  genus agrees with the submitted organism the GTDB tokens name the genome,
  otherwise the submitted name stays and the genome is flagged
  `classifier_disagrees` in `assembly_stats.tsv`. The sketch release is
  recorded in provenance.

- `repgenr assemble`, the second stage of the reads chain: fetches each
  selected run's FASTQ files from ENA over HTTPS with checksum verification,
  assembles them with the adapter that accepts the platform (`--assembler
  auto|skesa|shovill|flye`), filters contigs, names the genome from the
  reads stage's taxonomy and writes `genomes/`, `selection.tsv`, the
  manifest, `assembly_stats.tsv` and `excused_runs.tsv`. Finished runs carry
  a marker and are skipped on a re-run; `--jobs` bounds concurrency. `run
  --reads` runs the whole chain from a taxon or an accession file.

- `repgenr reads`, the first stage of a reads chain (`reads -> assemble ->
  dereplicate -> phylo -> tree2tax`): selects whole-genome sequencing runs
  from ENA (which mirrors SRA) by taxon or by run, sample or study accession,
  filters by platform and size, keeps the best run per sample, labels each
  run with its NCBI family, genus and species, and writes `reads.tsv`.
  `status` recognises the chain. The assemble stage follows.
- A sixth tool family, assemblers (`repgenr.assemblers` entry points), with
  `skesa` and `shovill` for Illumina and `flye` for Oxford Nanopore and
  PacBio reads, each pinning a BioContainers image, plus the shared contig
  filter and summary (contig count, total length, N50, largest contig) the
  assemble stage will use. `list-tools` lists the family.

- `vgenome --outgroup-accession` (and `run --viral --outgroup-accession`)
  pins the viral outgroup to a downloaded record on either back-end; the
  BV-BRC path previously asked the user to "specify one manually" with no
  way to do so. Segment-grouped runs (`--group-segments`) now search for an
  outgroup too, using the kept records' length span as the window; they were
  left unrooted before.

- `repgenr run` starts from local genomes with `--genomes-dir` (the `ingest`
  chain, with `--selection`, `--outgroup` and `--copy` passed through), runs
  the standalone SNP typing stage with `--with-snptype`, and forwards the
  stage options it did not before: `--metadata-path`, `--nodownload`,
  `--keep-files`, `--pre-primary-ani`, `--pre-secondary-ani`,
  `--node-basename`, `--root-name` and `--remove-outgroup`.

- `dereplicate-merge --reduce species|genus` and `--target-reps N`, the two
  representative-selection features that only the workdir `dereplicate`
  command had, so the Nextflow layer can use them through `--derep_reduce`
  and `--derep_target_reps`. The merge step takes the taxonomy from
  `selection.tsv` when one reaches it and from the genome filenames
  otherwise.

- `repgenr list-tools --check` runs every adapter's preflight and prints one
  line per tool: `ok` with the resolved versions, or `missing` with the
  binaries that are absent or below their version floor. An environment can
  be verified before a run, without a working directory.

- `derep/cluster_summary.tsv`: one row per representative with the member
  count, the species the cluster spans and the keeper's CheckM quality against
  its members, including which genome scores best. `dereplicate` writes it
  alongside `clusters.tsv` (the Nextflow chunk and merge steps too), and
  `repgenr cluster-summary -wd <workdir>` regenerates it for an existing
  working directory without rerunning the dereplicator.

### Changed
- A failed external tool prints one console line that names the tool and its
  exit status and points to `repgenr.log`. The command line and the output
  tail are written to the run log and shown on the console under `--verbose`;
  a data-channel step, which has no run log, still prints the tail.
- Resume now checks a stage's main outputs before skipping it. A stage whose
  parameters, inputs and environment are unchanged but whose deliverable was
  deleted (for example `derep/clusters.tsv` or `genomes_map.tsv`) logs
  `Stage 'X': deliverable <path> missing; re-running.` and runs again instead
  of being skipped, so `--force` is no longer needed for this case. `doctor`
  reports the same paths as warnings.
- One help string per shared flag: duplicated option texts are now shared
  constants, so `run` and the single commands read the same. `--mask` lists
  the registered maskers on snptype, phylo, run and phylo-build, and the three
  alignment commands say it needs `--msa-source snptype`. The known option-name
  differences (`--tool`, `--platform`, `--metadata-source`, `--viral-source`,
  and `--outgroup-accession` on the step commands) are described in the
  command help and the reference. No option was renamed.
- The sourmash classifier runs its per-genome gathers side by side within the
  thread budget and resolves all of them with one `tax genome` call (a gather
  against a GTDB-sized sketch is single-threaded and takes tens of seconds).
  The `assemble` record now stores the CheckM2 database, GTDB sketch and
  lineages as resolved from the flags or the environment; a re-run of
  `assemble` for quality or classification alone no longer requires the
  assembler binaries; the classifier disagreement warning names the
  compared genera.
- `assemble --jobs` defaults to 2, or to 1 as soon as a long-read run is
  pending, since memory rather than CPU bounds concurrent assemblies. The
  Flye path is verified through its pinned image on simulated Nanopore-like
  reads.
- CI lints and format-checks `scripts/` and `benchmarks/` as well as the
  package and the tests; `scripts/fasta_simulate_sequences.py`, which failed
  both, is rewritten with the same command line.
- The CLI exits with a distinct status per failure class instead of 1 for
  everything: 2 for invalid input, 3 for a bad working directory, 4 for a
  missing or outdated tool, 5 for an unknown adapter, 6 for a failed tool
  (the tool's own status under `REPGENR_PROPAGATE_TOOL_EXIT=1`, as before);
  1 stays for unexpected errors and for `doctor` finding failures. A
  Gubbins failure keeps its tool error and return code, with the divergence
  diagnosis appended, instead of being re-raised as a workdir error.
- Every version-checkable binary now carries a minimum-version floor at
  preflight, and `environment.yml` pins the same floors: dRep 3.0, galah
  0.4, cactus 2.5, parsnp 2.0, harvesttools 1.3, snippy 4.6, Gubbins 3.0,
  FastTree 2.1 and IQ-TREE 2.0 join the floors already declared. dRep's
  version is read from its help banner, as it answers no `--version` flag
  (the preflight reported it as unknown). progressiveMauve, SibeliaZ and
  hal2maf print no parseable version and stay lenient; a test guards the
  set.
- The single-package adapters pin a BioContainers image (galah 0.4.2,
  sourmash 4.9.4, dRep 3.7.1, snippy 4.6.0, ska2 0.5.1, Gubbins 3.4.3,
  IQ-TREE 3.1.3, FastTree 2.2.0, RAxML-NG 2.0.3, mashtree 1.4.6), each
  verified by running its stage through the image, so `--container
  docker|singularity` and the Nextflow container profiles run them in an
  image without Wave. Before, only progressiveMauve and cactus had an image
  and the other adapters fell back to the host with a warning. The
  multi-package `simple` SNP typer and parsnp still need `--wave`, and so
  do skder and SibeliaZ: their BioContainer images are BusyBox-based and
  the GNU-only calls in their shell wrappers (`sort --parallel`, `mktemp
  --suffix`) fail there silently.
- `--wave` now takes precedence over a pinned image when the adapter has a
  conda specification: the flag asks for an image minted from the spec (a
  native-architecture build), and the pin is the default without it.
- The Nextflow phylogeny processes publish the alignment they built
  (`phylo/align/` or `phylo/snp/`, with the reuse stamp) and the tree
  builder's own files under `phylo/tree/`; before, only `tree.nwk` left the
  task directory, so the SNP tables and the tree builder's logs were lost.
- The `cloud` Nextflow profile is removed and the `slurm` profile no longer
  names a queue: an executor's queue, account or region is site-specific and
  belongs in a site config (`-c`). The deprecated `nextflow.enable.dsl = 2`
  line is gone from the pipeline and its test harnesses.
- Documentation repairs from the gap analysis: the README counts five tool
  families and shows the chain `run` executes; `docs/output.md` lists
  `selection.tsv`, `outgroup_accession.txt`, `missing_accessions.txt`, the
  viral tables, the SNP distance matrix, the scratch directory, the aligner
  output under `align/` and the alignment stamp; `docs/verification.md`
  gains a per-adapter status table; `docs/developing.md` describes the layout
  as it is and the adapter checklist in full; the links broken by the
  consolidation are repaired; CITATIONS.md adds CheckM, Prodigal, Mash,
  FastANI and the HAL toolkit; `docs/usage.md` covers the inspection commands
  and names every SNP typer.
- `--mask` is refused before any SNP calling when the typer cannot feed a
  masker (SNP typers declare `produces_full_alignment`; ska2 does not), and
  the masker's binaries are checked before the typer runs. `repgenr run`
  preflights the masker with the other tools, so a missing Gubbins fails in
  the first second instead of after SNP typing.
- Adapters declare the standard parameters they do not pass to their tool
  (`ToolCapabilities.ignored_params`), and every stage warns by name when a
  user sets one of them to a non-default value. Before, `--primary-ani` was
  dropped without a message by skder and sourmash, `--aligned-fraction` by
  dRep and sourmash, and `--bootstrap` by mashtree and the sourmash tree
  builder.
- galah now receives `--primary-ani` as its pre-clustering ANI and
  `--aligned-fraction` as its minimum aligned fraction; both were ignored.
- The thread budget now reaches `sourmash compare` (dereplicator dense path
  and tree builder, `--processes`), FastTreeMP (`OMP_NUM_THREADS`) and
  cactus-pangenome (`--maxCores`); each previously ran on whatever the tool
  chose by itself.
- The documentation is consolidated from fifteen pages to six plus two audit
  records. `docs/usage.md` is the one how-to (command line, Nextflow,
  containers, troubleshooting; it absorbs `containers.md` and the README's
  longer sections), `docs/developing.md` merges `architecture.md` and
  `adding-tools.md`, and the scaling audit and its three SWOT pages are one
  record under `docs/audit/`. `docs/README.md` indexes the pages. The
  executed plans under `docs/superpowers/` and the unreferenced legacy figures
  under `docs/images/` are removed; git history keeps them.
- A stage's first run writes a provisional record to `repgenr.yaml` before its
  body starts, so a stage that fails or is killed shows as `[interrupted]` in
  `status` and as a failure in `doctor`, where it previously looked as if it
  had never started. A successful run replaces the record. A failure in
  parameter validation, a query-only invocation, and a stage that refuses its
  input (exit 2 or 3) without changing its main outputs write nothing.
- `phylo` and `phylo-build` compare the leaves of the built tree with the
  input genomes and exit 3, naming the missing or unexpected leaves and the
  tree builder, when they differ; mashtree can drop a degenerate genome and
  exit 0. The tree is kept for inspection and the stage is not recorded as
  completed.

### Fixed
- The `reads` no-match message names every active filter, including
  `--max-bases`, and the command reference shows list defaults (for example
  `reads --drop-selection`) as comma-separated values.
- Entrez taxonomy enrichment stops at the first batch of taxids that fails
  with a connection error (no HTTP status) and exits 3 with a message that
  names the host that gave no response. Before, every batch was
  retried three times, about 16 minutes for 1050 taxids with the network
  down. HTTP errors are still retried per batch.
- `glance` removes the previous dendrogram and plots once the comparison
  succeeds, so a plot with no similarity in the `--plot-min`/`--plot-max`
  range is absent instead of left from an earlier run. A dRep failure prints
  one console line; its command and output tail go to the run log (a test now
  covers this for glance).
- `derep-stock --action delete` of a run that is not stored exits 3 with a
  message naming the run and listing the stored runs. A repeat delete was
  skipped by the resume check and exited 0; delete is no longer recorded or
  skipped.
- `derep-unpack` warns about each cluster member that is missing from
  `genomes/` (one line listing them when there are more than ten) instead of
  leaving it out silently. `derep-stock --action unpack` of a stored run
  without `cluster_summary.tsv` rebuilds the summary from the restored
  clusters, and drops a live `genome_status.tsv` the stored run lacks,
  instead of keeping the files of the replaced dereplication.
- `assemble` and `assemble-run` with `--polisher auto` print one warning per
  platform when an adapter would polish the runs but its tool is not
  installed (for example medaka for ONT), naming the adapters and the remedy.
  The runs are still assembled, unpolished; before, this happened without a
  message.
- `reads --accession-file` treats text from a `#` to the end of the line as a
  comment, also after indentation or after an accession. `reads` (and `run
  --reads`) parse the accession file and check that a selection is given
  before the workdir is created, so a rejected call leaves no directory or
  log; entry stages in general build their parameters before creating the
  workdir.
- The `snippy` SNP typer names the reference record by its genome instead of
  snippy-core's `Reference`, so the tree has a leaf for the reference genome
  that `tree2tax` can map to an accession.
- `genome --accession-list-only` writes `ncbi_acc_download_list.txt` with a
  newline after the last accession, so `wc -l` and `while read` loops see
  every accession.
- `vgenome --outgroup-treebuilder` help names the accepted values, taken from
  the registered tree builders with distance-matrix support.
- `phylo` and `phylo-build` exit 3 with "A tree needs at least 3 genomes" when
  the ingroup has fewer than three genomes, before any aligner, SNP typer or
  tree builder runs, instead of failing inside the tree builder (for example
  quicktree on a single representative).
- A BV-BRC group download is written to a temporary file and renamed after
  the size check, so an interrupted transfer no longer leaves a partial
  `download.fa` that the next run reuses as complete. GTDB table downloads
  and the NCBI Virus path were already protected.
- `status` and `doctor` exit 3 with the shared "Workdir not found" message
  when `-wd` does not exist, instead of exiting 0; an existing directory
  without `repgenr.yaml` still prints the entry-stage hint (status) or a
  warning (doctor).
- Paired runs that ENA lists with a third, orphan FASTQ file (260 of the 424
  paired Illumina Wolbachia runs) no longer fail in skesa: the adapters find
  the `_1`/`_2` pair; skesa also takes the orphan file as unpaired input,
  shovill uses the pair only. The per-sample run choice prefers a long-read
  run only when it carries at least 100 Mb and a tenth of the sample's
  largest short-read run (a 45 kb PacBio run was chosen over 9 Gb of
  Illumina reads).
- Container mounts bind the real directory behind a symlinked path, so a
  containerised tool can open inputs that Nextflow staged through a linked
  directory (seen with sourmash `dereplicate-merge` under the docker
  backend). The stateless reads steps resolve their directories before any
  path reaches a tool.
- Segment-grouped viral isolates (`vgenome --group-segments`) reached
  `genomes_map.tsv` only as their synthetic `iso-` token; the member segment
  accessions were lost from the deliverable. `vgenome` now writes
  `segments.tsv` (token to member accessions) and `tree2tax` lists the
  members under the isolate's leaf.
- The `phylo` record in `repgenr.yaml` omitted the reference, the masker and
  the adapter tuning, so it could not tell a masked tree from an unmasked
  one; the three are recorded now (`reference`, `mask`, `extra`).
- `vgenome` (both back-ends) deleted `genomes/` before writing the new set
  and wrote `selection.tsv` last, so a crash mid-write left a partial set
  with no selection table and the next stage ran on it without a guard. The
  genomes are now built in a staging directory beside `genomes/` and swapped
  in only when every file is written; a failed run leaves the previous set
  and its table untouched.
- `tree2tax --root-name` was ignored unless `--remove-outgroup` was also
  given; the top node was always labelled `root`.
- The RAxML-NG tree builder published `.raxml.bestTree`, which carries no
  support values, so `tree2tax --collapse-support` never collapsed anything
  on its trees. It now publishes the support-annotated tree `--all` writes.
- dRep in virus mode passed `--S_algorithm` twice, so a `S_algorithm=`
  tool-arg was silently replaced by `ANImf`.
- `derep-stock --action unpack` restored the derep tables and representatives
  but left the `dereplicate` record with its previous fingerprint, the
  manifest with the previous per-genome status and `cluster_summary.tsv`
  absent. It now stores and restores the summary, refreshes the manifest and
  re-stamps the record without a fingerprint so the next `dereplicate` runs.
- The viral path (`vgenome`, both back-ends) never wrote the SQLite manifest;
  `dereplicate` then created an empty one, `--reduce species|genus` was a
  silent no-op and `doctor` reported every selected genome as missing.
- The contract TSVs (`clusters.tsv`, `genome_status.tsv`, `selection.tsv`,
  `tree2tax.tsv`, `genomes_map.tsv`) ended their rows with `\r\n`, the
  `csv` module's default, so the last column of every row carried a stray
  `\r`: `awk -F'\t' '$2=="x"'` on `clusters.tsv` never matched. They now
  end with `\n`. Readers were unaffected either way.
- Command audit (2026-10-07): `run --with-snptype --snptyper <unknown>` is
  rejected with exit 2 before any stage starts, and `run --with-snptype
  --mask gubbins` reaches the snptype stage instead of exiting 2.
- `versions` on a workdir without `repgenr.yaml` exits 3 and writes no
  fragment. `status` names `metadata`, `vmetadata`, `ingest` and `reads` as
  starting points, and the `-B/--bootstrap` help says what IQ-TREE needs.
- `list-tools` and the version records no longer store a tool's error line
  when the tool rejects its version flag (sibeliaz); the version is recorded
  as unknown.
- A missing workdir exits 3 on `dereplicate`, `tree2tax`, `glance`,
  `derep-unpack`, `cluster-summary`, `derep-stock`, `genome`, `vgenome` and
  `assemble`, instead of a traceback or a silently created workdir. `glance`
  reports a workdir without genomes before it checks for dRep, and records
  the dRep version.
- `genome-fetch` runs the `datasets` preflight without `--versions-out`
  (exit 4 instead of a traceback). A `selection.tsv` or `reads.tsv` that lacks
  required columns exits 3 and names them, instead of a KeyError.
- `excused_runs.tsv` keeps each excused run on one line; a multi-line failure
  reason no longer spans several physical lines.
- `snptype` drops optional outputs (`full_alignment.fasta`,
  `snp_distance_matrix.tsv`, `variants.vcf`) that the current typer does not write, and
  `phylo` clears the previous tree builder's files from `tree/` before a
  rebuild.
- The SibeliaZ wrapper on macOS skips AppleDouble `._*` files and leaves no
  empty block temp files in `align/`.
- `tree2tax` maps a dereplicated member that is also a leaf once, and does
  not log `--include-dereplicated` as a changed input file. A malformed Newick
  tree (`tree2tax`, `tree2tax-relations`) and a nonexistent `--clusters` path
  exit 3.
- `derep-stock`: `pack` refuses a workdir without dereplication outputs,
  `unpack` checks the stored run before it replaces anything, a repeat
  `unpack` reruns when the live dereplication changed, and every action exits
  3 on a nonexistent workdir.
- `docs/output.md` and `docs/usage.md` agree with the runs: `scratch/` is
  written by five stages, tool intermediates of `dereplicate` are under
  `scratch/`, and `glance` compares all genomes.
- `metadata` with a local GTDB table names `--gtdb-version` when it is missing
  (it named `--version`) and reruns when the `--metadata-path` table changes.
  `vmetadata --source bvbrc` reports an unreachable FTP server as exit 3.
- `snptype --tool simple` with no variable sites names the genome count and
  the reference, and suggests a closer reference, more divergent genomes, or
  an alignment-free tree.

### Changed
- A run whose assembler is not installed is excused as
  `assembler_not_installed` with a warning, and the other runs proceed. With
  `--assembler auto`, `assemble` exits 4 only when no run can be assembled.
- The `simple` SNP typer's core-SNP reduction is vectorised. It compared every
  pair of genomes character by character in Python, which cost about 8 ms per
  pair at 200000 sites and made the adapter's advertised limit of 2000 genomes
  unreachable. The output is byte for byte what it was; 40 genomes went from
  3.7 s to 0.13, and 500 genomes now take 35 s. The work is still quadratic in
  genomes, but each pair is now a vector comparison.
- Under a container backend, the `simple` typer runs each genome's chain of
  tools in one container instead of starting one per command. The live
  containerised test went from 21.2 s to 12.6 for 8 genomes, and its run log
  shows 9 engine invocations where it previously needed 58.
- The `simple` SNP typer runs genomes concurrently and passes a thread count to
  minimap2, samtools and bcftools; it previously mapped one genome at a time on
  one core whatever `--threads` said. On 69 Francisella genomes with 8 threads
  the typing stage went from 10.8 minutes to 2.8.
- The same typer writes its per-genome pileup and calls as compressed BCF and
  deletes each genome's intermediates once its consensus has been read. The
  Francisella run left 10 GB of scratch behind, 8.3 GB of it uncompressed
  pileup VCF; the same run now leaves 190 MB, nearly all of it the whole-genome
  alignment the stage publishes.
- `phylo` stamps the alignment it builds and reuses it when a later run changes
  only the tree builder, the bootstrap or the thread count. Trying a second
  tree builder repeated the whole alignment or SNP-calling step before, which
  on that same set was 11 minutes per attempt. The stamp records the genome
  set, the source settings and the alignment's digest, so any change to those,
  or `--force`, rebuilds it.

### Added
- The Nextflow layer can run the alignment and the tree as separate tasks
  (`--phylo_split_msa`, off by default). The phylogeny was one task, so
  changing the tree builder or the bootstrap repeated the alignment or the SNP
  calling, and both halves shared one resource label. The split gives the
  alignment its own cache entry and `process_high`, the tree `process_medium`.
- `phylo-build --msa-only` builds the alignment and stops, writing
  `msa.fasta`; `phylo-build --msa <file>` builds the tree from an alignment an
  earlier call produced. These are the two halves the Nextflow processes run.
- The Gubbins masker reports how much of the alignment is variable, warns above
  10%, and repeats the figure when Gubbins fails, instead of leaving a bare
  non-zero exit. Measured on Francisella subsets, a set of one species (1-8%
  variable) is masked in minutes while a two-species or genus-level set
  (13-39%) crashes Gubbins' recombination scan.
- `--mask gubbins` reads `--tool-arg gubbins_tree_builder=...`,
  `gubbins_first_tree_builder=...` and `gubbins_args="..."` (Gubbins'
  `--tree-builder`, `--first-tree-builder` and any further arguments). The
  masker's keys count as read in the unread-extras warning of the snptype
  and phylo stages.

### Fixed
- IQ-TREE refused a recombination-masked alignment with "Unknown sequence
  type": masking replaces recombinant bases with N, and a fifth of the
  characters was enough to defeat its guess. The adapter now states `-st DNA`,
  which every alignment this pipeline produces is.
- Gubbins exited before its first iteration on hosts whose RAxML package has
  no multi-threaded build (`raxmlHPC-PTHREADS*`; the osx-arm64 conda build
  among them) whenever `--threads` was above one, and the only trace was one
  stderr line. A native run without such a build now uses IQ-TREE as the
  Gubbins tree builder, or one thread when IQ-TREE is missing too, and logs
  the switch.
- Clearing a scratch or output directory on a volume without native extended
  attributes (an external exFAT disk) could fail with "No such file or
  directory": macOS drops a file's AppleDouble twin the moment the data file
  goes, and the removal walk then tried to unlink the vanished name. Every
  such removal now tolerates entries that disappear mid-walk (found on a
  genus-scale run whose workdir sat on an external volume).

### Documentation
- Help and wording from the command audit: `glance --tool` lists only
  dereplicators that support comparison, the `reads` taxon options say only the
  most specific is used, `tree2tax --node-basename` says internal nodes
  otherwise get hash-derived names, and `--metadata-path` says `-r` and
  `--gtdb-version` are still required. usage.md notes that variable-site-only
  alignments give inflated branch lengths and that `phylo-build --msa-only`
  leaves `snp/` in place; output.md notes the `ska2` k-mer files in
  `scratch/snptype/`.
- `docs/install.md` covers installation: the package and its extras, three
  ways to provide the external tools (one conda environment, several
  environments on `PATH`, containers), a recommendation per situation, a
  per-tool table with conda packages, container pins and databases, Cactus,
  Nextflow and how to verify an installation. The container image, storage and
  notes sections moved there from `docs/usage.md`.
- `docs/choosing-tools.md` helps choose tools by dataset type and size, with the
  measured runs behind each recommendation, the dereplicator and tree builder
  trade-offs, assemblers and polishers per platform, `run` versus Nextflow, and
  how `auto` chooses. It states which limits are declared and not benchmarked.
  A new unit test checks its table of declared genome limits against the
  registered adapters.
- The README installation section points to the install guide, its pipeline
  section shows the four stage chains, and its scalability claim states the
  largest recorded runs (1157 bacterial and 1256 viral genomes).
- `docs/output.md` lists `derep/unpacked/`, `derep/stock/<name>/`, the
  `glance` plots and `glance_wd/`, and `ncbi_acc_download_list.txt`. The
  exit-code table in `docs/usage.md` states that `assemble` and `reads-gather`
  exit 3 when every run was excused and nothing was produced.

## [2.1.0] - 2026-09-10

The first public release of the fork. It carries the 2.0.0 rewrite (kept
below as its own section, never published) plus everything since: the
quality keeper, ska2, the nf-core Nextflow layer, and the 2026-09-10
pipeline audit with its live suite.

### Added
- `run` exposes the stage flags it used to leave at their defaults (D-3):
  `--limit`; `--complete-only`, `--host`, `--released-after` for the viral
  chain; `--process-size`, `--num-processes`, `--reduce`, `--target-reps`,
  `--tool-arg`, `--allow-incomplete` for dereplication; `--all-genomes`,
  `--bootstrap`, `--reference`, `--aligner-arg`, `--mask` for the
  phylogeny; `--collapse-support` and `--collapse-length` for tree2tax.
  Each goes through the same builder as the per-stage command, so
  validation and the resume fingerprint are identical.
- `phylo-build --mask` (D-6): the stateless step accepts the same
  recombination masker as `phylo` for `--msa-source snptype`, so the
  Nextflow layer can mask through `phylo_args`.
- `docs/cli-reference.md`: every command with its options, defaults and
  help, generated from the command tree by `scripts/render_cli_matrix.py`
  and kept in sync by a test; a second test checks that every flag is
  mentioned in at least one document besides the audit matrix.
  `docs/verification.md` is now a how-to plus a results table that
  `scripts/live_report.py` renders from the live suite's junit output.
- The CLI matrix carries a `nextflow` column: for every flag, the
  `params.*` key or `ext.args` string that reaches the Nextflow layer, or
  the reason it does not; a test checks each key against
  `nextflow_schema.json`. A stub test asserts the dereplication `ext.args`
  composed from `params.derep_*` and the viral mode reaches the process.
  Live Nextflow tests (`tests/live/test_nextflow.py`): the local
  data-channel harness with two dereplicator and tree-builder pairs,
  `main.nf` in bacterial and viral mode, and the docker profile with
  progressiveMauve.
- Container live tests (`-m "live and container"`,
  `tests/live/test_container_runs.py`): skder, sourmash, dRep (`--virus`),
  the simple SNP typer and sibeliaz in Wave-minted images with
  `--container-cache`, the `REPGENR_CONTAINER*` variables and the rule that
  a native result is never reused by a container run; progressiveMauve and
  Cactus in their pinned amd64 images; `glance --tool drep`; the
  `--container-engine podman` negative path. The module skips itself without
  Docker, amd64 emulation and the Wave CLI.
- `dereplicate-chunk --selection-tsv --keeper` covered live.
- Live tests on the Francisella tularensis species set
  (`tests/live/test_species_set.py`): the simple, ska2 and parsnp SNP typers
  with `--reference`, `--all-genomes`, `--tool-arg` and `--allow-incomplete`;
  `--mask gubbins` changing the core alignment (and refused on ska2);
  IQ-TREE with `-B 1000` support labels and `--no-outgroup`, FastTree and
  RAxML-NG from the SNP source, `--mask` through `phylo`; every workdir
  `tree2tax` flag.
- Network live tests (`-m "live and network"`, `tests/live/test_network.py`):
  GTDB API selections at genus, species (`--limit`, `--outgroup-accession`)
  and family level; the GTDB TSV table with `--nodownload` and
  `--metadata-path`; `genome --accession-list-only` and `--keep-files`; NCBI
  Virus with `--complete-only`, `--released-after`, `--host`; every `vgenome`
  selection flag; BV-BRC `--list`, `--source bvbrc` and `--filter`; `run
  --dry-run` and the bacterial and viral chains end to end. Downloads are
  built once under the live cache directory and reused.
- Live tests for the offline stages (`tests/live/`, 29 tests, about three
  minutes): skder, sourmash and galah recover the synthetic partition; every
  `dereplicate` flag has an observable effect; `dereplicate-chunk` x3 plus
  `dereplicate-merge` by directory and by file list; `phylo-build` and
  `tree2tax-relations` with an outgroup, the collapse thresholds and
  `--versions-out`; the mashtree and sourmash tree builders; `derep-unpack`,
  `derep-stock`, `status`, `versions`, `doctor`; `--force`, `--quiet`,
  `--verbose`, `REPGENR_FORCE` and `REPGENR_LOG_LEVEL`; `ingest --outgroup`,
  `--copy` and `--selection`.
- A CLI matrix, `tests/audit/cli_matrix.yaml`, with one record per command
  and flag (aliases, the parameter it sets, how it is validated, the docs
  that mention it, its live test). `tests/unit/test_cli_matrix.py` checks the
  records against the real command tree on every run: flag and alias sets,
  help text, that each flag reaches its parameter, that a bad value is
  rejected naming the flag, and that live and docs references resolve.
  `scripts/render_cli_matrix.py` renders `docs/audit/cli-matrix.md`.
- `list-tools` now lists the maskers family.
- `repgenr ingest -wd WD --genomes-dir DIR [--selection TSV] [--outgroup NAME|FILE]
  [--copy]`: start a working directory from genomes already on disk. It
  writes the `selection.tsv` and manifest the download stages produce, links
  (or copies) the files into `genomes/`, stages an optional outgroup, records
  itself as a stage (resume, `status` with the local chain, `doctor`), and
  takes taxonomy and CheckM quality from the selection table or the canonical
  filename.
- A live verification suite under `tests/live/` (markers `live`, `network`,
  `container`; deselected by default via `addopts`; `--live-config` maps
  tools to bin directories). The first test runs ingest, sourmash, mashtree
  and tree2tax on a seeded synthetic set and checks the recovered partition
  against the generator's truth.
- `tree2tax --collapse-length L` and `--collapse-support S` (also on the
  `tree2tax-relations` step): internal nodes whose branch is shorter than L,
  or whose support is below the fraction S (percentage trees are normalised),
  merge into their parent before nodes are named, so weak splits do not
  become FlexTaxD nodes. The root, the outgroup/ingroup split and leaves
  never collapse. Provenance records the thresholds and the collapsed count.
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
- Recombination masking with an outgroup (D-10): Gubbins now scans the
  ingroup only and its predicted regions are masked in the whole-genome
  alignment for every genome, outgroup included, before the variable
  sites are extracted. A species-level outgroup used to abort Gubbins'
  scan ("gubbins ... exit 1"), which forced `--no-outgroup` on masked
  trees. The masked alignment for a run without an outgroup is Gubbins'
  own filtered file, as before.
- The Nextflow `standard` profile caps process resource requests to the
  host's CPUs and memory (D-11). A laptop run used to fail with "Process
  requirement exceeds available CPUs -- req: 32" unless `-profile test`
  was given; `slurm` and `cloud` stay uncapped.
- `glance`, `derep-unpack` and `derep-stock` (pack, unpack, delete) record
  a stage like every other command (D-7): `status` lists them under the
  optional stages, `doctor` sees them, and an identical repeat skips on
  the resume records (genomes and clusters are digested as inputs).
  `derep-stock --action list` stays a query and records nothing.
- `vmetadata --filter` is a BV-BRC option (D-4): it has no default on the
  command line, BV-BRC applies "complete genome" when it is unset, and
  passing it with the NCBI Virus source is an error that points at
  `--complete-only`. It used to be accepted and silently ignored there.
- Short options mean one thing everywhere (D-2): `-t` is `--threads` on
  every command that has threads, including `run`; `--target` on
  `vmetadata` and `run`, `--filter` and `--list` on `vmetadata`, and
  `--root-name` on `tree2tax` and `tree2tax-relations` have no short form
  any more (`-t`, `-f`, `-l` and `-r` used to collide with `--threads`,
  `--force`, `--level` and `--release`). Nextflow's `vmetadata_args`
  default and the run id derivation use `--target`.
- `tree2tax-relations --include-dereplicated` is on by default, as on
  `tree2tax` and `run`; `--no-include-dereplicated` turns it off. The
  Nextflow `tree2tax_args` default is therefore empty (D-1).
- With `--container` but no `--wave`, an adapter that only declares a conda
  spec runs on the host; the warning now says so and names the remedy
  (`pass --wave`), instead of claiming the tool declares no image.
- `Tree2taxParams.all_genomes`, a field no code read, is gone. It was part
  of the tree2tax resume fingerprint, so an existing workdir re-runs
  `tree2tax` once (seconds). `GlanceParams.threads` defaults to 16 like the
  CLI, and the three aligners take the MSA filename from the shared
  contract constant. Synthetic sets at very small n no longer carry an
  empty cluster. `benchmarks/run_bench.py` reads its storage root from
  `REPGENR_BENCH_STORAGE` (or `--storage`) instead of a fixed volume path.
- Closed-choice options are validated when the command is parsed rather
  than deep in the stage: `metadata -d/-l/--source` (and the same values
  under `run --dataset/--level/--metadata-source/--viral-source`, reported
  under those names), `vgenome --length-method` and `--outgroup-treebuilder`
  (any tree builder with distance-matrix support), `glance --tool`,
  `derep-stock --action`, and `--mask` on `snptype` and `phylo`, which is
  now checked against the masker registry instead of a fixed list. `phylo
  --mask` with `--msa-source aligner` is an error; before, the mask was
  silently dropped. `vgenome_params` has an explicit signature, so a
  misspelt keyword is a `TypeError` instead of a dataclass error.
- Every command-line option has help text (48 were blank: the ANI thresholds
  and `--threads` on the dereplication commands, the GTDB target flags, most
  `vgenome` selection flags, `glance` plot bounds, `derep-unpack
  --no-representant`). Shared texts live in `cli/base.py` so the same flag
  reads the same on every command. `run --snptyper` lists the registered SNP
  typers like `snptype --tool` does.
- `metadata --limit N` no longer keeps the first N genomes in GTDB file (or
  API) order. It round-robins over species, taking the best CheckM-scored
  genome of every species first, then each species' next best, until N;
  within a species unscored genomes rank last, then the GTDB representative
  flag, then accession. The same `--limit` therefore returns a different,
  better set than before. On the API path the per-genome quality cards are
  fetched for every candidate before the cut.
- **Nextflow layer in nf-core shape.** Every channel carries a run-level meta
  map (`id` from the selection target, `mode`), processes exchange
  `tuple val(meta), path(...)`, tool flags reach processes as
  `task.ext.args` mapped from the unchanged parameters in
  `nextflow/conf/modules.config` (where publishing now lives too), and
  resources and the retry window sit in `nextflow/conf/base.config`. The
  Nextflow floor is 26.04 (`!>=26.04.0`) with nf-schema 2.6.1, and the layer
  lints clean under the strict parser. Command lines and parameters are
  unchanged; per-stage `versions.yml` fragments are no longer copied under
  `--outdir` (the collected `pipeline_info/software_versions.yml` remains);
  code that included a RepGenR module or subworkflow directly must
  adapt to the tuple shapes.
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
- Cactus wrote Toil's `.toil/` state under the container's HOME, which
  the backend points at the working directory; the aligner now runs in its
  alignment directory, so nothing lands in the launch directory.
- Containers started in the temp mount when the adapter gave no working
  directory, so a stateless step run with relative paths (`phylo-build -o .`
  under the docker profile) could not open its outputs
  (`align/xmfa/...`). The container now starts in the host process's working
  directory, which is also bound.
- Nextflow: a run without an outgroup (viral `--no-outgroup`, or a GTDB
  selection with no outgroup candidate) aborted in the data-channel
  subworkflows with "Invalid method invocation `call`": the optional
  outgroup join emitted a bare meta map. The join now carries a placeholder
  tuple; the VACQUIRE stub honours `--no-outgroup` so a stub test covers it.
- mashtree received every genome path on its command line and failed with
  "Argument list too long" at about 9500 genomes (a viral `run` without a
  completeness filter). It now reads the paths from a file-of-files, and
  the genome directories are declared as container mounts.
- Container runs on a workdir populated by `repgenr ingest` (symlinked
  genomes) failed inside the container with dangling links: only the
  `genomes/` directory was bound, not the directories the links point to.
  Every mounted directory's symlink targets are bound as well.
- The simple SNP typer called variants with bcftools' diploid default, so
  heterozygous calls on haploid bacteria became IUPAC codes in the consensus
  and Gubbins refused the whole-genome alignment ("contains disallowed
  characters"). Calling is haploid now (`--ploidy 1`), and the Gubbins
  masker replaces any remaining non-ACGTN symbol with N before running.
- `snptype` without `--reference` picked the alphabetically first genome
  silently on the workdir path; it now logs the same warning the stateless
  path always did, naming the genome it chose.
- `phylo --msa-source snptype` typed only the ingroup, so the outgroup never
  reached the SNP alignment and the tree could not be rooted on it. The
  outgroup is now typed with the ingroup, as on the aligner path.
- `phylo --treebuilder fasttree --bootstrap N` was silently ignored; N is
  now FastTree's `-boot` resample count for its local support values.
- `vgenome --group-segments` concatenated every set of records that shared
  an isolate name, segmented or not; on hepatovirus (one segment) it turned
  367 complete genomes into 149 (216 records share the isolate name "RNA"
  and every one is labelled segment "ANONYMOUS"). An isolate is grouped
  only when its records carry at least two distinct real segment labels.
- `vgenome` on the BV-BRC source wrote the genomes but no `selection.tsv`,
  unlike the NCBI Virus path and the bacterial stages, so the downstream
  contract was incomplete. It now publishes the same table (taxonomy from
  the Entrez names, outgroup row included).
- The genome stage logged "Discarding non-FASTA download for <accession>
  (error page?)" for every accession after a successful batch on macOS
  volumes without native extended attributes: the extraction loop picked up
  the AppleDouble `._*.fna` twins. They are skipped now.
- `vmetadata --source bvbrc --target hepatitis_e_virus` failed with "Could
  not find virus group": the group name was capitalised on its first letter
  only, but BV-BRC names groups like `Hepatitis_E_virus`. The name is now
  resolved against the server listing, ignoring case.
- `--versions-out` on the stateless steps failed with "No such file or
  directory" when the fragment path sat in the not-yet-created output
  directory; the parent is now created first.
- Six small defects noted in the 2026-09-01 audit's self-review: `phylo`
  publishes `tree/tree.nwk` through the atomic copy used by every other
  deliverable; the Nextflow retry window no longer includes exit 130 and 131
  (a cancelled task is not resubmitted); an unknown `--container` value is a
  user-input error with the valid choices in the message; IQ-TREE refuses
  `--bootstrap` values below its floor of 1000 before running; the Wave image
  cache is keyed by platform as well as conda spec; and a failure inside the
  tool-output read loop kills the child process before the error propagates.
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

### Notes
- With the 26.04 floor, moving version reporting to `eval()` outputs on the
  `versions` topic costs nothing in compatibility and is the natural next step
  for the Nextflow layer.

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
