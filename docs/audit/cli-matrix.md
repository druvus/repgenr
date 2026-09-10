# CLI matrix

Generated from `tests/audit/cli_matrix.yaml` by `scripts/render_cli_matrix.py`;
`tests/unit/test_cli_matrix.py` keeps both in step with the command tree.

22 commands, 197 flags (19 with a live test or an n/a reason, 178 pending).

## Global flags

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--version` |  | n/a: prints the version and exits | none | todo: PR-E | docs/adding-tools.md |
| `--container` |  | container.backend | choice | todo: PR-H | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--container-engine` |  | container.engine | none | todo: PR-H | docs/containers.md |
| `--container-cache` |  | container.cache_dir | none | todo: PR-H | docs/containers.md, docs/verification.md |
| `--platform` |  | container.platform | none | todo: PR-H | README.md, docs/containers.md |
| `--wave` |  | container.wave_enabled | none | todo: PR-H | docs/adding-tools.md, docs/containers.md, docs/swot-viral.md |
| `--force` | -f | state.force | none | todo: PR-E | README.md |
| `--verbose` | -v | state.log_level | none | todo: PR-E | README.md |
| `--quiet` | -q | state.log_level | none | todo: PR-E |  |

## derep-stock

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--action` |  | DerepStockParams.action | choice | todo: PR-E |  |
| `--name` |  | DerepStockParams.name | none | todo: PR-E | docs/usage.md |

## derep-unpack

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--no-representant` |  | DerepUnpackParams.no_representant | none | todo: PR-E |  |

## dereplicate

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--tool` |  | DereplicateParams.tool | registry | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` | -pani | DereplicateParams.primary_ani | unit_interval | todo: PR-E | docs/verification.md |
| `--secondary-ani` | -sani | DereplicateParams.secondary_ani | unit_interval | todo: PR-E | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` | -af | DereplicateParams.aligned_fraction | unit_interval | todo: PR-E |  |
| `--threads` | -t | DereplicateParams.threads | range | todo: PR-E | docs/verification.md |
| `--process-size` | -s | DereplicateParams.process_size | none | todo: PR-E | docs/architecture.md, docs/scaling-audit.md, docs/swot-derep.md, docs/verification.md |
| `--num-processes` | -p | DereplicateParams.num_processes | none | todo: PR-E | docs/verification.md |
| `--pre-primary-ani` |  | DereplicateParams.pre_primary_ani | unit_interval | todo: PR-E | docs/verification.md |
| `--pre-secondary-ani` |  | DereplicateParams.pre_secondary_ani | unit_interval | todo: PR-E | docs/verification.md |
| `--reduce` |  | DereplicateParams.reduce | choice | todo: PR-E | docs/scaling-audit.md, docs/swot-derep.md |
| `--target-reps` |  | DereplicateParams.target_reps | range | todo: PR-E | docs/scaling-audit.md, docs/swot-derep.md |
| `--virus` |  | DereplicateParams.extra | none | todo: PR-E | README.md, docs/adding-tools.md, docs/swot-derep.md, docs/verification.md |
| `--tool-arg` |  | DereplicateParams.extra | callback | todo: PR-E | docs/adding-tools.md, docs/usage.md |
| `--allow-incomplete` |  | DereplicateParams.allow_incomplete | none | todo: PR-E |  |
| `--keeper` |  | DereplicateParams.keeper | choice | todo: PR-E | README.md, docs/swot-derep.md, docs/usage.md |

## dereplicate-chunk

dispatch: `step:repgenr.stages.derep_steps.dereplicate_chunk`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--genomes-fofn` |  | ChunkParams.genomes | stage | todo: PR-E |  |
| `--out` | -o | ChunkParams.out_dir | none | todo: PR-E |  |
| `--tool` |  | ChunkParams.tool | registry | todo: PR-E | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` | -pani | ChunkParams.primary_ani | unit_interval | todo: PR-E | docs/verification.md |
| `--secondary-ani` | -sani | ChunkParams.secondary_ani | unit_interval | todo: PR-E | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` | -af | ChunkParams.aligned_fraction | unit_interval | todo: PR-E |  |
| `--threads` | -t | ChunkParams.threads | range | todo: PR-E | docs/verification.md |
| `--virus` |  | ChunkParams.extra | none | todo: PR-E | README.md, docs/adding-tools.md, docs/swot-derep.md, docs/verification.md |
| `--tool-arg` |  | ChunkParams.extra | callback | todo: PR-E | docs/adding-tools.md, docs/usage.md |
| `--selection-tsv` |  | ChunkParams.selection_tsv | none | todo: PR-E |  |
| `--keeper` |  | ChunkParams.keeper | choice | todo: PR-E | README.md, docs/swot-derep.md, docs/usage.md |
| `--versions-out` |  | ChunkParams.versions_out | none | todo: PR-E |  |

## dereplicate-merge

dispatch: `step:repgenr.stages.derep_steps.dereplicate_merge`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--out` | -o | MergeParams.out_dir | none | todo: PR-E |  |
| `--chunk-dir` |  | MergeParams.chunk_dirs | none | todo: PR-E |  |
| `--chunk-fofn` |  | MergeParams.chunk_dirs | stage | todo: PR-E |  |
| `--tool` |  | MergeParams.tool | registry | todo: PR-E | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` | -pani | MergeParams.primary_ani | unit_interval | todo: PR-E | docs/verification.md |
| `--secondary-ani` | -sani | MergeParams.secondary_ani | unit_interval | todo: PR-E | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` | -af | MergeParams.aligned_fraction | unit_interval | todo: PR-E |  |
| `--threads` | -t | MergeParams.threads | range | todo: PR-E | docs/verification.md |
| `--virus` |  | MergeParams.extra | none | todo: PR-E | README.md, docs/adding-tools.md, docs/swot-derep.md, docs/verification.md |
| `--tool-arg` |  | MergeParams.extra | callback | todo: PR-E | docs/adding-tools.md, docs/usage.md |
| `--selection-tsv` |  | MergeParams.selection_tsv | none | todo: PR-E |  |
| `--keeper` |  | MergeParams.keeper | choice | todo: PR-E | README.md, docs/swot-derep.md, docs/usage.md |
| `--versions-out` |  | MergeParams.versions_out | none | todo: PR-E |  |

## doctor

dispatch: `query`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | todo: PR-E | docs/containers.md, docs/output.md |

## genome

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--accession-list-only` |  | GenomeParams.accession_list_only | none | todo: PR-F |  |
| `--keep-files` |  | GenomeParams.keep_files | none | todo: PR-F |  |

## genome-fetch

dispatch: `step:repgenr.stages.genome_steps.genome_fetch`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--selection` |  | GenomeFetchParams.selection_tsv | stage | todo: PR-F | README.md, docs/usage.md |
| `--out` | -o | GenomeFetchParams.out_dir | none | todo: PR-F |  |
| `--keep-files` |  | GenomeFetchParams.keep_files | none | todo: PR-F |  |
| `--versions-out` |  | GenomeFetchParams.versions_out | none | todo: PR-E |  |

## glance

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--tool` |  | GlanceParams.tool | registry | todo: PR-H | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--threads` | -t | GlanceParams.threads | range | todo: PR-H | docs/verification.md |
| `--plot-max` |  | GlanceParams.plot_max | none | todo: PR-H |  |
| `--plot-min` |  | GlanceParams.plot_min | none | todo: PR-H |  |
| `--keep-files` |  | GlanceParams.keep_files | none | todo: PR-H |  |

## ingest

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--genomes-dir` |  | IngestParams.genomes_dir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | README.md, docs/usage.md |
| `--selection` |  | IngestParams.selection | stage | todo: PR-E | README.md, docs/usage.md |
| `--outgroup` |  | IngestParams.outgroup | stage | todo: PR-E | README.md, docs/usage.md |
| `--copy` |  | IngestParams.copy | none | todo: PR-E | README.md, docs/usage.md |

## list-tools

dispatch: `query`

(no flags)

## metadata

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--dataset` | -d | MetadataParams.dataset | choice | todo: PR-F |  |
| `--level` | -l | MetadataParams.level | choice | todo: PR-F |  |
| `--source` |  | MetadataParams.source | choice | todo: PR-F | README.md, docs/usage.md, docs/verification.md |
| `--release` | -r | MetadataParams.release | none | todo: PR-F | README.md |
| `--gtdb-version` |  | MetadataParams.version | none | todo: PR-F | README.md, docs/usage.md |
| `--target-family` | -tf | MetadataParams.target_family | none | todo: PR-F |  |
| `--target-genus` | -tg | MetadataParams.target_genus | none | todo: PR-F | README.md |
| `--target-species` | -ts | MetadataParams.target_species | none | todo: PR-F |  |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | todo: PR-F |  |
| `--metadata-path` |  | MetadataParams.metadata_path | none | todo: PR-F |  |
| `--nodownload` |  | MetadataParams.nodownload | none | todo: PR-F |  |
| `--limit` |  | MetadataParams.limit | range | todo: PR-F | docs/scaling-audit.md, docs/swot-derep.md, docs/usage.md |

## phylo

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | todo: PR-G | README.md, docs/containers.md, docs/usage.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | todo: PR-G | README.md, docs/usage.md |
| `--aligner` |  | PhyloParams.aligner | registry | todo: PR-H | README.md, docs/containers.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | todo: PR-G |  |
| `--all-genomes` |  | PhyloParams.all_genomes | none | todo: PR-G |  |
| `--no-outgroup` |  | PhyloParams.no_outgroup | none | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree |  |
| `--bootstrap` | -B | PhyloParams.bootstrap | range | todo: PR-G | docs/usage.md |
| `--reference` |  | PhyloParams.reference | none | todo: PR-G | docs/swot-phylo.md |
| `--aligner-arg` |  | PhyloParams.extra | callback | todo: PR-H | docs/adding-tools.md |
| `--threads` | -t | PhyloParams.threads | range | todo: PR-G | docs/verification.md |
| `--mask` |  | PhyloParams.extra | registry | todo: PR-G | README.md, docs/adding-tools.md, docs/output.md, docs/usage.md, docs/verification.md |
| `--allow-incomplete` |  | PhyloParams.allow_incomplete | none | todo: PR-G |  |

## phylo-build

dispatch: `step:repgenr.stages.phylo.phylo_build`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--genomes-dir` |  | PhyloBuildParams.genomes_dir | stage | todo: PR-E | README.md, docs/usage.md |
| `--out` | -o | PhyloBuildParams.out_dir | none | todo: PR-E |  |
| `--outgroup-dir` |  | PhyloBuildParams.outgroup_dir | none | todo: PR-E |  |
| `--outgroup-accession` |  | PhyloBuildParams.outgroup_accession | none | todo: PR-E |  |
| `--treebuilder` |  | PhyloBuildParams.phylo.treebuilder | registry | todo: PR-E | README.md, docs/containers.md, docs/usage.md |
| `--msa-source` |  | PhyloBuildParams.phylo.msa_source | choice | todo: PR-E | README.md, docs/usage.md |
| `--aligner` |  | PhyloBuildParams.phylo.aligner | registry | todo: PR-H | README.md, docs/containers.md |
| `--snptyper` |  | PhyloBuildParams.phylo.snptyper | registry | todo: PR-E |  |
| `--no-outgroup` |  | PhyloBuildParams.phylo.no_outgroup | none | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree |  |
| `--bootstrap` | -B | PhyloBuildParams.phylo.bootstrap | range | todo: PR-E | docs/usage.md |
| `--reference` |  | PhyloBuildParams.phylo.reference | none | todo: PR-E | docs/swot-phylo.md |
| `--aligner-arg` |  | PhyloBuildParams.phylo.extra | callback | todo: PR-H | docs/adding-tools.md |
| `--threads` | -t | PhyloBuildParams.phylo.threads | range | todo: PR-E | docs/verification.md |
| `--versions-out` |  | PhyloBuildParams.versions_out | none | todo: PR-E |  |

## run

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--viral` |  | VmetadataParams | none | todo: PR-F | README.md |
| `--dataset` | -d | MetadataParams.dataset | choice | todo: PR-F |  |
| `--level` | -l | MetadataParams.level | choice | todo: PR-F |  |
| `--target-family` | -tf | MetadataParams.target_family | none | todo: PR-F |  |
| `--target-genus` | -tg | MetadataParams.target_genus | none | todo: PR-F | README.md |
| `--target-species` | -ts | MetadataParams.target_species | none | todo: PR-F |  |
| `--release` | -r | MetadataParams.release | none | todo: PR-F | README.md |
| `--gtdb-version` |  | MetadataParams.version | none | todo: PR-F | README.md, docs/usage.md |
| `--metadata-source` |  | MetadataParams.source | choice | todo: PR-F |  |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | todo: PR-F |  |
| `--target` | -t | VmetadataParams.target | none | todo: PR-F | README.md |
| `--viral-source` |  | VmetadataParams.source | choice | todo: PR-F |  |
| `--group-segments` |  | VgenomeParams.group_segments | none | todo: PR-F | README.md |
| `--tool` |  | DereplicateParams.tool | registry | todo: PR-F | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` |  | DereplicateParams.primary_ani | unit_interval | todo: PR-F | docs/verification.md |
| `--secondary-ani` |  | DereplicateParams.secondary_ani | unit_interval | todo: PR-F | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` |  | DereplicateParams.aligned_fraction | unit_interval | todo: PR-F |  |
| `--keeper` |  | DereplicateParams.keeper | choice | todo: PR-F | README.md, docs/swot-derep.md, docs/usage.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | todo: PR-F | README.md, docs/containers.md, docs/usage.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | todo: PR-F | README.md, docs/usage.md |
| `--aligner` |  | PhyloParams.aligner | registry | todo: PR-F | README.md, docs/containers.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | todo: PR-F |  |
| `--no-outgroup` |  | PhyloParams.no_outgroup | none | todo: PR-F |  |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | todo: PR-F | README.md, docs/usage.md |
| `--threads` |  | DereplicateParams.threads | range | todo: PR-F | docs/verification.md |
| `--dry-run` |  | n/a: prints the chain and exits before any stage | none | todo: PR-F |  |

## snptype

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--tool` |  | SnptypeParams.tool | registry | todo: PR-G | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--reference` |  | SnptypeParams.reference | none | todo: PR-G | docs/swot-phylo.md |
| `--all-genomes` |  | SnptypeParams.all_genomes | none | todo: PR-G |  |
| `--mask` |  | SnptypeParams.mask | registry | todo: PR-G | README.md, docs/adding-tools.md, docs/output.md, docs/usage.md, docs/verification.md |
| `--threads` | -t | SnptypeParams.threads | range | todo: PR-G | docs/verification.md |
| `--tool-arg` |  | SnptypeParams.extra | callback | todo: PR-G | docs/adding-tools.md, docs/usage.md |
| `--allow-incomplete` |  | SnptypeParams.allow_incomplete | none | todo: PR-G |  |

## status

dispatch: `query`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |

## tree2tax

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--node-basename` |  | Tree2taxParams.node_basename | none | todo: PR-E |  |
| `--root-name` | -r | Tree2taxParams.root_name | none | todo: PR-E |  |
| `--remove-outgroup` |  | Tree2taxParams.remove_outgroup | none | todo: PR-E |  |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | README.md, docs/usage.md |
| `--collapse-support` |  | Tree2taxParams.collapse_support | range | todo: PR-E | docs/swot-phylo.md, docs/usage.md |
| `--collapse-length` |  | Tree2taxParams.collapse_length | range | todo: PR-E | docs/swot-phylo.md, docs/usage.md |

## tree2tax-relations

dispatch: `step:repgenr.stages.tree2tax.tree2tax_relations`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--tree` |  | Tree2taxStepParams.tree | stage | todo: PR-E |  |
| `--out` | -o | Tree2taxStepParams.out_dir | none | todo: PR-E |  |
| `--clusters` |  | Tree2taxStepParams.clusters | none | todo: PR-E |  |
| `--outgroup-dir` |  | Tree2taxStepParams.outgroup_dir | none | todo: PR-E |  |
| `--outgroup-accession` |  | Tree2taxStepParams.outgroup_accession | none | todo: PR-E |  |
| `--node-basename` |  | Tree2taxStepParams.node_basename | none | todo: PR-E |  |
| `--root-name` | -r | Tree2taxStepParams.root_name | none | todo: PR-E |  |
| `--remove-outgroup` |  | Tree2taxStepParams.remove_outgroup | none | todo: PR-E |  |
| `--include-dereplicated` |  | Tree2taxStepParams.include_dereplicated | none | todo: PR-E | README.md, docs/usage.md |
| `--versions-out` |  | Tree2taxStepParams.versions_out | none | todo: PR-E |  |
| `--collapse-support` |  | Tree2taxStepParams.collapse_support | range | todo: PR-E | docs/swot-phylo.md, docs/usage.md |
| `--collapse-length` |  | Tree2taxStepParams.collapse_length | range | todo: PR-E | docs/swot-phylo.md, docs/usage.md |

## versions

dispatch: `query`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | todo: PR-E | docs/containers.md, docs/output.md |
| `--versions-out` |  | n/a: query command, no stage parameters | none | todo: PR-E |  |

## vgenome

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--target-genus` | -tg | VgenomeParams.target_genus | none | todo: PR-F | README.md |
| `--target-species` | -ts | VgenomeParams.target_species | none | todo: PR-F |  |
| `--target-serotype` | -tse | VgenomeParams.target_serotype | none | todo: PR-F |  |
| `--target-custom` | -tc | VgenomeParams.target_custom | none | todo: PR-F |  |
| `--length-all` |  | VgenomeParams.length_all | none | todo: PR-F | docs/usage.md |
| `--length-deviation` |  | VgenomeParams.length_deviation | range | todo: PR-F |  |
| `--length-method` |  | VgenomeParams.length_method | choice | todo: PR-F | docs/scaling-audit.md, docs/swot-viral.md, docs/usage.md |
| `--length-range` |  | VgenomeParams.length_range | none | todo: PR-F | docs/usage.md |
| `--discard` |  | VgenomeParams.discard | none | todo: PR-F |  |
| `--no-outgroup` |  | VgenomeParams.no_outgroup | none | todo: PR-F |  |
| `--group-segments` |  | VgenomeParams.group_segments | none | todo: PR-F | README.md |
| `--outgroup-candidates-taxid-min-genomes` |  | VgenomeParams.outgroup_candidates_taxid_min_genomes | none | todo: PR-F |  |
| `--outgroup-treebuilder` |  | VgenomeParams.outgroup_treebuilder | registry | todo: PR-F |  |
| `--glance` |  | VgenomeParams.glance | none | todo: PR-F |  |
| `--print-fasta-headers` |  | VgenomeParams.print_fasta_headers | none | todo: PR-F |  |
| `--ignore-duplicates` |  | VgenomeParams.ignore_duplicates | none | todo: PR-F |  |
| `--keep-files` |  | VgenomeParams.keep_files | none | todo: PR-F |  |

## vmetadata

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--target` | -t | VmetadataParams.target | none | todo: PR-F | README.md |
| `--source` |  | VmetadataParams.source | choice | todo: PR-F | README.md, docs/usage.md, docs/verification.md |
| `--filter` | -f | VmetadataParams.filter | none | todo: PR-F |  |
| `--host` |  | VmetadataParams.host | none | todo: PR-F |  |
| `--complete-only` |  | VmetadataParams.complete_only | none | todo: PR-F |  |
| `--released-after` |  | VmetadataParams.released_after | callback | todo: PR-F |  |
| `--list` | -l | VmetadataParams.list_targets | none | todo: PR-F |  |

## Short-alias collisions

- `-t`: `--target`, `--threads`
- `-r`: `--release`, `--root-name`
- `-l`: `--level`, `--list`
- `-f`: `--filter`, `--force`
