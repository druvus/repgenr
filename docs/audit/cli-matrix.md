# CLI matrix

Generated from `tests/audit/cli_matrix.yaml` by `scripts/render_cli_matrix.py`;
`tests/unit/test_cli_matrix.py` keeps both in step with the command tree.

22 commands, 197 flags (144 with a live test or an n/a reason, 53 pending).

## Global flags

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--version` |  | n/a: prints the version and exits | none | n/a: unit test test_version_flag_prints_version | docs/adding-tools.md |
| `--container` |  | container.backend | choice | todo: PR-H | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--container-engine` |  | container.engine | none | todo: PR-H | docs/containers.md |
| `--container-cache` |  | container.cache_dir | none | todo: PR-H | docs/containers.md, docs/verification.md |
| `--platform` |  | container.platform | none | todo: PR-H | README.md, docs/containers.md |
| `--wave` |  | container.wave_enabled | none | todo: PR-H | docs/adding-tools.md, docs/containers.md, docs/swot-viral.md |
| `--force` | -f | state.force | none | tests/live/test_aux_commands.py::test_second_run_skips_and_force_reruns | README.md |
| `--verbose` | -v | state.log_level | none | tests/live/test_aux_commands.py::test_logging_flags_and_env | README.md |
| `--quiet` | -q | state.log_level | none | tests/live/test_aux_commands.py::test_logging_flags_and_env |  |

## derep-stock

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/containers.md, docs/output.md |
| `--action` |  | DerepStockParams.action | choice | tests/live/test_aux_commands.py::test_derep_stock_round_trip |  |
| `--name` |  | DerepStockParams.name | none | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/usage.md |

## derep-unpack

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_aux_commands.py::test_derep_unpack_with_and_without_representant | docs/containers.md, docs/output.md |
| `--no-representant` |  | DerepUnpackParams.no_representant | none | tests/live/test_aux_commands.py::test_derep_unpack_with_and_without_representant |  |

## dereplicate

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--tool` |  | DereplicateParams.tool | registry | tests/live/test_dereplicators.py::test_adapter_recovers_the_synthetic_partition | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` | -pani | DereplicateParams.primary_ani | unit_interval | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--secondary-ani` | -sani | DereplicateParams.secondary_ani | unit_interval | tests/live/test_dereplicators.py::test_secondary_ani_sweep_changes_representative_count | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` | -af | DereplicateParams.aligned_fraction | unit_interval | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path |  |
| `--threads` | -t | DereplicateParams.threads | range | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--process-size` | -s | DereplicateParams.process_size | none | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/architecture.md, docs/scaling-audit.md, docs/swot-derep.md, docs/verification.md |
| `--num-processes` | -p | DereplicateParams.num_processes | none | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--pre-primary-ani` |  | DereplicateParams.pre_primary_ani | unit_interval | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--pre-secondary-ani` |  | DereplicateParams.pre_secondary_ani | unit_interval | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--reduce` |  | DereplicateParams.reduce | choice | tests/live/test_dereplicators.py::test_reduce_species_keeps_one_representative_per_species | docs/scaling-audit.md, docs/swot-derep.md |
| `--target-reps` |  | DereplicateParams.target_reps | range | tests/live/test_dereplicators.py::test_target_reps_lands_on_the_requested_count | docs/scaling-audit.md, docs/swot-derep.md |
| `--virus` |  | DereplicateParams.extra | none | tests/live/test_dereplicators.py::test_virus_flag_on_auto_tool_is_reported_when_ignored | README.md, docs/adding-tools.md, docs/swot-derep.md, docs/verification.md |
| `--tool-arg` |  | DereplicateParams.extra | callback | tests/live/test_dereplicators.py::test_tool_arg_reaches_the_tool_command_line | docs/adding-tools.md, docs/usage.md |
| `--allow-incomplete` |  | DereplicateParams.allow_incomplete | none | tests/live/test_dereplicators.py::test_allow_incomplete_gates_a_missing_genome |  |
| `--keeper` |  | DereplicateParams.keeper | choice | tests/live/test_dereplicators.py::test_keeper_quality_promotes_the_best_scored_member | README.md, docs/swot-derep.md, docs/usage.md |

## dereplicate-chunk

dispatch: `step:repgenr.stages.derep_steps.dereplicate_chunk`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--genomes-fofn` |  | ChunkParams.genomes | stage | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions |  |
| `--out` | -o | ChunkParams.out_dir | none | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions |  |
| `--tool` |  | ChunkParams.tool | registry | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` | -pani | ChunkParams.primary_ani | unit_interval | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/verification.md |
| `--secondary-ani` | -sani | ChunkParams.secondary_ani | unit_interval | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` | -af | ChunkParams.aligned_fraction | unit_interval | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) |  |
| `--threads` | -t | ChunkParams.threads | range | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/verification.md |
| `--virus` |  | ChunkParams.extra | none | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | README.md, docs/adding-tools.md, docs/swot-derep.md, docs/verification.md |
| `--tool-arg` |  | ChunkParams.extra | callback | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/adding-tools.md, docs/usage.md |
| `--selection-tsv` |  | ChunkParams.selection_tsv | none | todo: PR-G |  |
| `--keeper` |  | ChunkParams.keeper | choice | todo: PR-G | README.md, docs/swot-derep.md, docs/usage.md |
| `--versions-out` |  | ChunkParams.versions_out | none | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions |  |

## dereplicate-merge

dispatch: `step:repgenr.stages.derep_steps.dereplicate_merge`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--out` | -o | MergeParams.out_dir | none | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition |  |
| `--chunk-dir` |  | MergeParams.chunk_dirs | none | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition |  |
| `--chunk-fofn` |  | MergeParams.chunk_dirs | stage | tests/live/test_steps.py::test_merge_by_chunk_fofn |  |
| `--tool` |  | MergeParams.tool | registry | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` | -pani | MergeParams.primary_ani | unit_interval | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/verification.md |
| `--secondary-ani` | -sani | MergeParams.secondary_ani | unit_interval | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` | -af | MergeParams.aligned_fraction | unit_interval | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) |  |
| `--threads` | -t | MergeParams.threads | range | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | docs/verification.md |
| `--virus` |  | MergeParams.extra | none | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | README.md, docs/adding-tools.md, docs/swot-derep.md, docs/verification.md |
| `--tool-arg` |  | MergeParams.extra | callback | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/adding-tools.md, docs/usage.md |
| `--selection-tsv` |  | MergeParams.selection_tsv | none | todo: PR-G |  |
| `--keeper` |  | MergeParams.keeper | choice | todo: PR-G | README.md, docs/swot-derep.md, docs/usage.md |
| `--versions-out` |  | MergeParams.versions_out | none | tests/live/test_steps.py::test_merge_by_chunk_fofn |  |

## doctor

dispatch: `query`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | tests/live/test_aux_commands.py::test_doctor_passes_then_fails_on_a_corrupt_genome | docs/containers.md, docs/output.md |

## genome

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_network.py::test_api_genus_representatives | docs/containers.md, docs/output.md |
| `--accession-list-only` |  | GenomeParams.accession_list_only | none | tests/live/test_network.py::test_genome_accession_list_only_is_a_pure_query |  |
| `--keep-files` |  | GenomeParams.keep_files | none | tests/live/test_network.py::test_genome_keep_files_retains_the_download_scratch |  |

## genome-fetch

dispatch: `step:repgenr.stages.genome_steps.genome_fetch`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--selection` |  | GenomeFetchParams.selection_tsv | stage | todo: PR-F | README.md, docs/usage.md |
| `--out` | -o | GenomeFetchParams.out_dir | none | todo: PR-F |  |
| `--keep-files` |  | GenomeFetchParams.keep_files | none | todo: PR-F |  |
| `--versions-out` |  | GenomeFetchParams.versions_out | none | todo: PR-G |  |

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
| `--selection` |  | IngestParams.selection | stage | tests/live/test_ingest_flags.py::test_selection_table_drives_taxonomy_and_subset | README.md, docs/usage.md |
| `--outgroup` |  | IngestParams.outgroup | stage | tests/live/test_ingest_flags.py::test_outgroup_and_copy | README.md, docs/usage.md |
| `--copy` |  | IngestParams.copy | none | tests/live/test_ingest_flags.py::test_outgroup_and_copy | README.md, docs/usage.md |

## list-tools

dispatch: `query`

(no flags)

## metadata

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_network.py::test_api_genus_representatives | docs/containers.md, docs/output.md |
| `--dataset` | -d | MetadataParams.dataset | choice | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup |  |
| `--level` | -l | MetadataParams.level | choice | tests/live/test_network.py::test_api_target_family_widens_the_selection |  |
| `--source` |  | MetadataParams.source | choice | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md, docs/usage.md, docs/verification.md |
| `--release` | -r | MetadataParams.release | none | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md |
| `--gtdb-version` |  | MetadataParams.version | none | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md, docs/usage.md |
| `--target-family` | -tf | MetadataParams.target_family | none | tests/live/test_network.py::test_api_target_family_widens_the_selection |  |
| `--target-genus` | -tg | MetadataParams.target_genus | none | tests/live/test_network.py::test_api_genus_representatives | README.md |
| `--target-species` | -ts | MetadataParams.target_species | none | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup |  |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup |  |
| `--metadata-path` |  | MetadataParams.metadata_path | none | tests/live/test_network.py::test_tsv_nodownload_and_metadata_path_reuse_the_table |  |
| `--nodownload` |  | MetadataParams.nodownload | none | tests/live/test_network.py::test_tsv_nodownload_and_metadata_path_reuse_the_table |  |
| `--limit` |  | MetadataParams.limit | range | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/scaling-audit.md, docs/swot-derep.md, docs/usage.md |

## phylo

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | tests/live/test_treebuilders_offline.py::test_alignment_free_builder_on_representatives | README.md, docs/containers.md, docs/usage.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | todo: PR-G | README.md, docs/usage.md |
| `--aligner` |  | PhyloParams.aligner | registry | todo: PR-H | README.md, docs/containers.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | todo: PR-G |  |
| `--all-genomes` |  | PhyloParams.all_genomes | none | tests/live/test_treebuilders_offline.py::test_all_genomes_puts_every_genome_in_the_tree |  |
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
| `--genomes-dir` |  | PhyloBuildParams.genomes_dir | stage | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/usage.md |
| `--out` | -o | PhyloBuildParams.out_dir | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--outgroup-dir` |  | PhyloBuildParams.outgroup_dir | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--outgroup-accession` |  | PhyloBuildParams.outgroup_accession | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--treebuilder` |  | PhyloBuildParams.phylo.treebuilder | registry | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/containers.md, docs/usage.md |
| `--msa-source` |  | PhyloBuildParams.phylo.msa_source | choice | todo: PR-G | README.md, docs/usage.md |
| `--aligner` |  | PhyloBuildParams.phylo.aligner | registry | todo: PR-H | README.md, docs/containers.md |
| `--snptyper` |  | PhyloBuildParams.phylo.snptyper | registry | todo: PR-G |  |
| `--no-outgroup` |  | PhyloBuildParams.phylo.no_outgroup | none | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree |  |
| `--bootstrap` | -B | PhyloBuildParams.phylo.bootstrap | range | todo: PR-G | docs/usage.md |
| `--reference` |  | PhyloBuildParams.phylo.reference | none | todo: PR-G | docs/swot-phylo.md |
| `--aligner-arg` |  | PhyloBuildParams.phylo.extra | callback | todo: PR-H | docs/adding-tools.md |
| `--threads` | -t | PhyloBuildParams.phylo.threads | range | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/verification.md |
| `--versions-out` |  | PhyloBuildParams.versions_out | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |

## run

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/containers.md, docs/output.md |
| `--viral` |  | VmetadataParams | none | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md |
| `--dataset` | -d | MetadataParams.dataset | choice | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--level` | -l | MetadataParams.level | choice | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--target-family` | -tf | MetadataParams.target_family | none | todo: PR-F |  |
| `--target-genus` | -tg | MetadataParams.target_genus | none | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md |
| `--target-species` | -ts | MetadataParams.target_species | none | todo: PR-F |  |
| `--release` | -r | MetadataParams.release | none | todo: PR-F | README.md |
| `--gtdb-version` |  | MetadataParams.version | none | todo: PR-F | README.md, docs/usage.md |
| `--metadata-source` |  | MetadataParams.source | choice | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--target` | -t | VmetadataParams.target | none | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md |
| `--viral-source` |  | VmetadataParams.source | choice | tests/live/test_network.py::test_run_viral_chain_end_to_end |  |
| `--group-segments` |  | VgenomeParams.group_segments | none | todo: PR-F | README.md |
| `--tool` |  | DereplicateParams.tool | registry | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` |  | DereplicateParams.primary_ani | unit_interval | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/verification.md |
| `--secondary-ani` |  | DereplicateParams.secondary_ani | unit_interval | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` |  | DereplicateParams.aligned_fraction | unit_interval | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--keeper` |  | DereplicateParams.keeper | choice | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/swot-derep.md, docs/usage.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/containers.md, docs/usage.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | todo: PR-F | README.md, docs/usage.md |
| `--aligner` |  | PhyloParams.aligner | registry | todo: PR-F | README.md, docs/containers.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | todo: PR-F |  |
| `--no-outgroup` |  | PhyloParams.no_outgroup | none | tests/live/test_network.py::test_run_viral_chain_end_to_end |  |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md, docs/usage.md |
| `--threads` |  | DereplicateParams.threads | range | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/verification.md |
| `--dry-run` |  | n/a: prints the chain and exits before any stage | none | tests/live/test_network.py::test_run_dry_run_prints_the_chain_without_network |  |

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
| `--node-basename` |  | Tree2taxParams.node_basename | none | todo: PR-G |  |
| `--root-name` | -r | Tree2taxParams.root_name | none | todo: PR-G |  |
| `--remove-outgroup` |  | Tree2taxParams.remove_outgroup | none | todo: PR-G |  |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | README.md, docs/usage.md |
| `--collapse-support` |  | Tree2taxParams.collapse_support | range | todo: PR-G | docs/swot-phylo.md, docs/usage.md |
| `--collapse-length` |  | Tree2taxParams.collapse_length | range | todo: PR-G | docs/swot-phylo.md, docs/usage.md |

## tree2tax-relations

dispatch: `step:repgenr.stages.tree2tax.tree2tax_relations`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--tree` |  | Tree2taxStepParams.tree | stage | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--out` | -o | Tree2taxStepParams.out_dir | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--clusters` |  | Tree2taxStepParams.clusters | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--outgroup-dir` |  | Tree2taxStepParams.outgroup_dir | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--outgroup-accession` |  | Tree2taxStepParams.outgroup_accession | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--node-basename` |  | Tree2taxStepParams.node_basename | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--root-name` | -r | Tree2taxStepParams.root_name | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--remove-outgroup` |  | Tree2taxStepParams.remove_outgroup | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--include-dereplicated` |  | Tree2taxStepParams.include_dereplicated | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/usage.md |
| `--versions-out` |  | Tree2taxStepParams.versions_out | none | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--collapse-support` |  | Tree2taxStepParams.collapse_support | range | tests/live/test_steps.py::test_tree2tax_relations_collapse_flags | docs/swot-phylo.md, docs/usage.md |
| `--collapse-length` |  | Tree2taxStepParams.collapse_length | range | tests/live/test_steps.py::test_tree2tax_relations_collapse_flags | docs/swot-phylo.md, docs/usage.md |

## versions

dispatch: `query`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | tests/live/test_aux_commands.py::test_status_and_versions | docs/containers.md, docs/output.md |
| `--versions-out` |  | n/a: query command, no stage parameters | none | tests/live/test_aux_commands.py::test_status_and_versions |  |

## vgenome

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_network.py::test_vgenome_selection_flags | docs/containers.md, docs/output.md |
| `--target-genus` | -tg | VgenomeParams.target_genus | none | tests/live/test_network.py::test_vgenome_selection_flags | README.md |
| `--target-species` | -ts | VgenomeParams.target_species | none | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--target-serotype` | -tse | VgenomeParams.target_serotype | none | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--target-custom` | -tc | VgenomeParams.target_custom | none | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--length-all` |  | VgenomeParams.length_all | none | tests/live/test_network.py::test_vgenome_selection_flags | docs/usage.md |
| `--length-deviation` |  | VgenomeParams.length_deviation | range | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--length-method` |  | VgenomeParams.length_method | choice | tests/live/test_network.py::test_vgenome_selection_flags | docs/scaling-audit.md, docs/swot-viral.md, docs/usage.md |
| `--length-range` |  | VgenomeParams.length_range | none | tests/live/test_network.py::test_vgenome_selection_flags | docs/usage.md |
| `--discard` |  | VgenomeParams.discard | none | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |
| `--no-outgroup` |  | VgenomeParams.no_outgroup | none | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--group-segments` |  | VgenomeParams.group_segments | none | tests/live/test_network.py::test_vgenome_selection_flags | README.md |
| `--outgroup-candidates-taxid-min-genomes` |  | VgenomeParams.outgroup_candidates_taxid_min_genomes | none | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--outgroup-treebuilder` |  | VgenomeParams.outgroup_treebuilder | registry | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |
| `--glance` |  | VgenomeParams.glance | none | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |
| `--print-fasta-headers` |  | VgenomeParams.print_fasta_headers | none | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |
| `--ignore-duplicates` |  | VgenomeParams.ignore_duplicates | none | tests/live/test_network.py::test_vgenome_bvbrc_needs_ignore_duplicates |  |
| `--keep-files` |  | VgenomeParams.keep_files | none | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |

## vmetadata

dispatch: `stage`

| flag | aliases | param | validated | live | docs |
|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | docs/containers.md, docs/output.md |
| `--target` | -t | VmetadataParams.target | none | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | README.md |
| `--source` |  | VmetadataParams.source | choice | tests/live/test_network.py::test_vmetadata_bvbrc_source_and_filter | README.md, docs/usage.md, docs/verification.md |
| `--filter` | -f | VmetadataParams.filter | none | tests/live/test_network.py::test_vmetadata_bvbrc_source_and_filter |  |
| `--host` |  | VmetadataParams.host | none | tests/live/test_network.py::test_vmetadata_released_after_and_host_narrow_the_set |  |
| `--complete-only` |  | VmetadataParams.complete_only | none | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only |  |
| `--released-after` |  | VmetadataParams.released_after | callback | tests/live/test_network.py::test_vmetadata_released_after_and_host_narrow_the_set |  |
| `--list` | -l | VmetadataParams.list_targets | none | tests/live/test_network.py::test_vmetadata_list_targets_reaches_bvbrc |  |

## Short-alias collisions

- `-t`: `--target`, `--threads`
- `-r`: `--release`, `--root-name`
- `-l`: `--level`, `--list`
- `-f`: `--filter`, `--force`
