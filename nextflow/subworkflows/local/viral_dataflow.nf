// Viral pipeline, data-channel form:
//   VACQUIRE (vmetadata -> vgenome) -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
//
// Reuses the same data-channel dereplication, phylo and tree2tax as the bacterial
// path; only the acquisition front differs (BV-BRC instead of GTDB/NCBI).

include { VACQUIRE            } from '../../modules/local/dataflow/vacquire'
include { DEREPLICATE_SCATTER } from './dereplicate_scatter'
include { PHYLO               } from '../../modules/local/dataflow/phylo'
include { TREE2TAX            } from '../../modules/local/dataflow/tree2tax'

workflow VIRAL_DATAFLOW {
    take:
    ch_meta

    main:
    def ch_versions = channel.empty()

    VACQUIRE(ch_meta)
    ch_versions = ch_versions.mix(VACQUIRE.out.versions)

    def ch_genomes = VACQUIRE.out.genomes
        .map { meta, files -> tuple(meta, files instanceof List ? files : [files]) }
    // vmetadata writes no selection.tsv and viral genomes carry no CheckM
    // quality, so the keeper input is an empty list under the run meta.
    def ch_no_selection = ch_meta.map { meta -> tuple(meta, []) }
    DEREPLICATE_SCATTER(ch_genomes, ch_no_selection)
    ch_versions = ch_versions.mix(DEREPLICATE_SCATTER.out.versions)

    // Glue until Task 6: phylo and tree2tax still take bare paths.
    def ch_reps     = DEREPLICATE_SCATTER.out.reps.map { _meta, dir -> dir }
    def ch_outgroup = VACQUIRE.out.outgroup.map { _meta, files -> files }.ifEmpty([])
    def ch_og_acc   = VACQUIRE.out.outgroup_accession.map { _meta, acc -> acc }

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
