// Bacterial pipeline, data-channel form:
//   ACQUIRE (metadata -> genome) -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
//
// Genomes, representatives, the tree and the taxonomy flow between processes as
// staged channel files; there is no shared working directory. This is the
// data-channel replacement for the legacy done-signal BACTERIAL subworkflow.

include { ACQUIRE             } from './acquire'
include { DEREPLICATE_SCATTER } from './dereplicate_scatter'
include { PHYLO               } from '../../modules/local/dataflow/phylo'
include { TREE2TAX            } from '../../modules/local/dataflow/tree2tax'

workflow BACTERIAL_DATAFLOW {
    take:
    ch_meta

    main:
    def ch_versions = channel.empty()

    ACQUIRE(ch_meta)
    ch_versions = ch_versions.mix(ACQUIRE.out.versions)

    DEREPLICATE_SCATTER(ACQUIRE.out.genomes, ACQUIRE.out.selection)
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
