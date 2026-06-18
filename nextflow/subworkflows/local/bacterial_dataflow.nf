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
    main:
    ch_versions = Channel.empty()

    ACQUIRE()
    ch_versions = ch_versions.mix(ACQUIRE.out.versions)

    DEREPLICATE_SCATTER(ACQUIRE.out.genomes)
    ch_versions = ch_versions.mix(DEREPLICATE_SCATTER.out.versions)

    ch_reps      = DEREPLICATE_SCATTER.out.reps.map { meta, dir -> dir }
    ch_outgroup  = ACQUIRE.out.outgroup.collect().ifEmpty([])
    ch_selection = ACQUIRE.out.selection

    PHYLO(ch_reps, ch_outgroup, ch_selection)
    ch_versions = ch_versions.mix(PHYLO.out.versions)

    TREE2TAX(PHYLO.out.tree, ch_reps, ch_outgroup, ch_selection)
    ch_versions = ch_versions.mix(TREE2TAX.out.versions)

    emit:
    tree        = PHYLO.out.tree
    tree2tax    = TREE2TAX.out.tree2tax
    genomes_map = TREE2TAX.out.genomes_map
    versions    = ch_versions
}
