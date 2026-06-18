// Acquire genomes (data-channel form): metadata -> genome.
//
// METADATA selects accessions and emits selection.tsv; GENOME downloads them and
// emits the genome FASTAs. The genomes are flattened to individual paths so the
// output drops straight into DEREPLICATE_SCATTER, which takes a channel of
// individual genome files.

include { METADATA } from '../../modules/local/dataflow/metadata'
include { GENOME   } from '../../modules/local/dataflow/genome'

workflow ACQUIRE {
    main:
    ch_versions = Channel.empty()

    METADATA()
    ch_versions = ch_versions.mix(METADATA.out.versions)

    GENOME(METADATA.out.selection)
    ch_versions = ch_versions.mix(GENOME.out.versions)

    emit:
    genomes   = GENOME.out.genomes.flatten()   // individual genome FASTA paths
    outgroup  = GENOME.out.outgroup.flatten()
    selection = METADATA.out.selection
    versions  = ch_versions
}
