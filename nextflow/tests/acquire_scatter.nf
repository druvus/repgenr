#!/usr/bin/env nextflow

// Data-channel front of the pipeline: ACQUIRE (metadata -> genome) feeding the
// scatter-gather dereplication. Demonstrates that the genome channel emitted by
// ACQUIRE drops straight into DEREPLICATE_SCATTER. Run with `-stub` for a quick
// wiring check, or for real with a GTDB target and a dereplicator.

nextflow.enable.dsl = 2

include { ACQUIRE             } from '../subworkflows/local/acquire'
include { DEREPLICATE_SCATTER } from '../subworkflows/local/dereplicate_scatter'

workflow {
    ACQUIRE()
    DEREPLICATE_SCATTER(ACQUIRE.out.genomes)

    DEREPLICATE_SCATTER.out.reps
        .map { meta, dir -> dir }
        .collectFile(name: 'merged_path.txt', storeDir: params.outdir) { dir -> "${dir}\n" }
}
