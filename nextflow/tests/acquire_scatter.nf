#!/usr/bin/env nextflow

// Data-channel front of the pipeline: ACQUIRE (metadata -> genome) feeding the
// scatter-gather dereplication. Demonstrates that the genome channel emitted by
// ACQUIRE drops straight into DEREPLICATE_SCATTER. Run with `-stub` for a quick
// wiring check, or for real with a GTDB target and a dereplicator.

nextflow.enable.dsl = 2

include { ACQUIRE             } from '../subworkflows/local/acquire'
include { DEREPLICATE_SCATTER } from '../subworkflows/local/dereplicate_scatter'
include { run_meta            } from '../subworkflows/local/run_meta'

workflow {
    ACQUIRE(channel.value(run_meta(params)))

    // Glue until Task 5: the scatter still takes bare paths.
    def ch_genome_files = ACQUIRE.out.genomes.flatMap { _meta, files -> files }
    def ch_selection    = ACQUIRE.out.selection.map { _meta, sel -> sel }
    DEREPLICATE_SCATTER(ch_genome_files, ch_selection)

    DEREPLICATE_SCATTER.out.reps
        .map { _meta, dir -> dir }
        .collectFile(name: 'merged_path.txt', storeDir: params.outdir) { dir -> "${dir}\n" }
}
