#!/usr/bin/env nextflow

// RepGenR pipeline entry point (data-channel form).
//
// Stages exchange typed channel files -- metadata selection, genome FASTAs,
// per-chunk and merged representatives, the tree and the taxonomy -- rather than
// sharing one working directory. Nextflow owns the fan-out (scatter-gather
// dereplication), and results are published under --outdir.

nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { BACTERIAL_DATAFLOW } from './subworkflows/local/bacterial_dataflow'
include { VIRAL_DATAFLOW     } from './subworkflows/local/viral_dataflow'

workflow {
    // Print usage and exit before validation, so `--help` needs no other params.
    if (params.help) {
        log.info """
        RepGenR -- ${workflow.manifest.description} (v${workflow.manifest.version})

        Usage:
          nextflow run main.nf --outdir <DIR> [-profile standard|slurm|cloud|docker|singularity|wave|test]

        Key parameters (see nextflow_schema.json and docs/usage.md for all):
          --outdir <DIR>         Published results (default: results)
          --mode bacterial|viral Lineage pipeline to run (default: bacterial)
          --metadata_args '<str>'  GTDB selection (bacterial)
          --vmetadata_args / --vgenome_args '<str>'  BV-BRC selection (viral)
          --derep_tool <tool>    Dereplicator for the scatter-gather step
          --derep_process_size N Genomes per dereplication chunk (large inputs)
          --phylo_args '<str>'   Aligner or tree builder for the phylogeny
        """.stripIndent()
        return
    }

    validateParameters()
    log.info paramsSummaryLog(workflow)

    if (params.mode == 'viral') {
        VIRAL_DATAFLOW()
    }
    else {
        BACTERIAL_DATAFLOW()
    }
}
