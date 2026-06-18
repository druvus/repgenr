#!/usr/bin/env nextflow

// RepGenR pipeline entry point.
//
// RepGenR manages a single working directory (params.workdir) that every stage
// reads and writes. The Nextflow layer orders the stages and assigns resource
// labels/profiles for HPC and cloud execution; it does not stage genome files
// between processes, since the working directory is the shared state.
//
// (The data-channel rewrite that replaces this shared-workdir contract with
// typed tuple(meta, path) channels lands in later Phase 4 increments.)

nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { BACTERIAL } from './subworkflows/local/bacterial'
include { VIRAL } from './subworkflows/local/viral'

workflow {
    // Print usage and exit before validation, so `--help` needs no other params.
    if (params.help) {
        log.info """
        RepGenR -- ${workflow.manifest.description} (v${workflow.manifest.version})

        Usage:
          nextflow run main.nf --workdir <DIR> [-profile standard|slurm|cloud|docker|singularity|wave|test]

        Key parameters (see nextflow_schema.json and docs/usage.md for all):
          --workdir <DIR>        Shared RepGenR working directory (required)
          --outdir <DIR>         Reports / published results (default: results)
          --mode bacterial|viral Lineage pipeline to run (default: bacterial)
          --*_args '<str>'       Per-stage repgenr arguments (metadata_args, etc.)
          --derep_process_size N Two-stage dereplication chunk size (large inputs)
        """.stripIndent()
        return
    }

    // Validate --params against nextflow_schema.json and echo the resolved set.
    validateParameters()
    log.info paramsSummaryLog(workflow)

    if (!params.workdir) {
        error "Provide --workdir <path> (the shared RepGenR working directory)."
    }
    if (params.mode == 'viral') {
        VIRAL()
    }
    else {
        BACTERIAL()
    }
}
