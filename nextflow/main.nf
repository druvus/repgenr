#!/usr/bin/env nextflow

// RepGenR pipeline entry point.
//
// RepGenR manages a single working directory (params.workdir) that every stage
// reads and writes. The Nextflow layer orders the stages and assigns resource
// labels/profiles for HPC and cloud execution; it does not stage genome files
// between processes, since the working directory is the shared state.

nextflow.enable.dsl = 2

include { BACTERIAL } from './subworkflows/local/bacterial'

workflow {
    if (!params.workdir) {
        error "Provide --workdir <path> (the shared RepGenR working directory)."
    }
    BACTERIAL()
}
