#!/usr/bin/env nextflow

// Standalone harness for the reads data-channel pipeline
// (ACQUIRE_READS -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX). Run with `-stub`
// for a wiring check, or for real with an ENA selection and the assemblers.

include { READS_DATAFLOW   } from '../subworkflows/local/reads_dataflow'
include { PUBLISH_VERSIONS } from '../subworkflows/local/publish_versions'
include { run_meta         } from '../subworkflows/local/run_meta'

workflow {
    READS_DATAFLOW(channel.value(run_meta(params)))
    PUBLISH_VERSIONS(READS_DATAFLOW.out.versions)
}
