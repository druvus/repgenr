#!/usr/bin/env nextflow

// Standalone harness for the full viral data-channel pipeline
// (VACQUIRE -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX). Run with `-stub` for a
// wiring check, or for real with a BV-BRC target and a dereplicator/tree builder.

nextflow.enable.dsl = 2

include { VIRAL_DATAFLOW } from '../subworkflows/local/viral_dataflow'

workflow {
    VIRAL_DATAFLOW()
}
