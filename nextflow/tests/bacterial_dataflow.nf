#!/usr/bin/env nextflow

// Standalone harness for the full bacterial data-channel pipeline
// (ACQUIRE -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX). Run with `-stub` for a
// wiring check, or for real with a GTDB target and a dereplicator/tree builder.

nextflow.enable.dsl = 2

include { BACTERIAL_DATAFLOW } from '../subworkflows/local/bacterial_dataflow'

workflow {
    BACTERIAL_DATAFLOW()
}
