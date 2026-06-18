#!/usr/bin/env nextflow

// Standalone harness for the scatter-gather dereplication subworkflow.
//
// Feeds a directory of genome FASTAs straight into DEREPLICATE_SCATTER, with no
// GTDB/NCBI front-end, so the scatter-gather mechanism can be exercised on local
// files (real tools or `-stub`). The merged representative set is published to
// `--outdir`.
//
//   nextflow run nextflow/tests/dereplicate_scatter.nf \
//       --genomes_dir <DIR> --derep_tool sourmash -profile test

nextflow.enable.dsl = 2

include { DEREPLICATE_SCATTER } from '../subworkflows/local/dereplicate_scatter'

params.genomes_dir = null

workflow {
    if (!params.genomes_dir) {
        error "Provide --genomes_dir <DIR> with genome FASTAs."
    }

    ch_genomes = Channel
        .fromPath("${params.genomes_dir}/*.{fasta,fa,fna,fas}")
        .filter { !it.name.startsWith('._') }   // skip macOS AppleDouble files

    DEREPLICATE_SCATTER(ch_genomes)

    DEREPLICATE_SCATTER.out.reps
        .map { meta, dir -> dir }
        .collectFile(name: 'merged_path.txt', storeDir: params.outdir) { dir ->
            "${dir}\n"
        }
}
