// Scatter-gather dereplication (data-channel style).
//
// Genomes arrive as a channel of individual FASTA paths. They are grouped into
// chunks of `derep_process_size`, each chunk is dereplicated independently and
// in parallel (scatter, one task per chunk -- on HPC these land on separate
// nodes), and the union of the chunk representatives is dereplicated once more
// (gather) to produce the final representative set. This is the two-stage
// reduce-tree of the in-process dereplicate stage, expressed as typed channels
// so Nextflow owns the fan-out instead of Python threads.
//
// A single round (chunk -> merge) is run here. For very large inputs whose
// merged union still exceeds a chunk, feed `reps` back through another round.

include { DEREP_CHUNK } from '../../modules/local/derep_chunk'
include { DEREP_MERGE } from '../../modules/local/derep_merge'

workflow DEREPLICATE_SCATTER {
    take:
    ch_genomes   // channel: individual genome FASTA paths

    main:
    ch_versions = Channel.empty()

    // Group genomes into chunks; a null/zero process size means a single chunk.
    // Coerce because command-line params arrive as strings.
    def requested_size = params.derep_process_size ? (params.derep_process_size as Integer) : 0
    def chunk_size = requested_size > 0 ? requested_size : 1000000

    ch_chunks = ch_genomes
        .collect()
        .flatMap { files ->
            files.collate(chunk_size).withIndex().collect { chunk, i ->
                tuple([id: "chunk_${i}"], chunk)
            }
        }

    DEREP_CHUNK(ch_chunks)
    ch_versions = ch_versions.mix(DEREP_CHUNK.out.versions.first())

    // Gather every chunk result directory into a single merge input.
    ch_merge_in = DEREP_CHUNK.out.chunk
        .map { meta, dir -> dir }
        .collect()
        .map { dirs -> tuple([id: 'merged'], dirs) }

    DEREP_MERGE(ch_merge_in)
    ch_versions = ch_versions.mix(DEREP_MERGE.out.versions)

    emit:
    reps     = DEREP_MERGE.out.reps   // tuple(meta, path to merged result dir)
    versions = ch_versions
}
