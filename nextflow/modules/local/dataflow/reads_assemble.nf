// Fetch and assemble one sequencing run (data-channel form).
//
// One task per run: `repgenr assemble-run` downloads the run's FASTQ files,
// assembles them with the adapter that accepts the platform, and writes a
// directory named by the run accession holding contigs.fasta and the
// assembly.ok marker -- or excused_runs.tsv when the run cannot be fetched or
// assembled. The directory is emitted as a typed output for the QC and gather
// steps. Memory follows the platform (conf/base.config, process_assembly,
// through ext.platform); the meta carries the run accession and platform.

process READS_ASSEMBLE {
    label 'process_assembly'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads_tsv)

    output:
    tuple val(meta), path("${meta.accession}"), emit: assembly
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} assemble-run \\
        --reads-tsv ${reads_tsv} \\
        --run ${meta.accession} \\
        --out ${meta.accession} \\
        ${args} \\
        --threads ${task.cpus} \\
        --memory-gb ${task.memory.toGiga()} \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    echo "ext.platform: ${task.ext.platform ?: ''}"
    mkdir -p ${meta.accession}
    if [ "${meta.platform}" = "ILLUMINA" ] || [ "${meta.platform}" = "OXFORD_NANOPORE" ] || [ "${meta.platform}" = "PACBIO_SMRT" ]; then
        printf '>${meta.accession}_contig1\\nACGTACGTACGT\\n' > ${meta.accession}/contigs.fasta
        printf '{"assembler": "stub", "version": "stub", "stats": {"n_contigs": 1, "total_length": 12, "n50": 12, "largest_contig": 12}, "tool_stats": {}}\\n' > ${meta.accession}/assembly.ok
    else
        printf 'run_accession\\tstep\\treason\\n${meta.accession}\\tassemble\\tunsupported_platform\\n' > ${meta.accession}/excused_runs.tsv
    fi
    touch versions.yml
    """
}
