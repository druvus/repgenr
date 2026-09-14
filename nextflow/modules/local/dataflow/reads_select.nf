// Select sequencing runs from ENA/SRA (data-channel form).
//
// Runs the reads stage in a task-local working directory and emits the
// portable reads.tsv (run, sample, organism, platform, layout, FASTQ URLs) as a
// channel output -- the hand-off the per-run assembly step consumes. No shared
// workdir. Tool flags arrive as task.ext.args from conf/modules.config;
// publishing is configured there too.

process READS_SELECT {
    label 'process_low'
    tag "${meta.id}"

    input:
    val meta

    output:
    tuple val(meta), path("reads.tsv"), emit: reads
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} reads -wd reads_wd ${args}
    cp reads_wd/reads.tsv reads.tsv

    repgenr versions -wd reads_wd --versions-out tool_versions.yml
    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    printf 'run_accession\\tbiosample\\tbioproject\\torganism\\ttaxid\\tplatform\\tinstrument_model\\tlayout\\tbases\\tread_count\\tfamily\\tgenus\\tspecies\\tfastq_urls\\tfastq_md5\\tfastq_bytes\\n' > reads.tsv
    printf 'SRR000001\\tSAMN000001\\tPRJNA000001\\tFrancisella tularensis\\t263\\tILLUMINA\\tIllumina MiSeq\\tPAIRED\\t600000000\\t2000000\\tFrancisellaceae\\tFrancisella\\ttularensis\\thttps://example.org/SRR000001_1.fastq.gz;https://example.org/SRR000001_2.fastq.gz\\ta;b\\t100;100\\n' >> reads.tsv
    printf 'SRR000002\\tSAMN000002\\tPRJNA000001\\tFrancisella tularensis\\t263\\tOXFORD_NANOPORE\\tGridION\\tSINGLE\\t400000000\\t40000\\tFrancisellaceae\\tFrancisella\\ttularensis\\thttps://example.org/SRR000002_1.fastq.gz\\tc\\t100\\n' >> reads.tsv
    printf 'SRR000003\\tSAMN000003\\tPRJNA000001\\tFrancisella tularensis\\t263\\tION_TORRENT\\tIon Torrent S5\\tSINGLE\\t100000000\\t500000\\tFrancisellaceae\\tFrancisella\\ttularensis\\thttps://example.org/SRR000003.fastq.gz\\td\\t100\\n' >> reads.tsv
    touch versions.yml
    """
}
