// Dereplicate one chunk of genomes (scatter step).
//
// Data-channel module: genome FASTAs are staged in as channel inputs (not read
// from a shared workdir), and the chunk result directory is emitted as a typed
// output other processes consume. Wraps `repgenr dereplicate-chunk`. Tool
// flags arrive as task.ext.args from conf/modules.config.

process DEREP_CHUNK {
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(genomes, stageAs: 'inputs/*')
    path selection, stageAs: 'selection.tsv'

    output:
    tuple val(meta), path("${meta.id}"), emit: chunk
    path 'versions.yml'                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    # Build a file-of-filenames from the staged genomes (never argv -- ARG_MAX).
    ls -1 inputs/* > genomes.fofn

    # selection.tsv is optional: no bacterial ACQUIRE selection (viral path, or
    # a harness with no metadata front) stages nothing under that name, so the
    # quality-aware keeper is skipped rather than passed a missing file. Every
    # genome in the chunk has a real file here (unlike at the merge step), so
    # any promotion the keeper makes is always resolvable.
    sel=""
    [ -e selection.tsv ] && sel="--selection-tsv selection.tsv"

    repgenr ${opts} dereplicate-chunk \\
        --genomes-fofn genomes.fofn \\
        --out ${meta.id} \\
        ${args} \\
        \$sel \\
        --threads ${task.cpus} \\
        --versions-out tool_versions.yml

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    repgenr: \$(repgenr --version | sed 's/repgenr //')
END_VERSIONS
    cat tool_versions.yml >> versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    mkdir -p ${meta.id}/representatives
    printf 'representative\\tmember\\n' > ${meta.id}/clusters.tsv
    printf 'genome\\tstatus\\n' > ${meta.id}/genome_status.tsv
    for f in inputs/*; do
        b=\$(basename \$f)
        cp \$f ${meta.id}/representatives/\$b
        printf '%s\\t%s\\n' "\$b" "\$b" >> ${meta.id}/clusters.tsv
        printf '%s\\trepresentative\\n' "\$b" >> ${meta.id}/genome_status.tsv
    done
    touch versions.yml
    """
}
