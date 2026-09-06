// Dereplicate the union of chunk representatives (gather step).
//
// Data-channel module: every chunk result directory produced by DEREP_CHUNK is
// staged in, and the merged representative set is emitted as a typed output.
// Wraps `repgenr dereplicate-merge`. Tool flags arrive as task.ext.args from
// conf/modules.config; publishing is configured there too.

process DEREP_MERGE {
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(chunks, stageAs: 'chunks/*')
    path selection, stageAs: 'selection.tsv'

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}"), emit: reps
    path 'versions.yml'                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    # One --chunk-dir per staged chunk directory.
    args=""
    for d in chunks/*; do
        args="\$args --chunk-dir \$d"
    done

    # selection.tsv is optional: no bacterial ACQUIRE selection (viral path, or
    # a harness with no metadata front) stages nothing under that name, so the
    # quality-aware keeper is skipped rather than passed a missing file.
    sel=""
    [ -e selection.tsv ] && sel="--selection-tsv selection.tsv"

    repgenr ${opts} dereplicate-merge \\
        \$args \\
        --out ${prefix} \\
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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "ext.args: ${args}"
    mkdir -p ${prefix}/representatives
    printf 'representative\\tmember\\n' > ${prefix}/clusters.tsv
    printf 'genome\\tstatus\\n' > ${prefix}/genome_status.tsv
    for d in chunks/*; do
        for f in \$d/representatives/*; do
            [ -e "\$f" ] || continue
            b=\$(basename \$f)
            cp \$f ${prefix}/representatives/\$b
            printf '%s\\t%s\\n' "\$b" "\$b" >> ${prefix}/clusters.tsv
            printf '%s\\trepresentative\\n' "\$b" >> ${prefix}/genome_status.tsv
        done
    done
    touch versions.yml
    """
}
