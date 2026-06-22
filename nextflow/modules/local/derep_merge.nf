// Dereplicate the union of chunk representatives (gather step).
//
// Data-channel module: every chunk result directory produced by DEREP_CHUNK is
// staged in, and the merged representative set is emitted as a typed output.
// Wraps `repgenr dereplicate-merge`.

process DEREP_MERGE {
    label 'process_high'
    tag "${meta.id}"
    publishDir "${params.outdir}/dereplicate", mode: 'copy'

    input:
    tuple val(meta), path(chunks, stageAs: 'chunks/*')

    output:
    tuple val(meta), path("${meta.id}"), emit: reps
    path 'versions.yml'                , emit: versions

    script:
    """
    # One --chunk-dir per staged chunk directory.
    args=""
    for d in chunks/*; do
        args="\$args --chunk-dir \$d"
    done

    repgenr ${params.repgenr_opts} dereplicate-merge \\
        \$args \\
        --out ${meta.id} \\
        --tool ${params.derep_tool} \\
        --primary-ani ${params.derep_primary_ani} \\
        --secondary-ani ${params.derep_secondary_ani} \\
        --aligned-fraction ${params.derep_aligned_fraction} \\
        --threads ${task.cpus} \\
        --versions-out tool_versions.yml

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    repgenr: \$(repgenr --version | sed 's/repgenr //')
END_VERSIONS
    cat tool_versions.yml >> versions.yml
    """

    stub:
    """
    mkdir -p ${meta.id}/representatives
    printf 'representative\\tmember\\n' > ${meta.id}/clusters.tsv
    printf 'genome\\tstatus\\n' > ${meta.id}/genome_status.tsv
    for d in chunks/*; do
        for f in \$d/representatives/*; do
            [ -e "\$f" ] || continue
            b=\$(basename \$f)
            cp \$f ${meta.id}/representatives/\$b
            printf '%s\\t%s\\n' "\$b" "\$b" >> ${meta.id}/clusters.tsv
            printf '%s\\trepresentative\\n' "\$b" >> ${meta.id}/genome_status.tsv
        done
    done
    touch versions.yml
    """
}
