// Turn the rooted tree into FlexTaxD relations (data-channel form).
//
// Calls the stateless `repgenr tree2tax-relations` step directly on the staged
// channel files: the tree, the merged representatives' clusters.tsv (for
// --include-dereplicated) and the outgroup directory. Emits tree2tax.tsv and
// genomes_map.tsv as channel outputs; there is no shared working directory.
// Tool flags arrive as task.ext.args from conf/modules.config; publishing is
// configured there too.

process TREE2TAX {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(tree), path(reps_dir), path(outgroup, stageAs: 'outgroup/*'), path(outgroup_accession)

    output:
    tuple val(meta), path("tree2tax.tsv")   , emit: tree2tax
    tuple val(meta), path("genomes_map.tsv"), emit: genomes_map
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} tree2tax-relations \\
        --tree ${tree} \\
        --clusters ${reps_dir}/clusters.tsv \\
        --outgroup-dir outgroup \\
        --outgroup-accession ${outgroup_accession} \\
        -o . ${args} \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    printf 'child\\tparent\\n' > tree2tax.tsv
    for f in ${reps_dir}/representatives/*; do
        leaf=\$(basename \$f | sed 's/\\.[^.]*\$//')
        printf '%s\\troot\\n' "\$leaf" >> tree2tax.tsv
    done
    : > genomes_map.tsv
    for f in ${reps_dir}/representatives/*; do
        leaf=\$(basename \$f | sed 's/\\.[^.]*\$//')
        acc=\$(echo \$leaf | awk -F_ '{print \$(NF-1)"_"\$NF}')
        printf '%s\\t%s\\n' "\$acc" "\$leaf" >> genomes_map.tsv
    done
    touch versions.yml
    """
}
