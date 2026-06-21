// Turn the rooted tree into FlexTaxD relations (data-channel form).
//
// Calls the stateless `repgenr tree2tax-relations` step directly on the staged
// channel files: the tree, the merged representatives' clusters.tsv (for
// --include-dereplicated) and the outgroup directory. Emits tree2tax.tsv and
// genomes_map.tsv as channel outputs; there is no shared working directory.

process TREE2TAX {
    label 'process_low'
    tag "tree2tax"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tree
    path reps_dir
    path outgroup, stageAs: 'outgroup/*'
    path outgroup_accession

    output:
    path "tree2tax.tsv"   , emit: tree2tax
    path "genomes_map.tsv", emit: genomes_map
    path "versions.yml"   , emit: versions

    script:
    """
    repgenr ${params.repgenr_opts} tree2tax-relations \\
        --tree ${tree} \\
        --clusters ${reps_dir}/clusters.tsv \\
        --outgroup-dir outgroup \\
        --outgroup-accession ${outgroup_accession} \\
        -o . ${params.tree2tax_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repgenr: \$(repgenr --version | sed 's/repgenr //')
    END_VERSIONS
    """

    stub:
    """
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
