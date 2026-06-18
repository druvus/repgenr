// Turn the rooted tree into FlexTaxD relations (data-channel form).
//
// Bridges like PHYLO: stages the tree, the merged representatives (for
// --include-dereplicated), the outgroup and selection into a task-local workdir,
// runs the existing `repgenr tree2tax`, and emits tree2tax.tsv + genomes_map.tsv.

process TREE2TAX {
    label 'process_low'
    tag "tree2tax"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tree
    path reps_dir
    path outgroup
    path selection

    output:
    path "tree2tax.tsv"   , emit: tree2tax
    path "genomes_map.tsv", emit: genomes_map
    path "versions.yml"   , emit: versions

    script:
    """
    mkdir -p wd/tree wd/derep wd/outgroup
    cp ${tree} wd/tree/tree.nwk
    [ -f ${reps_dir}/clusters.tsv ] && cp ${reps_dir}/clusters.tsv wd/derep/clusters.tsv

    for f in ${outgroup}; do
        [ -e "\$f" ] && cp "\$f" wd/outgroup/
    done
    awk -F'\\t' 'NR>1 && \$5==1 {print \$1}' ${selection} | head -1 > wd/outgroup_accession.txt

    repgenr ${params.repgenr_opts} tree2tax -wd wd ${params.tree2tax_args}
    cp wd/tree2tax.tsv .
    cp wd/genomes_map.tsv .

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
