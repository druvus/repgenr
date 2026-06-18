// Download selected genomes (data-channel form).
//
// Consumes the selection.tsv channel, downloads the genomes with
// `repgenr genome-fetch` (no manifest, no shared workdir), and emits the genome
// FASTAs (and the outgroup, if any) as channel outputs the dereplication
// scatter consumes.

process GENOME {
    label 'process_medium'
    tag "genome"
    publishDir "${params.outdir}/genomes", mode: 'copy'

    input:
    path selection

    output:
    path "out/genomes/*" , emit: genomes
    path "out/outgroup/*", emit: outgroup, optional: true
    path "versions.yml"  , emit: versions

    script:
    """
    repgenr ${params.repgenr_opts} genome-fetch --selection ${selection} --out out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repgenr: \$(repgenr --version | sed 's/repgenr //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out/genomes out/outgroup
    tail -n +2 ${selection} | while IFS=\$'\\t' read -r acc fam gen sp og fname; do
        [ -z "\$fname" ] && continue
        if [ "\$og" = "1" ]; then
            printf '>x\\nACGT\\n' > out/outgroup/\$fname
        else
            printf '>x\\nACGT\\n' > out/genomes/\$fname
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repgenr: stub
    END_VERSIONS
    """
}
