// Score and classify a batch of assemblies (data-channel form).
//
// Runs `repgenr genome-qc` once over every assemble-run directory of a run:
// CheckM2 when params.checkm2_db is set, the sourmash classifier when
// params.gtdb_sketch is set (both reach the task as task.ext.args from
// conf/modules.config). Emits the qc directory (quality.tsv,
// classification.tsv) the gather step reads. The subworkflow skips this
// process when neither database is configured.

process GENOME_QC {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(assemblies, stageAs: 'assemblies/*')

    output:
    tuple val(meta), path("qc"), emit: qc
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} genome-qc \\
        --assemblies assemblies \\
        --out qc \\
        ${args} \\
        --threads ${task.cpus} \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    mkdir -p qc
    printf 'run_accession\\tcompleteness\\tcontamination\\n' > qc/quality.tsv
    printf 'run_accession\\ttaxonomy\\trank\\tscore\\tdb_version\\n' > qc/classification.tsv
    for d in assemblies/*/; do
        run=\$(basename \$d)
        [ -e "\$d/assembly.ok" ] || continue
        printf '%s\\t99.00\\t0.50\\n' "\$run" >> qc/quality.tsv
        printf '%s\\td__Bacteria;f__Francisellaceae;g__Francisella;s__Francisella tularensis\\tspecies\\t0.9\\tstub\\n' "\$run" >> qc/classification.tsv
    done
    touch versions.yml
    """
}
