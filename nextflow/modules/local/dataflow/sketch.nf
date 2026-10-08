// Sketch every genome of a run with sourmash (data-channel form).
//
// One sourmash call per genome writes sketches/<record name>.sig.zip, a zip
// collection of three DNA signatures (k=21, k=31, k=51; scaled=1000; no
// abundance), the same layout as the sketches/ directory that the repgenr
// genome-writing stages fill. The record name is the file name without its
// FASTA suffix (gzip suffixes are stripped first), as in
// repgenr.core.contracts.record_name. Each signature is written to a hidden
// partial file and renamed when complete, so a stopped task leaves no partial
// sketch under its final name. The calls run in parallel, up to task.cpus.
//
// Unlike the other dataflow modules this process calls sourmash directly, not
// through repgenr, so it declares the image and conda package itself. They
// match repgenr.core.sourmash (SOURMASH_IMAGE, SOURMASH_CONDA); a unit test
// keeps the two in step.

process SKETCH {
    label 'process_medium'
    tag "${meta.id}"

    conda 'bioconda::sourmash'
    container 'quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0'

    input:
    tuple val(meta), path(genomes, stageAs: 'genomes/*')

    output:
    tuple val(meta), path('sketches'), emit: sketches
    path 'versions.yml'              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p sketches

    # One line per genome: the FASTA path and its record name. Dotfiles
    # (for example macOS ._ companions) are not genomes.
    : > genomes.list
    for f in genomes/*; do
        b=\$(basename "\$f")
        case "\$b" in
            .*) continue ;;
        esac
        name=""
        for s in .fasta.gz .fna.gz .fa.gz .fasta .fa .fna .fas; do
            case "\$b" in
                *"\$s") name="\${b%"\$s"}"; break ;;
            esac
        done
        [ -n "\$name" ] || continue
        printf '%s\\t%s\\n' "\$f" "\$name" >> genomes.list
    done

    # One sourmash process per genome, up to task.cpus at a time.
    tr '\\t' '\\n' < genomes.list | xargs -n 2 -P ${task.cpus} sh -c '
        sourmash sketch dna -p k=21,k=31,k=51,scaled=1000 ${args} \\
            --name "\$2" -o "sketches/.\$2.partial.sig.zip" "\$1" \\
        && mv "sketches/.\$2.partial.sig.zip" "sketches/\$2.sig.zip"
    ' sh

    rm -f genomes.list

    printf '"%s":\n    sourmash: "%s"\n' "${task.process}" \
        "\$(sourmash --version 2>&1 | sed 's/^sourmash //')" > versions.yml
    """

    stub:
    """
    mkdir -p sketches
    for f in genomes/*; do
        b=\$(basename "\$f")
        case "\$b" in
            .*) continue ;;
        esac
        for s in .fasta.gz .fna.gz .fa.gz .fasta .fa .fna .fas; do
            case "\$b" in
                *"\$s") touch "sketches/\${b%"\$s"}.sig.zip"; break ;;
            esac
        done
    done
    touch versions.yml
    """
}
