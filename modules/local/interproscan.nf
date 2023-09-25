process INTERPROSCAN {
    tag "${protein_fasta.baseName}"
    label 'process_single'

    conda "bioconda::interproscan=5.55_88.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/interproscan:5.55_88.0--hec16e2b_1':
        'biocontainers/interproscan:5.55_88.0--hec16e2b_1' }"

    input:
    tuple val(meta), path(protein_fasta)

    output:
    tuple val(meta), path('*.tsv'), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${protein_fasta.baseName}"
    """
    interproscan.sh \\
        $args \\
        -i $protein_fasta \\
        -o ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$( interproscan.sh --version | sed '1!d; s/.*version //' )
    END_VERSIONS
    """
}
