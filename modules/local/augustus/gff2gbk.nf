process AUGUSTUS_GFF2GBK {
    tag "${gff.baseName}"
    label 'process_single'

    conda "bioconda::augustus=3.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augustus:3.4.0--pl5321h5f9f3d9_6':
        'biocontainers/augustus:3.4.0--pl5321h5f9f3d9_6' }"

    input:
    tuple val(meta), path(gff)
    path genome

    output:
    tuple val(meta), path ("*.gbk"), emit: gbk
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gff.baseName}"
    """
    gff2gbSmallDNA.pl \\
        $gff \\
        $genome \\
        $args \\
        ${prefix}.gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        augustus: \$( augustus | sed '1!d; s/.*(//; s/).*//' )
    END_VERSIONS
    """
}
