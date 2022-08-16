process AUGUSTUS_GFF2GBK {
    tag "${gff.baseName}"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::augustus=3.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augustus:3.4.0--pl5321h5f9f3d9_6':
        'quay.io/biocontainers/augustus:3.4.0--pl5321h5f9f3d9_6' }"

    input:
    path gff
    path genome

    output:
    path "*.gbk"       , emit: gbk
    path "versions.yml", emit: versions

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
