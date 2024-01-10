process AGAT_MANAGEFUNCTIONALANNOTATION {
    tag "${gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::agat=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.2.0--pl5321hdfd78af_0':
        'biocontainers/agat:1.2.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)
    path merged_blast_results
    path merged_interproscan_results
    path blast_db

    output:
    tuple val(meta), path("*_plus-functional-annotation.gff"), emit: gff
    tuple val(meta), path("*.tsv", includeInputs: true)      , emit: tsv
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gff.baseName}"
    def VERSION = '1.2.0'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_manage_functional_annotation.pl \\
        $args \\
        -f ${gff} \\
        -b ${merged_blast_results} \\
        -i ${merged_interproscan_results} \\
        -db ${blast_db} \\
        -o ${prefix}_plus-functional-annotation.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}
