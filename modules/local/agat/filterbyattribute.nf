process AGAT_FILTERBYATTRIBUTE {
    tag "${mrna_gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::agat=0.9.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(mrna_gff)

    output:
    tuple val(meta), path("*.filter.gff"), emit: selected_models
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${mrna_gff.baseName}"
    def VERSION = '0.9.2'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_filter_feature_by_attribute_value.pl \\
        --gff ${mrna_gff} \\
        $args \\
        -o ${prefix}.filter.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}
