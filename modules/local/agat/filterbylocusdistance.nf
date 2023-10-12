process AGAT_FILTERBYLOCUSDISTANCE {
    tag "${coding_gene_features_gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::agat=0.9.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path (coding_gene_features_gff)

    output:
    tuple val(meta), path ("*.good_distance.gff"), emit: distanced_models
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${coding_gene_features_gff.baseName}"
    def VERSION = '0.9.2'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_filter_by_locus_distance.pl \\
        $args \\
        --gff ${coding_gene_features_gff} \\
        -o ${prefix}.good_distance.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}
