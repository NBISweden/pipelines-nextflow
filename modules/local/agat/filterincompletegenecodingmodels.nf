process AGAT_FILTERINCOMPLETEGENECODINGMODELS {
    tag "${coding_gene_features_gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path coding_gene_features_gff
    path genome

    output:
    path "*.complete.gff", emit: complete_gene_models
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${coding_gene_features_gff.baseName}"
    def VERSION = '0.9.2'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_filter_incomplete_gene_coding_models.pl \\
        --gff ${coding_gene_features_gff} \\
        -f ${genome} \\
        -o ${prefix}.complete.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}
