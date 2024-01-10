process AGAT_FILTERBYMRNABLASTVALUE {
    tag "${gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::agat=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.2.0--pl5321hdfd78af_0':
        'biocontainers/agat:1.2.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)
    path blast_tbl

    output:
    tuple val(meta), path("*_blast-filtered.gff3"), emit: blast_filtered
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gff.baseName}"
    def VERSION = '1.2.0'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_filter_by_mrnaBlastValue.pl \\
        --gff $gff \\
        --blast $blast_tbl \\
        --outfile ${prefix}_blast-filtered.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}
