process AGAT_EXTRACTSEQUENCES {
    tag "${gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path (gff)
    path genome

    output:
    tuple val(meta), path ("${gff.baseName}_proteins.fasta"), emit: proteins
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${gff.baseName}"
    def VERSION = '0.9.2'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_extract_sequences.pl \\
        $args \\
        --g $gff \\
        -f $genome \\
        -o ${prefix}_proteins.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}
