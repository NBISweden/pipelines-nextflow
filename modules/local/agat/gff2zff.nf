process AGAT_GFF2ZFF {
    tag "${annotation}"
    label 'process_single'

    conda "bioconda::agat=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.2.0--pl5321hdfd78af_0':
        'biocontainers/agat:1.2.0--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path (annotation)
    path genome

    output:
    tuple val(meta), path ("*.{ann,dna}"), emit: zff
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${genome.baseName}"
    def VERSION = '1.2.0'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_convert_sp_gff2zff.pl \\
        --gff $annotation \\
        --fasta $genome \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}
