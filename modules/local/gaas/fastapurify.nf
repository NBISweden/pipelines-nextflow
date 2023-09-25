process GAAS_FASTAPURIFY {
    tag "${fasta.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::gaas=1.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gaas:1.2.0--pl526r35_0':
        'biocontainers/gaas:1.2.0--pl526r35_0' }"

    input:
    path fasta

    output:
    path "*_purified/*_purified.fa", emit: fasta
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def VERSION = '1.2.0'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gaas_fasta_purify.pl \\
        $args \\
        --infile $fasta \\
        --output ${prefix}_purified

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gaas: $VERSION
    END_VERSIONS
    """
}
