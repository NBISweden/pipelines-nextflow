process HISAT2_BUILD {

    tag "$fasta"
    label 'process_high'
    label 'process_high_memory'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::hisat2=2.2.0 bioconda::samtools=1.15.1" : null)
    container "nbisweden/hisat2:2.1.0"

    input:
    path( fasta )

    output:
    path "hisat2"       , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def VERSION = '2.2.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir hisat2
    hisat2-build \\
        -p $task.cpus \\
        $args \\
        $fasta \\
        hisat2/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: $VERSION
    END_VERSIONS
    """
}
