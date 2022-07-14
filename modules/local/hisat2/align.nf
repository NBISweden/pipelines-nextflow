process HISAT2_ALIGN {

    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::hisat2=2.2.0 bioconda::samtools=1.15.1" : null)
    container "nbisweden/hisat2:2.1.0"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.bam")          , emit: bam
    tuple val(meta), path("*.log")          , emit: summary
    tuple val(meta), path("*splicesite.txt"), emit: splicesites, optional: true
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.2.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    if (meta.single_end) {
        """
        INDEX=\$( find -L ./ -name "*.1.ht2" | sed 's/.1.ht2//' )
        hisat2 \\
            -x \$INDEX \\
            -U $reads \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $args \\
            | samtools sort --threads task.cpus -o ${prefix}.bam -

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: $VERSION
            samtools: \$( samtools --version | sed '1!d; s/samtools //' )
        END_VERSIONS
        """
    } else {
        """
        INDEX=\$( find -L ./ -name "*.1.ht2" | sed 's/.1.ht2//' )
        hisat2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $args \\
            | samtools sort --threads $task.cpus -o ${prefix}.bam -

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: $VERSION
            samtools: \$( samtools --version | sed '1!d; s/samtools //' )
        END_VERSIONS
        """
    }
}
