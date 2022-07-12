process SNAP_TRAINING {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::snap=2013_11_29" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snap:2013_11_29--hec16e2b_4':
        'quay.io/biocontainers/snap:2013_11_29--hec16e2b_4' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path "*.hmm", emit: training_model
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    ann_file = training_files.find { it =~ /.ann$/ }
    dna_file = training_files.find { it =~ /.dna$/ }
    """
    fathom -categorize ${params.flank_region_size} ${ann_file} ${dna_file}
    fathom -export ${params.flank_region_size} -plus uni.ann uni.dna
    forge export.ann export.dna
    hmm-assembler.pl "$species_label" . > "${species_label}.hmm"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snap: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
