process SNAP_TRAINING {
    tag "$species_label"
    label 'process_single'

    conda "bioconda::snap=2013_11_29"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snap:2013_11_29--hec16e2b_4':
        'biocontainers/snap:2013_11_29--hec16e2b_4' }"

    input:
    tuple val(meta), path (training_files)
    val species_label

    output:
    tuple val(meta), path ("*.hmm"), emit: training_model
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${species_label}"
    ann_file = training_files.find { it =~ /.ann$/ }
    dna_file = training_files.find { it =~ /.dna$/ }
    """
    fathom \\
        $args \\
        ${ann_file} \\
        ${dna_file}
    fathom \\
        $args2 \\
        uni.ann \\
        uni.dna
    forge export.ann export.dna
    hmm-assembler.pl "$species_label" . > "${species_label}.hmm"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fathom: \$( fathom |& sed '2!d ; s/.*version //; s/)//' )
        forge: \$( forge |& sed '2!d; s/.*version //; s/)//' )
        snap: \$( snap |& sed '2!d; s/.*version //; s/)//' )
    END_VERSIONS
    """
}
