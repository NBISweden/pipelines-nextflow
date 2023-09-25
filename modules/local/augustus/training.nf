process AUGUSTUS_TRAINING {
    tag "$species_label"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::augustus=3.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augustus:3.4.0--pl5321h5f9f3d9_6':
        'biocontainers/augustus:3.4.0--pl5321h5f9f3d9_6' }"

    input:
    tuple val(meta), path (training_file)
    path test_file
    val species_label

    output:
    tuple val(meta), path ("${species_label}"), emit: training_model
    tuple val(meta), path ("*_run.log")       , emit: log
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${species_label}"
    """
    : \${AUGUSTUS_CONFIG_PATH:=/usr/local/config}
    cp -rv \${AUGUSTUS_CONFIG_PATH} .
    export AUGUSTUS_CONFIG_PATH="\$PWD/config"
    new_species.pl --species=$species_label
    etraining --species=$species_label $training_file
    augustus --species=$species_label $test_file | tee ${prefix}_run.log
    mv config/species/${species_label} .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        augustus: \$( augustus | sed '1!d; s/.*(//; s/).*//' )
    END_VERSIONS
    """
}
