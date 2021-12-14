process MAKEBLASTDB {

    conda "${task.ext.enable_conda ? 'bioconda::blast:2.12.0' : '' }"
    container "${workflow.containerEngine == 'singularity' &&
                  !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
              'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    label 'blast'

    input:
    path genome
    val state

    output:
    path "*.fna*"

    when:
    state == 'DBFILES_ABSENT'

    script:
    """
    makeblastdb -in $genome -dbtype nucl
    """

}
