process EXTRACT_FINAL {

    conda "${task.ext.enable_conda ? 'bioconda::seqtk:1.3' : '' }"
    container "${workflow.containerEngine == 'singularity' &&
                  !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
              'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    path assembly
    path matched_accessions
    val organelle

    output:
    path "${assembly.baseName}_${organelle}.fna"
    path "${assembly.baseName}_nuclear.fna"

    script: 
    """
    grep '>' $assembly | cut -c2- | grep -v -f $matched_accessions > nuclear_accessions.lst
    seqtk subseq $assembly $matched_accessions > ${assembly.baseName}_${organelle}.fna
    seqtk subseq $assembly nuclear_accessions.lst > ${assembly.baseName}_nuclear.fna
    """
}
