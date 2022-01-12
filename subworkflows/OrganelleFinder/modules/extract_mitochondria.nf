process EXTRACT_MITOCHONDRIA {

    conda "${task.ext.enable_conda ? 'bioconda::seqtk:1.3' : '' }"
    container "${workflow.containerEngine == 'singularity' &&
                  !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
              'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"
    input:
    path assembly
    path matched_accessions

    output:
    path "${assembly.baseName}_mitochondria.fna"
    path "${assembly.baseName}.fna", emit: no_mitochondria

    script: 
    """
    grep '>' $assembly | cut -c2- | grep -v -f $matched_accessions > nuclear_accessions.lst
    seqtk subseq $assembly $matched_accessions > ${assembly.baseName}_mitochondria.fna
    seqtk subseq $assembly nuclear_accessions.lst > ${assembly.baseName}_no_mit.fna
    """
}
