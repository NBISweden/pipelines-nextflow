process EXTRACT_MITOCHONDRIA {

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
