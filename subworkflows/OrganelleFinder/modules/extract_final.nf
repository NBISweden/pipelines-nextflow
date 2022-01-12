process EXTRACT_FINAL {

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
