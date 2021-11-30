process STATISTICS {

    publishDir "${outdir}/statistics", mode: 'copy', pattern: "${organelle}_significant_matches_statistics.tsv"

    input:
    path statistics_bitfiltered
    path accessions_matchfiltered
    path outdir
    val organelle

    output:
    path "${organelle}_significant_matches_statistics.tsv"

    script:
    """
    grep -f $accessions_matchfiltered $statistics_bitfiltered > ${organelle}_significant_matches_statistics.tsv
    
    """

}


//LINES=\$(cat $accessions_matchfiltered)
//    for line in \$LINES
//    do
//        grep \$line $statistics_bitfiltered >> ${organelle}_significant_matches_statistics.tsv
//    done