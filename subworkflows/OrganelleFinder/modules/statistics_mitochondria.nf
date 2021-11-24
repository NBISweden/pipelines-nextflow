process STATISTICS_MITOCHONDRIA {

    publishDir "${params.outdir}/statistics", mode: 'copy', pattern: "mitochondria_significant_matches_statistics.tsv"

    input:
    path statistics_bitfiltered
    path accessions_matchfiltered

    output:
    path "mitochondria_significant_matches_statistics.tsv"

    script:
    """
    LINES=\$(cat $accessions_matchfiltered)
    for line in \$LINES
    do
        grep \$line $statistics_bitfiltered >> mitochondria_significant_matches_statistics.tsv
    done
    """

}
