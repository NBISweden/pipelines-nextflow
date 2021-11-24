process STATISTICS_CHLOROPLAST {

    publishDir "${params.outdir}/statistics", mode: 'copy', pattern: "chloroplast_significant_matches_statistics.tsv"

    input:
    path statistics_bitfiltered
    path accessions_matchfiltered

    output:
    path "chloroplast_significant_matches_statistics.tsv"

    script:
    """
    LINES=\$(cat $accessions_matchfiltered)
    for line in \$LINES
    do
        grep \$line $statistics_bitfiltered >> chloroplast_significant_matches_statistics.tsv
    done
    """

}