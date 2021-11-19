process STATISTICS {

    publishDir "${params.outdir}/statistics", mode: 'copy', pattern: "statistics_significant_matches.tsv"

    input:
    path statistics_bitfiltered
    path accessions_matchfiltered

    output:
    path "statistics_significant_matches.tsv"

    script:
    """
    LINES=\$(cat $accessions_matchfiltered)
    for line in \$LINES
    do
        grep \$line $statistics_bitfiltered >> statistics_significant_matches.tsv
    done
    """

}