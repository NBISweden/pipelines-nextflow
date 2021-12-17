process STATISTICS {

    input:
    path statistics_bitfiltered
    path accessions_matchfiltered
    path accessions_suspicious
    val organelle

    output:
    path "${organelle}_organelle_matches_statistics.tsv"
    path "${organelle}_suspicious_matches_statistics.tsv"

    script:
    """
    grep -f $accessions_matchfiltered $statistics_bitfiltered | column -t > ${organelle}_organelle_matches_statistics.tsv
    grep -f $accessions_suspicious $statistics_bitfiltered | column -t > ${organelle}_suspicious_matches_statistics.tsv
    
    """

}
