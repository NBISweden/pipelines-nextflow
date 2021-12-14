process DEPTH_FILTER_PLANT {

    input:
    path depth_file
    path accessions_mit 
    path accessions_chl

    output:
    path "*depth.tsv*"

    script:
    """
    cat $accessions_mit >> accessions.tsv
    cat $accessions_chl >> accessions.tsv    
    sort accessions.tsv | uniq > unique_accessions.tsv
    LINES=\$(cat unique_accessions.tsv)
    for line in \$LINES
    do
        grep \$line $depth_file | sort -k 2,2 -n > \${line}_depth.tsv
    done
    """    
}


