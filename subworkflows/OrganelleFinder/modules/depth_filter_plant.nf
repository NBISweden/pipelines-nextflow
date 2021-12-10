process DEPTH_FILTER_PLANT {

    input:
    path depth_file
    path accessions_mit 
    path accessions_chl

    output:
    path "*depth.tsv*"

    script:
    """
    LINES=\$(cat $accessions_mit)
    for line in \$LINES
    do
        echo \$line >> accessions.tsv
    done
    LINES=\$(cat $accessions_chl)
    for line in \$LINES
    do
        echo \$line >> accessions.tsv
    done
    awk '{print \$1}' accessions.tsv | sort | uniq > unique_accessions.tsv
    LINES=\$(cat unique_accessions.tsv)
    for line in \$LINES
    do
        grep \$line $depth_file | sort -k 2,2 -n > \${line}_depth.tsv
    done
    """    
}


