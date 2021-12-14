process DEPTH_FILTER_ANIMAL {

    input:
    path depth_file
    path accessions
    path accessions_suspicious

    output:
    path "*depth.tsv*"

    script:
    """
    cat $accessions_suspicious>>$accessions
    LINES=\$(cat $accessions)
    for line in \$LINES
    do
        grep \$line $depth_file | sort -k 2,2 -n > \${line}_depth.tsv
    done
    """    
}
