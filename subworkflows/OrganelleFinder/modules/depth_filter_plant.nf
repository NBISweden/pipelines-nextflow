process DEPTH_FILTER_PLANT {

    input:
    path depth_file
    path accessions_mit 
    path accessions_chl
    path accessions_suspicious_mit
    path accessions_suspicious_chl

    output:
    path "*depth.tsv*"

    script:
    """
    cat $accessions_mit \\
        $accessions_chl \\
        $accessions_suspicious_mit \\
        $accessions_suspicious_chl | \\
        sort -u | \\
        while read -r line; do
            grep \$line $depth_file | sort -k 2,2 -n > \${line}_depth.tsv
        done
    """    
}


