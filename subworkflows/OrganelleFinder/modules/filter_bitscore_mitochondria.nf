process FILTER_BITSCORE_MITOCHONDRIA {
    
    publishDir "${params.outdir}/statistics", mode: 'copy', pattern: "mitochondria_statistics_summary.tsv"

    input:
    path blast_file

    output:
    path "statistics_mitochondria.tsv", emit: statistics
    path "accessions_matchfiltered.tsv", emit: accessions
    path "mitochondria_statistics_summary.tsv"

    script:
    """
    echo -e "Accession\\tUnique_matches" >> mitochondria_statistics_summary.tsv
    awk '\$12>${params.mit_bitscore} {print}' $blast_file > statistics_mitochondria.tsv
    awk '{print \$2}' statistics_mitochondria.tsv | sort | uniq > unique_bitscore.tsv
    LINES=\$(cat unique_bitscore.tsv)
    for line in \$LINES
    do
        unique_count=\$(grep \$line statistics_mitochondria.tsv | awk '{print \$1}' | sort | uniq | wc -l)
        echo -e "\$line\\t\$unique_count" >> mitochondria_statistics_summary.tsv
        if [ \$unique_count -gt ${params.mit_significant_gene_matches} ]
        then
            echo \$line >> accessions_matchfiltered.tsv
        fi
    done
    """    
}



//awk '{print \$2}' statistics_bitfiltered.tsv | sort | uniq -c | sort -r | awk '\$1>${params.significant_gene_matches} {print}' | awk '{print \$2}' > accessions_matchfiltered.tsv


//| awk '{print \$1}' | sort | uniq | wc -l
