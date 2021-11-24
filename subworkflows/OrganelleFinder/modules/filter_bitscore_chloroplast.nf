process FILTER_BITSCORE_CHLOROPLAST {

    publishDir "${params.outdir}/statistics", mode: 'copy', pattern: "chloroplast_statistics_summary.tsv"

    input:
    path blast_file

    output:
    path "statistics_chloroplast.tsv", emit: statistics
    path "accessions_matchfiltered.tsv", emit: accessions
    path "chloroplast_statistics_summary.tsv"

    script:
    """
    echo -e "Accession\\tUnique_matches" >> chloroplast_statistics_summary.tsv
    awk '\$12>${params.chl_bitscore} {print}' $blast_file > statistics_chloroplast.tsv
    awk '{print \$2}' statistics_chloroplast.tsv | sort | uniq > unique_bitscore.tsv
    LINES=\$(cat unique_bitscore.tsv)
    for line in \$LINES
    do
        unique_count=\$(grep \$line statistics_chloroplast.tsv | awk '{print \$1}' | sort | uniq | wc -l)
        echo -e "\$line\\t\$unique_count" >> chloroplast_statistics_summary.tsv
        if [ \$unique_count -gt ${params.chl_significant_gene_matches} ]
        then
            echo \$line >> accessions_matchfiltered.tsv
        fi
    done
    """    
}



//awk '{print \$2}' statistics_bitfiltered.tsv | sort | uniq -c | sort -r | awk '\$1>${params.significant_gene_matches} {print}' | awk '{print \$2}' > accessions_matchfiltered.tsv


//| awk '{print \$1}' | sort | uniq | wc -l
