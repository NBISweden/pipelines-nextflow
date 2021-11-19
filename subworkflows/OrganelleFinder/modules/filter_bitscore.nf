process FILTER_BITSCORE {

    input:
    path blast_file

    output:
    path "statistics_bitfiltered.tsv", emit: statistics
    path "accessions_matchfiltered.tsv", emit: accessions

    script:
    """
    awk '\$12>${params.bitscore} {print}' $blast_file > statistics_bitfiltered.tsv
    awk '{print \$2}' statistics_bitfiltered.tsv | sort | uniq -c | sort -r | awk '\$1>${params.significant_gene_matches} {print}' | awk '{print \$2}' > accessions_matchfiltered.tsv
    """    
}