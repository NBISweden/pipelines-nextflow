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
    echo -e "Accession\\tUnique_matches\\tCoding_fraction\\tSpan_fraction" >> chloroplast_statistics_summary.tsv
    awk '\$12>${params.chl_bitscore} {print}' $blast_file > statistics_chloroplast.tsv
    awk '{print \$2}' statistics_chloroplast.tsv | sort | uniq > unique_bitscore.tsv
    LINES=\$(cat unique_bitscore.tsv)
    for line in \$LINES
    do
        grep \$line statistics_chloroplast.tsv | sort -k 9,9 -n > line_file.tsv
        unique_count=\$(awk '{print \$1}' line_file.tsv | sort | uniq | wc -l)
        coding_length=\$(awk '{print \$4*3}' line_file.tsv | awk '{s+=\$1}END{print s}')
        tot_length=\$(awk 'FNR == 1 {print \$14}' line_file.tsv)
        coding_fraction=\$(awk "BEGIN {print \$coding_length/\$tot_length}")
        first_match_pos=\$(awk 'FNR == 1 {print \$9}' line_file.tsv)
        last_match_pos=\$(sort -k 9,9 -n -r line_file.tsv | awk 'FNR == 1 {print \$10}')
        length_span=\$(((\$last_match_pos-\$first_match_pos)))
        span_fraction=\$(awk "BEGIN {print \$length_span/\$tot_length}")
        echo -e "\$line\\t\$unique_count\\t\$coding_fraction\\t\$span_fraction" >> chloroplast_statistics_summary.tsv
        if [ \$unique_count -gt ${params.mit_significant_gene_matches} ]
        then
            echo \$line >> accessions_matchfiltered.tsv
        fi
    done
    """    
}



//awk '{print \$2}' statistics_bitfiltered.tsv | sort | uniq -c | sort -r | awk '\$1>${params.significant_gene_matches} {print}' | awk '{print \$2}' > accessions_matchfiltered.tsv


//| awk '{print \$1}' | sort | uniq | wc -l
