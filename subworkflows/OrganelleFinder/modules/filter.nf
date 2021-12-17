process FILTER {

    input:
    path blast_file
    val bitscore
    val organelle_gene_matches
    val suspicious_gene_matches
    val max_scaffold_length
    val min_span_fraction
    val organelle

    output:
    path "statistics_${organelle}.tsv", emit: statistics
    path "accessions_matchfiltered.tsv", emit: accessions
    path  "accessions_suspicious.tsv", emit: accessions_suspicious
    path "${organelle}_statistics_summary.tsv"

    script:
    """
    touch accessions_matchfiltered.tsv
    touch accessions_suspicious.tsv
    echo -e "Accession\\tUnique_matches\\tSpan_fraction\\tScaffold_length\\tClass" >> ${organelle}_draft_statistics_summary.tsv
    awk '\$12>$bitscore {print}' $blast_file > statistics_${organelle}.tsv
    awk '{print \$2}' statistics_${organelle}.tsv | sort | uniq > ${organelle}_unique_bitscore.tsv
    LINES=\$(cat ${organelle}_unique_bitscore.tsv)
    for line in \$LINES
    do
        grep \$line statistics_${organelle}.tsv | sort -k 9,9 -n > line_file.tsv
        unique_count=\$(awk '{print \$1}' line_file.tsv | sort | uniq | wc -l)
        coding_length=\$(awk '{print \$4*3}' line_file.tsv | awk '{s+=\$1}END{print s}')
        tot_length=\$(awk 'FNR == 1 {print \$14}' line_file.tsv)
        awk '{print \$9}' line_file.tsv > positions.tsv
        awk '{print \$10}' line_file.tsv >> positions.tsv
        first_match_pos=\$(sort -n positions.tsv|awk 'FNR == 1 {print \$1}')
        last_match_pos=\$(sort  -n -r positions.tsv | awk 'FNR == 1 {print \$1}')
        raw_length_span=\$(((\$last_match_pos-\$first_match_pos)))
        length_span=\${raw_length_span#-}
        span_fraction=\$(awk "BEGIN {print \$length_span/\$tot_length}")

        if [ \$unique_count -gt $organelle_gene_matches ] && [ $max_scaffold_length -gt \$tot_length ] && [[ \$span_fraction > $min_span_fraction ]]
        then
            echo \$line >> accessions_matchfiltered.tsv
            echo -e "\$line\\t\$unique_count\\t\$span_fraction\\t\$tot_length\\t${organelle}" >> ${organelle}_draft_statistics_summary.tsv
        elif [ \$unique_count -gt $suspicious_gene_matches ]
        then
            echo \$line >> accessions_suspicious.tsv
            echo -e "\$line\\t\$unique_count\\t\$span_fraction\\t\$tot_length\\tsuspicious" >> ${organelle}_draft_statistics_summary.tsv
        fi
    done
    column -t ${organelle}_draft_statistics_summary.tsv > "${organelle}_statistics_summary.tsv"
    """    
}
