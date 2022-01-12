process FILTER {

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

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
    touch accessions_matchfiltered.tsv  # Create output files
    touch accessions_suspicious.tsv
    echo -e "Accession\\tUnique_matches\\tSpan_fraction\\tScaffold_length\\tClass" >> ${organelle}_draft_statistics_summary.tsv # Create header of summary file
    awk '\$12>$bitscore {print}' $blast_file > statistics_${organelle}.tsv  # Filter matches by bitscore
    awk '{print \$2}' statistics_${organelle}.tsv | sort | uniq > ${organelle}_unique_bitscore.tsv  # Extract all unique scaffold headers
    LINES=\$(cat ${organelle}_unique_bitscore.tsv)  # Initialize loop over unique headers
    for line in \$LINES
    do
        grep \$line statistics_${organelle}.tsv | sort -k 9,9 -n > line_file.tsv    # Grep each row containing the header, and sort for start position
        unique_count=\$(awk '{print \$1}' line_file.tsv | sort | uniq | wc -l)  # Count unique gene matches 
        tot_length=\$(awk 'FNR == 1 {print \$14}' line_file.tsv)    # Extract scaffold length
        awk '{print \$9}' line_file.tsv > positions.tsv # Extract first match position
        awk '{print \$10}' line_file.tsv >> positions.tsv   # Extract last match position
        first_match_pos=\$(sort -n positions.tsv|awk 'FNR == 1 {print \$1}')    # Find first match position
        last_match_pos=\$(sort  -n -r positions.tsv | awk 'FNR == 1 {print \$1}')   # Find last match position
        raw_length_span=\$(((\$last_match_pos-\$first_match_pos)))  # Calculate span of matches 
        length_span=\${raw_length_span#-}   # Remove negative sign if present
        span_fraction=\$(awk "BEGIN {print \$length_span/\$tot_length}")    # Calculate span fraction of total scaffold

        if [ \$unique_count -gt $organelle_gene_matches ] && [ $max_scaffold_length -gt \$tot_length ] && [[ \$span_fraction > $min_span_fraction ]]    # Find organellar scaffolds
        then
            echo \$line >> accessions_matchfiltered.tsv
            echo -e "\$line\\t\$unique_count\\t\$span_fraction\\t\$tot_length\\t${organelle}" >> ${organelle}_draft_statistics_summary.tsv
        elif [ \$unique_count -gt $suspicious_gene_matches ]    # Find suspicious scaffolds
        then
            echo \$line >> accessions_suspicious.tsv
            echo -e "\$line\\t\$unique_count\\t\$span_fraction\\t\$tot_length\\tsuspicious" >> ${organelle}_draft_statistics_summary.tsv
        fi
    done
    """    
}
