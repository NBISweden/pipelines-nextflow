process EXTRACT_MITOCHONDRIA {

    //publishDir "${outdir}", mode: 'copy', pattern: "${assembly.baseName}_mitochondria.fna"

    input:
    path assembly
    path matched_accessions
    //path outdir

    output:
    path "${assembly.baseName}_mitochondria.fna"
    path "${assembly.baseName}.fna", emit: no_mitochondria

    script: 
    """
    cp $assembly ${assembly.baseName}_no_mit.fna
    LINES=\$(cat $matched_accessions)
    for line in \$LINES
    do
        grep -n '>' ${assembly.baseName}_no_mit.fna > row_number
        grep -A 1 \$line row_number > header_rows
        grep -oP  '.*?(?=:>)' header_rows > numbers_file
        start_index=\$(head -n 1 numbers_file)
        next_index=\$(tail -n 1 numbers_file)
        if [ \$(wc -l numbers_file | awk '{print \$1}') -eq 1 ]
        then
            end_of_file=\$(wc -l ${assembly.baseName}_no_mit.fna | awk '{print \$1}')
            next_index=\$(((\$end_of_file+1)))
        fi
        row_count=\$(((\$next_index-\$start_index)))
        end_index=\$(((\$next_index-1)))
        head -n \$end_index ${assembly.baseName}_no_mit.fna | tail -n \$row_count >> ${assembly.baseName}_mitochondria.fna
        sed -e \$start_index,\$end_index\\d ${assembly.baseName}_no_mit.fna > intermediate_file.fna
        mv intermediate_file.fna ${assembly.baseName}_no_mit.fna
    done
    mv ${assembly.baseName}_no_mit.fna ${assembly.baseName}.fna
    """
}