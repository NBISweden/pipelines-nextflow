process CUSTOM_RANKMODELS {
    tag "$prefix"
    label 'process_single'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0':
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(augustus_logs)

    output:
    path "*sweep_summary.tsv", emit: summary
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: meta.id
    """
    ( 
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
            "locus_distance" \\
            "model_selection_value" \\
            "exon_sensitivity" \\
            "exon_specificity" \\
            "nucleotide_sensitivity" \\
            "nucleotide_specificity" \\
            "gene_sensitivity" \\
            "gene_specificity" \\
            "genes"
        for LOG in $augustus_logs; do
            printf "%d\\t%f\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%d\\n" \\
                \$( grep -Eo "LD[0-9]+" <<< \$LOG | cut -c 3- ) \\
                \$( grep -Eo "AED[0-9.]+" <<< \$LOG  | cut -c 4- ) \\
                \$( grep "^exon level" \$LOG | grep -Eo "[0-9.]+" | sed -n "4p" ) \\
                \$( grep "^exon level" \$LOG | grep -Eo "[0-9.]+" | sed -n "5p" ) \\
                \$( grep "^nucleotide level" \$LOG | grep -Eo "[0-9.]+" | sed -n "1p" ) \\
                \$( grep "^nucleotide level" \$LOG | grep -Eo "[0-9.]+" | sed -n "2p" ) \\
                \$( grep "^gene level" \$LOG | grep -Eo "[0-9.]+" | sed -n "6p" ) \\
                \$( grep "^gene level" \$LOG | grep -Eo "[0-9.]+" | sed -n "7p" ) \\
                \$( grep "Training gene count:" \$LOG | grep -Eo "[0-9.]+" )
        done
    ) > ${prefix}_sweep_summary.tsv

    # genes=\$( grep -c "LOCUS" codingGeneFeatures.filter.longest_cds.complete.good_distance_blast-filtered.gbk.train )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$( grep -V |& sed '2!d;s/.*v//;s/ .*//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_sweep_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$( grep -V |& sed '2!d;s/.*v//;s/ .*//' )
    END_VERSIONS
    """
}
