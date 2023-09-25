process AGAT_SEPARATEBYRECORD {
    tag "${gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path (gff)

    output:
    tuple val(meta), path("maker_results_noAbinitio_clean/mrna.gff"), emit: transcripts
    tuple val(meta), path("maker_results_noAbinitio_clean/*")       , emit: all
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.9.2'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_separate_by_record_type.pl \\
        -g ${gff} \\
        -o maker_results_noAbinitio_clean
    if test -f maker_results_noAbinitio_clean/mrna.gff && test -f maker_results_noAbinitio_clean/transcript.gff; then
        agat_sp_merge_annotations.pl \\
            --gff maker_results_noAbinitio_clean/mrna.gff \\
            --gff maker_results_noAbinitio_clean/transcript.gff \\
            --out merged_transcripts.gff
        mv merged_transcripts.gff maker_results_noAbinitio_clean/mrna.gff
    elif test -f maker_results_noAbinitio_clean/transcript.gff; then
        cp maker_results_noAbinitio_clean/transcript.gff maker_results_noAbinitio_clean/mrna.gff
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}
