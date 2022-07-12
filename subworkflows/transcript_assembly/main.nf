// include { FASTQC              } from "$projectDir/modules/nf-core/modules/fastqc/main"
// include { HISAT2_ALIGN        } from "$projectDir/modules/nf-core/modules/hisat2/align/main"
// include { HISAT2_BUILD        } from "$projectDir/modules/nf-core/modules/hisat2/build/main"
// include { FASTP               } from "$projectDir/modules/nf-core/modules/fastp/main"
// include { STRINGTIE_STRINGTIE } from "$projectDir/modules/nf-core/modules/stringtie/stringtie/main"
// include { MULTIQC             } from "$projectDir/modules/nf-core/modules/multiqc/main"

workflow TRANSCRIPT_ASSEMBLY {

    main:
    log.info """
        Transcript assembly using Hisat2/Stringtie workflow
        ===================================================
    """
    Channel.fromFilePairs( params.reads, size: params.single_end ? 1 : 2, checkIfExists: true )
        // .ifEmpty { error "Cannot find reads matching ${params.reads}!\n" }
        .set { reads }
    Channel.fromPath( params.genome, checkIfExists: true )
        // .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" }
        .set { genome }

    FASTQC ( reads )
    HISAT2_INDEX ( genome )
    if ( !params.skip_trimming ){
            FASTP( reads )
            FASTP.out.trimmed_reads.set { trimmed_reads }
            FASTP.out.json.set { trimming_logs }
    } else {
        reads.set { trimmed_reads }
    }
    HISAT2 ( trimmed_reads, HISAT2_INDEX.out.indices.collect() )
    STRINGTIE ( HISAT2.out.trimmed_reads )
    MULTIQC( 
        FASTQC.out.logs.mix( trimming_logs, HISAT2.out.json ).collect(),
        params.multiqc_config
    )

}

process FASTQC {

    tag "$sample_id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"

    input:
    tuple val(sample_id), path(reads)

    output:
    path ("fastqc_${sample_id}_logs"), emit: logs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_id}"
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -t ${task.cpus} -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process FASTP {

    tag "$sample"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::fastp=0.23.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0' :
        'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0' }"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('*fastp-trimmed*.fastq.gz'), emit: trimmed_reads
    path "${sample}_fastp.json"                        , emit: json

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample}"
    if (params.single_end) {
        """
        fastp ${params.fastp_options} -w ${task.cpus} -i ${reads} \\
            -o ${sample}_fastp-trimmed.fastq.gz \\
            --json ${sample}_fastp.json
        """
    } else {
        """
        fastp ${params.fastp_options} -w ${task.cpus} -i ${reads[0]} -I ${reads[1]} \\
            -o ${sample}_fastp-trimmed_R1.fastq.gz \\
            -O ${sample}_fastp-trimmed_R2.fastq.gz \\
            --json ${sample}_fastp.json
        """
    }
}

process HISAT2_INDEX {

    tag "$genome_fasta"
    label 'process_high'
    label 'process_high_memory'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::hisat2=2.2.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:0e773bb207600fcb4d38202226eb20a33c7909b6-0' :
        'quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:0e773bb207600fcb4d38202226eb20a33c7909b6-0' }"

    input:
    path(genome_fasta)

    output:
    path('*.ht2'), emit: indices

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${genome_fasta.baseName}"
    """
    hisat2-build -p ${task.cpus} $genome_fasta ${genome_fasta.baseName}.hisat2_index
    """
}

process HISAT2 {

    tag "$sample"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::hisat2=2.2.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:0e773bb207600fcb4d38202226eb20a33c7909b6-0' :
        'quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:0e773bb207600fcb4d38202226eb20a33c7909b6-0' }"

    input:
    tuple val(sample), path(reads)
    path hisat2_index_files

    output:
    path "${sample}_sorted.bam"    , emit: bam
    path "${sample}_splicesite.txt", emit: splicesites
    path "*hisat2-summary.txt"     , emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample}"
    index_basename = hisat2_index_files[0].toString() - ~/.\d.ht2l?/
    if (params.single_end){
        """
        hisat2 ${params.hisat2_options} --novel-splicesite-outfile ${sample}_splicesite.txt \\
            --new-summary --summary-file ${sample}.hisat2-summary.txt \\
            -p ${task.cpus} -x $index_basename -U $reads | \\
            samtools sort -@ ${task.cpus} -o ${sample}_sorted.bam -
        """
    } else {
        """
        hisat2 ${params.hisat2_options} --novel-splicesite-outfile ${sample}_splicesite.txt \\
            --new-summary --summary-file ${sample}.hisat2-summary.txt \\
            -p ${task.cpus} -x $index_basename -1 ${reads[0]} -2 ${reads[1]} | \\
            samtools sort -@ ${task.cpus} -o ${sample}_sorted.bam -
        """
    }
}

process STRINGTIE {

    tag "$bam"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2' :
        'quay.io/biocontainers/stringtie:2.2.1--hecb563c_2' }"

    input:
    path bam

    output:
    path "${bam.baseName}-transcripts.gtf", emit gtf_transcripts

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bam.baseName}"
    """
    stringtie \\
        $bam \\
        -l ${bam.baseName} \\
        -o ${bam.baseName}-transcripts.gtf \\
        -p ${task.cpus} \\
        ${params.stringtie_options}
    """
}

process MULTIQC {

    label 'process_single'

    conda (params.enable_conda ? 'bioconda::multiqc=1.13a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13a--pyhdfd78af_1' :
        'quay.io/biocontainers/multiqc:1.13a--pyhdfd78af_1' }"

    input:
    path log_files
    path multiqc_config

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "multiqc"
    """
    multiqc . -c $multiqc_config
    """
}
