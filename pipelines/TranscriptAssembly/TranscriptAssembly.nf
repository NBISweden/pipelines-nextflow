#! /usr/bin/env nextflow

nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.reads = "/path/to/reads_{1,2}.fastq.gz"
params.genome = "/path/to/genome.fa"
params.single_end = false
params.outdir = "results"

trimming_tools = [ 'fastp', 'trimmomatic' ]
params.skip_trimming = false
params.trimmer = 'fastp'

params.fastp_options = ' -Q -L'

params.trimmomatic_adapter_path = '/path/to/trimmomatic/adapters.fasta'
params.trimmomatic_clip_options = 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

params.hisat2_options = ''

params.stringtie_options = ''

params.multiqc_config = "$baseDir/config/multiqc_conf.yml"

log.info """
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Transcript assembly using Hisat2/Stringtie workflow
 ===================================================

 General Parameters
     genome                     : ${params.genome}
     reads                      : ${params.reads}
     single_end                 : ${params.single_end}
     outdir                     : ${params.outdir}
     skip_trimming              : ${params.skip_trimming}
     trimmer                    : ${params.trimmer}

 Fastp parameters
     fastp_options              : ${params.fastp_options}

 Trimmomatic parameters
     trimmomatic_adapter_path   : ${params.trimmomatic_adapter_path}
     trimmomatic_clip_options   : ${params.trimmomatic_clip_options}

 Hisat2 parameters
     hisat2_options             : ${params.hisat2_options}

 StringTie parameters
     stringtie_options          : ${params.stringtie_options}

 MultiQC parameters
     multiqc_config             : ${params.multiqc_config}

 """

if( !params.skip_trimming && !(params.trimmer in trimming_tools) ){
    exit 1, "Error: ${params.trimmer} is not a valid trimming tool option.\n Please provide a valid option from ${trimming_tools}.\n"
}
if( !params.skip_trimming && params.trimmer == 'trimmomatic' && !file(params.trimmomatic_adapter_path).exists() ){
    exit 1, "The adapter file '${params.trimmomatic_adapter_path}' does not exist!\n"
}

workflow {

    main:
        reads = Channel.fromFilePairs(params.reads, size: params.single_end ? 1 : 2, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find reads matching ${params.reads}!\n" }
        genome = Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
        transcript_assembly(reads,genome)
}

workflow transcript_assembly {

    get:
        reads
        genome

    main:
        fastqc(reads)
        hisat2_index(genome)
        if(!params.skip_trimming){
            if (params.trimmer == 'fastp') {
                fastp(reads)
                fastp.out[0].set {trimmed_reads}
                fastp.out[1].set {trimming_logs}
            } else if (params.trimmer == 'trimmomatic') {
                trimmomatic(reads)
                trimmomatic.out[0].set {trimmed_reads}
                trimmomatic.out[2].set {trimming_logs}
            }
        } else {
            reads.set {trimmed_reads}
        }
        hisat2(trimmed_reads,hisat2_index.out.collect())
        stringtie(hisat2.out[0])
        fastqc.out.mix(trimming_logs).mix(hisat2.out[2]).set {logs}
        multiqc(logs.collect(),params.multiqc_config)

}

process fastqc {

    tag "$sample_id"
    publishDir "${params.outdir}/FastQC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path ("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -t ${task.cpus} -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """

}

process fastp {

    tag "$sample"
    publishDir "${params.outdir}/Fastp", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('*fastp-trimmed*.fastq.gz')
    path "${sample}_fastp.json"

    script:
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


process trimmomatic {

    tag "$sample"
    publishDir "${params.outdir}/Trimmomatic", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('*trimmomatic-trimmed*.fastq.gz')
    tuple val(sample), path('*trimmomatic-unpaired*.fastq.gz') optional true
    path 'trimmomatic.log'

    script:
    if (params.single_end) {
    """
    trimmomatic SE -threads ${task.cpus} $reads \\
        ${sample}_trimmomatic-trimmed.fastq.gz \\
        ILLUMINACLIP:${params.trimmomatic_adapter_path}:2:30:10 \\
        ${params.trimmomatic_clip_options} 2> ${sample}_trimmomatic.log
    """
    } else {
    """
    trimmomatic PE -threads ${task.cpus} $reads \\
        ${sample}_trimmomatic-trimmed_R1.fastq.gz ${sample}_trimmomatic-unpaired_R1.fastq.gz \\
        ${sample}_trimmomatic-trimmed_R2.fastq.gz ${sample}_trimmomatic-unpaired_R2.fastq.gz \\
        ILLUMINACLIP:${params.trimmomatic_adapter_path}:2:30:10 \\
        ${params.trimmomatic_clip_options} 2> ${sample}_trimmomatic.log
    """
    }

}

process hisat2_index {

    tag "$genome_fasta"
    publishDir "${params.outdir}/Hisat2_indicies", mode: 'copy'

    input:
    path(genome_fasta)

    output:
    path('*.ht2')

    script:
    """
    hisat2-build -p ${task.cpus} $genome_fasta ${genome_fasta.baseName}.hisat2_index
    """
}

process hisat2 {

    tag "$sample"
    publishDir "${params.outdir}/Hisat2_alignments", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path hisat2_index_files

    output:
    path "${sample}_sorted.bam"
    path "${sample}_splicesite.txt"
    path "*hisat2-summary.txt"

    script:
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

process stringtie {

    tag "$bam"
    publishDir "${params.outdir}/Stringtie_transcripts", mode: 'copy'

    input:
    path bam

    output:
    path "${bam.baseName}-transcripts.gtf"

    script:
    """
    stringtie $bam -l ${bam.baseName} -o ${bam.baseName}-transcripts.gtf \\
        -p ${task.cpus} ${params.stringtie_options}
    """

}

process multiqc {

    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path log_files
    path multiqc_config
    // path('fastqc/*')
    // path('trimmomatic/trimmomatic_log*')
    // path('hisat2/*')
    // path('stringtie/stringtie_log*')

    output:
    path "*multiqc_report.html"
    path "*_data"

    script:
    """
    multiqc . -c $multiqc_config
    """
}


workflow.onComplete {
    log.info ( workflow.success ? "\nTranscript assembly complete!\n" : "Oops .. something went wrong\n" )
}
