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

params.trimming_tool = 'trimmomatic'

params.trimmomatic_adapter_path = '/path/to/trimmomatic/adapters.fasta'
params.trimmomatic_clip_options = 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

params.hisat2_options = ''

params.stringtie_options = ''

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
     trimming_tool              : ${params.trimming_tool}

 Trimmomatic parameters
     trimmomatic_adapter_path   : ${params.trimmomatic_adapter_path}
     trimmomatic_clip_options   : ${params.trimmomatic_clip_options}

 Hisat2 parameters
     hisat2_options             : ${params.hisat2_options}

 StringTie parameters
     stringtie_options          : ${params.stringtie_options}

 """

if( params.trimming_tool == 'trimmomatic' && !file(params.trimmomatic_adapter_path).exists() ){
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
        trimmomatic(reads)
        hisat2_index(genome)
        hisat2(trimmomatic.out[0].mix(trimmomatic.out[2]),
                hisat2_index.out.collect())
        stringtie(hisat2.out[0])
        multiqc(fastqc.out.collect(),
                trimmomatic.out[3].collect(),
                hisat2.out[2].collect(),
                stringtie.out[1].collect())

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

process trimmomatic {

    tag "Adapter-trimming: $sample_id"
    publishDir "${params.outdir}/Trimmomatic", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*_paired_*.fastq.gz') optional true
    tuple val(sample_id), path('*_unpaired_*.fastq.gz') optional true
    tuple val(sample_id), path('*_trimmed.fastq.gz') optional true
    path 'trimmomatic.log'

    script:
    if (params.single_end) {
    """
    trimmomatic SE -threads ${task.cpus} $reads \\
        ${sample_id}_trimmed.fastq.gz \\
        ILLUMINACLIP:${params.trimmomatic_adapter_path}:2:30:10 \\
        ${params.trimmomatic_clip_options} 2> trimmomatic.log
    """
    } else {
    """
    trimmomatic PE -threads ${task.cpus} $reads \\
        ${sample_id}_paired_R1.fastq.gz ${sample_id}_unpaired_R1.fastq.gz \\
        ${sample_id}_paired_R2.fastq.gz ${sample_id}_unpaired_R2.fastq.gz \\
        ILLUMINACLIP:${params.trimmomatic_adapter_path}:2:30:10 \\
        ${params.trimmomatic_clip_options} 2> trimmomatic.log
    """
    }

}

process hisat2_index {

    tag "Indexing $genome_fasta"
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

    tag "Aligning reads (${sample_id}) to genome"
    publishDir "${params.outdir}/Hisat2_alignments", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path hisat2_index_files

    output:
    path "${sample_id}_sorted_alignment.bam"
    path 'splicesite.txt'
    path "*hisat2_summary.txt"

    script:
    hisat2_basename = hisat2_index_files[0].toString() - ~/.\d.ht2l?/
    if (params.single_end){
    """
    hisat2 ${params.hisat2_options} --novel-splicesite-outfile splicesite.txt \\
        --new-summary --summary-file ${sample_id}.hisat2_summary.txt \\
        -p ${task.cpus} -x $hisat2_basename -U $reads | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}_sorted_alignment.bam -
    """
    } else {
    """
    hisat2 ${params.hisat2_options} --novel-splicesite-outfile splicesite.txt \\
        --new-summary --summary-file ${sample_id}.hisat2_summary.txt \\
        -p ${task.cpus} -x $hisat2_basename -1 ${reads[0]} -2 ${reads[1]} | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}_sorted_alignment.bam -
    """
    }

}

process stringtie {

    tag "${sorted_bam_file.name}"
    publishDir "${params.outdir}/Stringtie_transcripts", mode: 'copy'

    input:
    path sorted_bam_file

    output:
    path "${sorted_bam_file.name}_transcripts.gtf"
    path ".command.log"

    script:
    """
    stringtie ${sorted_bam_file} -l ${sorted_bam_file.name} -o ${sorted_bam_file.name}_transcripts.gtf -p ${task.cpus} ${params.stringtie_options}
    """

}

process multiqc {

    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path('fastqc/*')
    path('trimmomatic/trimmomatic_log*')
    path('hisat2/*')
    path('stringtie/stringtie_log*')

    output:
    path "*multiqc_report.html"
    path "*_data"

    script:
    """
    multiqc . -m fastqc -m trimmomatic -m hisat2
    """
}


workflow.onComplete {
    log.info ( workflow.success ? "\nTranscript assembly complete!\n" : "Oops .. something went wrong\n" )
}
