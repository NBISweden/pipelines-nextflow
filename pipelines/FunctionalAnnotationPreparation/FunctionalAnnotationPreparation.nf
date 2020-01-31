#! /usr/bin/env nextflow

nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.gff_annotation = "/path/to/annotation.gff"
params.genome = "/path/to/genome.fasta"
params.outdir = "results"

params.codon_table = 1

params.records_per_file = 1000

params.blast_db_fasta = '/path/to/protein/database.fasta'

params.interproscan_db = 'all'

params.merge_annotation_identifier = 'ID'

log.info("""
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Functional annotation input preparation workflow
 ===================================

 General parameters
     gff_annotation                 : ${params.gff_annotation}
     genome                         : ${params.genome}
     outdir                         : ${params.outdir}

 Parallelisation parameters
     records_per_file               : ${params.records_per_file}

 Gff2Protein parameters
     codon_table                    : ${params.codon_table}

 Blast parameters
     blast_db_fasta                 : ${params.blast_db_fasta}

 Interproscan parameters
     interproscan_db                : ${params.interproscan_db}

 Merge functional annotation parameters
     merge_annotation_identifier    : ${params.merge_annotation_identifier}

 """)

workflow {

    main:
        annotation = Channel.fromPath(params.gff_annotation, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find gff file matching ${params.gff_annotation}!\n" }
        genome = Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
        blastdb = Channel.fromPath("${params.blast_db_fasta}{,.p*}", checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find blast database files matching ${params.blast_db_fasta}{,.p*}" }
        functional_annotation_input_preparation(annotation,genome,blastdb)

}

workflow functional_annotation_input_preparation {

    get:
        gff_file
        genome
        blastdb

    main:
        gff2protein(gff_file,genome.collect())
        blastp(gff2protein.out.splitFasta(by: params.records_per_file, file: true),
            blastdb.collect())
        interproscan(gff2protein.out.splitFasta(by: params.records_per_file, file: true))
        merge_functional_annotation(gff_file,
            blastp.out.collectFile(name:'blast_merged.tsv').collect(),
            interproscan.out.collectFile(name:'interproscan_merged.tsv').collect(),
            blastdb.collect())
}

process gff2protein {

    label 'AGAT'

    input:
    path gff_file
    path genome_fasta

    output:
    path "${gff_file.baseName}_proteins.fasta"

    script:
    """
    agat_sp_extract_sequences.pl -o ${gff_file.baseName}_proteins.fasta -f $genome_fasta \\
        -p -cfs -cis -ct ${params.codon_table} --gff $gff_file
    """
    // agat_sp_extract_sequences.pl is a script from AGAT

}

process blastp {

    input:
    path fasta_file
    path blastdb

    output:
    path "${fasta_file.baseName}_blast.tsv"

    script:
    database = blastdb[0].toString() - ~/.p\w\w$/
    """
    blastp -query $fasta_file -db ${database} -num_threads ${task.cpus} \\
        -outfmt 6 -out ${fasta_file.baseName}_blast.tsv
    """

}

process interproscan {

    input:
    path protein_fasta

    output:
    path '*.tsv'

    script:
    applications = { params.interproscan_db ? "-appl ${params.interproscan_db}" : '' }
    """
    interproscan.sh ${applications} -i $protein_fasta -o ${protein_fasta.baseName}.tsv \\
        -f TSV --iprlookup --goterms -pa -dp -t p
    """

}

process merge_functional_annotation {

    publishDir "${params.outdir}/blast_tsv", mode: 'copy', pattern: 'blast_merged.tsv'
    publishDir "${params.outdir}/interproscan_tsv", mode: 'copy', pattern: 'interproscan_merged.tsv'
    publishDir "${params.outdir}/final_annotation", mode: 'copy', pattern: "${gff_annotation.baseName}_plus-functional-annotation.gff"
    label 'AGAT'

    input:
    path gff_annotation
    path merged_blast_results
    path merged_interproscan_results
    path blast_files

    output:
    path "${gff_annotation.baseName}_plus-functional-annotation.gff"

    script:
    """
    agat_sp_manage_functional_annotation.pl -f ${gff_annotation} \\
        -b ${merged_blast_results} -i ${merged_interproscan_results} \\
        -db ${params.blast_db_fasta} -id ${params.merge_annotation_identifier} \\
        -o ${gff_annotation.baseName}_plus-functional-annotation.gff
    """
    // agat_sp_manage_functional_annotation.pl is a script from AGAT

}

workflow.onComplete {
    log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}
