#! /usr/bin/env nextflow

nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
Vendela, William and Viktor
 */
// Test5 by Vendela
// Test of conflicts
// Test by William

params.genome = "../../../GCF_003254395.2_Amel_HAv3.1_genomic.fna"
params.outdir = "results"

// blastx parameters
params.blast_db_fasta = '../../../c_picta_ref_mito.fna'
params.blast_evalue = '1e-21'

log.info("""
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Functional annotation workflow
 ===================================

 General parameters
     genome                         : ${params.genome}
     outdir                         : ${params.outdir}

 Blast parameters
     blast_db_fasta                 : ${params.blast_db_fasta}
     params.blast_evalue            : ${params.blast_evalue}
 """)

workflow {

    main:
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
            .set {genome}
        Channel.fromPath(params.blast_db_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find blast database files matching ${params.blast_db_fasta}!\n" }
            .set {blastdb}
        organelle_finder(genome,blastdb)

}

workflow organelle_finder {

    take:
        genome
        blastdb

    main:
        makeblastdb(blastdb,blastdb.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        makeblastdb.out.mix(blastdb.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{blastfiles}
        blastx(genome,
            blastfiles)
    emit:
        blast_result = blastx.out.collect()
}



process makeblastdb {

    label 'blast'

    input:
    path reference
    val state

    output:
    path "*.fna*"

    when:
    state == 'DBFILES_ABSENT'

    script:
    """
    makeblastdb -in $reference -dbtype prot
    """

}

process blastx {

    publishDir "${params.outdir}", mode: 'copy', pattern: "${fasta_file.baseName}_blast.tsv"

    label 'blast'

    input:
    path fasta_file
    path blastdb

    output:
    path "${fasta_file.baseName}_blast.tsv"

    script:
    database = blastdb.find { it =~ /\.f(ast|n)?a$/ }
    """
    blastx -query $fasta_file -db ${database} -evalue ${params.blast_evalue} -outfmt 6 -out ${fasta_file.baseName}_blast.tsv
    """

}


workflow.onComplete {
    log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}
