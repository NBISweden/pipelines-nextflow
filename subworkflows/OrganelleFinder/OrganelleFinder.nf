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


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { MAKEBLASTDB   } from './modules/makeblastdb.nf'
include { TBLASTN   } from './modules/tblastn.nf'
include { FILTER_BITSCORE  } from './modules/filter_bitscore.nf'
include { STATISTICS  } from './modules/statistics.nf'
include { EXTRACT  } from './modules/extract.nf'

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
     reference_fasta                 : ${params.reference_fasta}
     params.blast_evalue            : ${params.blast_evalue}
 Filter parameters
     params.bitscore                : ${params.bitscore}
     params.significant_gene_matches: ${params.significant_gene_matches}
 """)

workflow {

    main:
        Channel.fromPath(params.genome_assembly, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome_assembly}!\n" }
            .set {genome_assembly}
        Channel.fromPath(params.reference_organelle, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find blast database files matching ${params.reference_organelle}!\n" }
            .set {reference_organelle}
        ORGANELLE_FINDER(genome_assembly,reference_organelle)

}

workflow ORGANELLE_FINDER {

    take:
        genome_assembly
        reference_organelle

    main:
        MAKEBLASTDB(genome_assembly,genome_assembly.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        MAKEBLASTDB.out.mix(genome_assembly.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{blastfiles}
        TBLASTN(reference_organelle,
            blastfiles)
        FILTER_BITSCORE(TBLASTN.out.output_blast)
        STATISTICS(FILTER_BITSCORE.out.statistics, FILTER_BITSCORE.out.accessions)
        EXTRACT(genome_assembly,
            FILTER_BITSCORE.out.accessions)
    emit:
        blast_result = EXTRACT.out[0]
}

workflow.onComplete {
    log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}