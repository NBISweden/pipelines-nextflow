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
include { MAKEBLASTDB_CHLOROPLAST   } from './modules/makeblastdb_chloroplast.nf'
include { MAKEBLASTDB_MITOCHONDRIA   } from './modules/makeblastdb_mitochondria.nf'
include { TBLASTN_CHLOROPLAST   } from './modules/tblastn_chloroplast.nf'
include { TBLASTN_MITOCHONDRIA   } from './modules/tblastn_mitochondria.nf'
include { FILTER_BITSCORE_CHLOROPLAST  } from './modules/filter_bitscore_chloroplast.nf'
include { FILTER_BITSCORE_MITOCHONDRIA  } from './modules/filter_bitscore_mitochondria.nf'
include { STATISTICS_CHLOROPLAST  } from './modules/statistics_chloroplast.nf'
include { STATISTICS_MITOCHONDRIA  } from './modules/statistics_mitochondria.nf'
include { EXTRACT  } from './modules/extract.nf'
include { EXTRACT_CHLOROPLAST  } from './modules/extract_chloroplast.nf'
include { EXTRACT_MITOCHONDRIA  } from './modules/extract_mitochondria.nf'

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
     params.genome_assembly              : ${params.genome_assembly}
     params.reference_mitochondria       : ${params.reference_mitochondria}
     params.reference_chloroplast        : ${params.reference_chloroplast}
     params.input_type                   : ${params.input_type}

     // Output directory
     params.outdir                       : ${params.outdir}
    
     // Mitochondria parameters
     params.mit_blast_evalue             : ${params.mit_blast_evalue}
     params.mit_bitscore                 : ${params.mit_bitscore}
     params.mit_significant_gene_matches : ${params.mit_significant_gene_matches}

     // Chloroplast parameters
     params.chl_blast_evalue             : ${params.chl_blast_evalue}
     params.chl_bitscore                 : ${params.chl_bitscore}
     params.chl_significant_gene_matches : ${params.chl_significant_gene_matches}
 """)

workflow {

    main:
        Channel.fromPath(params.genome_assembly, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome_assembly}!\n" }
            .set {genome_assembly}
        Channel.fromPath(params.reference_mitochondria, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find blast database files matching ${params.reference_mitochondria}!\n" }
            .set {reference_mitochondria}
        if (params.input_type == 'animal') {
            ANIMAL_ORGANELLE_FINDER(genome_assembly,reference_mitochondria)
        } else if (params.input_type == 'plant') {
            Channel.fromPath(params.reference_chloroplast, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find blast database files matching ${params.reference_chloroplast}!\n" }
            .set {reference_chloroplast}
            PLANT_ORGANELLE_FINDER(genome_assembly,reference_mitochondria,reference_chloroplast)
        }

}

workflow ANIMAL_ORGANELLE_FINDER {

    take:
        genome_assembly
        reference_mitochondria

    main:
        MAKEBLASTDB_MITOCHONDRIA(genome_assembly,genome_assembly.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        MAKEBLASTDB_MITOCHONDRIA.out.mix(genome_assembly.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{blastfiles}
        TBLASTN_MITOCHONDRIA(reference_mitochondria,
            blastfiles)
        FILTER_BITSCORE_MITOCHONDRIA(TBLASTN_MITOCHONDRIA.out.output_blast)
        STATISTICS_MITOCHONDRIA(FILTER_BITSCORE_MITOCHONDRIA.out.statistics, FILTER_BITSCORE_MITOCHONDRIA.out.accessions)
        EXTRACT(genome_assembly,
            FILTER_BITSCORE_MITOCHONDRIA.out.accessions)
    emit:
        blast_result = EXTRACT.out[0]
}

workflow PLANT_ORGANELLE_FINDER {

    take:
        genome_assembly
        reference_mitochondria
        reference_chloroplast

    main:
        MAKEBLASTDB_MITOCHONDRIA(genome_assembly,genome_assembly.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        MAKEBLASTDB_MITOCHONDRIA.out.mix(genome_assembly.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{mitochondria_blastfiles}
        TBLASTN_MITOCHONDRIA(reference_mitochondria,
            mitochondria_blastfiles)
        FILTER_BITSCORE_MITOCHONDRIA(TBLASTN_MITOCHONDRIA.out.output_blast)
        STATISTICS_MITOCHONDRIA(FILTER_BITSCORE_MITOCHONDRIA.out.statistics, FILTER_BITSCORE_MITOCHONDRIA.out.accessions)
        EXTRACT_MITOCHONDRIA(genome_assembly,
            FILTER_BITSCORE_MITOCHONDRIA.out.accessions)
        MAKEBLASTDB_CHLOROPLAST(EXTRACT_MITOCHONDRIA.out.no_mitochondria,EXTRACT_MITOCHONDRIA.out.no_mitochondria.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        MAKEBLASTDB_CHLOROPLAST.out.mix(EXTRACT_MITOCHONDRIA.out.no_mitochondria.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{chloroplast_blastfiles}
        TBLASTN_CHLOROPLAST(reference_chloroplast,chloroplast_blastfiles)
        FILTER_BITSCORE_CHLOROPLAST(TBLASTN_CHLOROPLAST.out.output_blast)
        STATISTICS_CHLOROPLAST(FILTER_BITSCORE_CHLOROPLAST.out.statistics, FILTER_BITSCORE_CHLOROPLAST.out.accessions)
        EXTRACT_CHLOROPLAST(EXTRACT_MITOCHONDRIA.out.no_mitochondria, FILTER_BITSCORE_CHLOROPLAST.out.accessions)
    emit:
        blast_result = EXTRACT_CHLOROPLAST.out[0]
}

workflow.onComplete {
    log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}