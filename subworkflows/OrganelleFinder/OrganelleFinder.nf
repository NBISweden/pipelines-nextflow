#! /usr/bin/env nextflow

nextflow.enable.dsl=2

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
include { MAKEBLASTDB as MAKEBLASTDB_CHLOROPLAST; MAKEBLASTDB as MAKEBLASTDB_MITOCHONDRIA } from './modules/makeblastdb.nf'
include { TBLASTN as TBLASTN_CHLOROPLAST; TBLASTN as TBLASTN_MITOCHONDRIA } from './modules/tblastn.nf'
include { FILTER as FILTER_CHLOROPLAST; FILTER as FILTER_MITOCHONDRIA } from './modules/filter.nf'
include { STATISTICS as STATISTICS_CHLOROPLAST; STATISTICS as STATISTICS_MITOCHONDRIA  } from './modules/statistics.nf'
include { EXTRACT_FINAL  } from './modules/extract_final.nf'
include { EXTRACT_MITOCHONDRIA  } from './modules/extract_mitochondria.nf'
include { PLOTING } from './modules/ploting.nf'
include { MAPPING } from './modules/mapping.nf'
include { DEPTH_FILTER_ANIMAL } from './modules/depth_filter_animal'
include { DEPTH_FILTER_PLANT } from './modules/depth_filter_plant'

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
     params.reads_file                   : ${params.reads_file}
     params.input_type                   : ${params.input_type}

     // Output directory
     params.outdir                       : ${params.outdir}

     // Mitochondria parameters
     params.mit_blast_evalue             : ${params.mit_blast_evalue}
     params.mit_bitscore                 : ${params.mit_bitscore}
     params.mit_significant_gene_matches : ${params.mit_significant_gene_matches}
     params.mit_max_contig_length        : ${params.mit_max_contig_length}
     params.mit_min_span_fraction        : ${params.mit_min_span_fraction}

     // Chloroplast parameters
     params.chl_blast_evalue             : ${params.chl_blast_evalue}
     params.chl_bitscore                 : ${params.chl_bitscore}
     params.chl_significant_gene_matches : ${params.chl_significant_gene_matches}
     params.chl_max_contig_length        : ${params.chl_max_contig_length}
     params.chl_min_span_fraction        : ${params.chl_min_span_fraction}
 """)

workflow {

    main:
        Channel.fromPath(params.genome_assembly, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome_assembly}!\n" }
            .set {genome_assembly}
        Channel.fromPath(params.reference_mitochondria, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find mitochondria file matching ${params.reference_mitochondria}!\n" }
            .set {reference_mitochondria}
        reads_file = Channel.empty()
        if (params.reads_file) {
            Channel.fromPath(params.reads_file, checkIfExists: true)
                .ifEmpty { exit 1, "Cannot find read files matching ${params.reads_file}!\n" }
                .set {reads_file}
        }
        if (params.input_type == 'animal') {
            ANIMAL_ORGANELLE_FINDER(genome_assembly,reference_mitochondria,reads_file)
        } else if (params.input_type == 'plant') {
            Channel.fromPath(params.reference_chloroplast, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find chloroplast file matching ${params.reference_chloroplast}!\n" }
            .set {reference_chloroplast}
            PLANT_ORGANELLE_FINDER(genome_assembly,reference_mitochondria,reference_chloroplast,reads_file)
        }

}


workflow ANIMAL_ORGANELLE_FINDER {

    take:
        genome_assembly
        reference_mitochondria
        reads_file

    main:
        MAKEBLASTDB_MITOCHONDRIA(genome_assembly,genome_assembly.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        MAKEBLASTDB_MITOCHONDRIA.out.mix(genome_assembly.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{blastfiles}
        TBLASTN_MITOCHONDRIA(reference_mitochondria,
            blastfiles, params.mit_blast_evalue)
        FILTER_MITOCHONDRIA(TBLASTN_MITOCHONDRIA.out.output_blast, params.outdir, params.mit_bitscore, params.mit_significant_gene_matches, params.mit_suspicious_gene_matches, params.mit_max_contig_length, params.mit_min_span_fraction, "mitochondria")
        MAPPING(genome_assembly, reads_file)
        DEPTH_FILTER_ANIMAL(MAPPING.out.depth_file, FILTER_MITOCHONDRIA.out.accessions, FILTER_MITOCHONDRIA.out.accessions_suspicious)
        PLOTING(DEPTH_FILTER_ANIMAL.out.collect(), params.outdir)
        STATISTICS_MITOCHONDRIA(FILTER_MITOCHONDRIA.out.statistics, FILTER_MITOCHONDRIA.out.accessions, FILTER_MITOCHONDRIA.out.accessions_suspicious, params.outdir, "mitochondria")
        EXTRACT_FINAL(genome_assembly,
            FILTER_MITOCHONDRIA.out.accessions, params.outdir, "mitochondria")
    emit:
        blast_result = EXTRACT_FINAL.out[0]
}

workflow PLANT_ORGANELLE_FINDER {

    take:
        genome_assembly
        reference_mitochondria
        reference_chloroplast
        reads_file

    main:
        MAKEBLASTDB_MITOCHONDRIA(genome_assembly,genome_assembly.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        MAKEBLASTDB_MITOCHONDRIA.out.mix(genome_assembly.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{mitochondria_blastfiles}
        TBLASTN_MITOCHONDRIA(reference_mitochondria,
            mitochondria_blastfiles, params.mit_blast_evalue)
        FILTER_MITOCHONDRIA(TBLASTN_MITOCHONDRIA.out.output_blast, params.outdir, params.mit_bitscore, params.mit_significant_gene_matches, params.mit_suspicious_gene_matches, params.mit_max_contig_length, params.mit_min_span_fraction, "mitochondria")
        STATISTICS_MITOCHONDRIA(FILTER_MITOCHONDRIA.out.statistics, FILTER_MITOCHONDRIA.out.accessions, FILTER_MITOCHONDRIA.out.accessions_suspicious, params.outdir, "mitochondria")
        EXTRACT_MITOCHONDRIA(genome_assembly,
            FILTER_MITOCHONDRIA.out.accessions, params.outdir)
        MAKEBLASTDB_CHLOROPLAST(EXTRACT_MITOCHONDRIA.out.no_mitochondria,EXTRACT_MITOCHONDRIA.out.no_mitochondria.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        MAKEBLASTDB_CHLOROPLAST.out.mix(EXTRACT_MITOCHONDRIA.out.no_mitochondria.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{chloroplast_blastfiles}
        TBLASTN_CHLOROPLAST(reference_chloroplast,chloroplast_blastfiles, params.chl_blast_evalue)
        FILTER_CHLOROPLAST(TBLASTN_CHLOROPLAST.out.output_blast, params.outdir, params.chl_bitscore, params.chl_significant_gene_matches, params.chl_suspicious_gene_matches, params.chl_max_contig_length, params.chl_min_span_fraction, "chloroplast")
        MAPPING(genome_assembly, reads_file)
        DEPTH_FILTER_PLANT(MAPPING.out.depth_file, FILTER_MITOCHONDRIA.out.accessions, FILTER_CHLOROPLAST.out.accessions, FILTER_MITOCHONDRIA.out.accessions_suspicious, FILTER_CHLOROPLAST.out.accessions_suspicious )
        PLOTING(DEPTH_FILTER_PLANT.out.collect(), params.outdir)
        STATISTICS_CHLOROPLAST(FILTER_CHLOROPLAST.out.statistics, FILTER_CHLOROPLAST.out.accessions, FILTER_CHLOROPLAST.out.accessions_suspicious, params.outdir, "chloroplast")
        EXTRACT_FINAL(EXTRACT_MITOCHONDRIA.out.no_mitochondria, FILTER_CHLOROPLAST.out.accessions, params.outdir, "chloroplast")
    emit:
        blast_result = EXTRACT_FINAL.out[0]
}

workflow.onComplete {
    log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}
