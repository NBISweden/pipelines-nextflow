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

// tblastn parameters
params.reference_fasta = '../../../c_picta_ref_mito.fna'
params.blast_evalue = '1e-21'

// filter parameters
params.bitscore = 150
params.significant_gene_matches = 3

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
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
            .set {genome_assembly}
        Channel.fromPath(params.reference_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find blast database files matching ${params.reference_fasta}!\n" }
            .set {reference_organelle}
        organelle_finder(genome_assembly,reference_organelle)

}

workflow organelle_finder {

    take:
        genome_assembly
        reference_organelle

    main:
        makeblastdb(genome_assembly,genome_assembly.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        makeblastdb.out.mix(genome_assembly.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{blastfiles}
        tblastn(reference_organelle,
            blastfiles)
        filter_bitscore(tblastn.out.collect())
        statistics(filter_bitscore.out.statistics, filter_bitscore.out.accessions)
        extract(genome_assembly,
            filter_bitscore.out.accessions)
    emit:
        blast_result = extract.out[0]
}

process makeblastdb {

    label 'blast'

    input:
    path genome
    val state

    output:
    path "*.fna*"

    when:
    state == 'DBFILES_ABSENT'

    script:
    """
    makeblastdb -in $genome -dbtype nucl
    """

}

process tblastn {

    label 'blast'

    input:
    path reference_organelle
    path blastdb

    output:
    path "output_blast.tsv"

    script:
    database = blastdb.find { it =~ /\.f(ast|n)?a$/ }
    """
    tblastn -query $reference_organelle -db ${database} -evalue ${params.blast_evalue} -outfmt 6 -out output_blast.tsv
    """

}

process filter_bitscore {

    input:
    path blast_file

    output:
    path "statistics_bitfiltered.tsv", emit: statistics
    path "accessions_matchfiltered.tsv", emit: accessions

    script:
    """
    awk '\$12>${params.bitscore} {print}' $blast_file > statistics_bitfiltered.tsv
    awk '{print \$2}' statistics_bitfiltered.tsv | sort | uniq -c | sort -r | awk '\$1>${params.significant_gene_matches} {print}' | awk '{print \$2}' > accessions_matchfiltered.tsv
    """    
}

process statistics {

    publishDir "${params.outdir}/statistics", mode: 'copy', pattern: "statistics_significant_matches.tsv"

    input:
    path statistics_bitfiltered
    path accessions_matchfiltered

    output:
    path "statistics_significant_matches.tsv"

    script:
    """
    LINES=\$(cat $accessions_matchfiltered)
    for line in \$LINES
    do
        grep \$line $statistics_bitfiltered >> statistics_significant_matches.tsv
    done
    """

}

process extract {

    publishDir "${params.outdir}", mode: 'copy', pattern: "${assembly.baseName}_mitochondria.fna"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${assembly.baseName}_nuclear.fna"

    input:
    path assembly
    path matched_accessions

    output:
    path "${assembly.baseName}_mitochondria.fna"
    path "${assembly.baseName}_nuclear.fna"

    script: 
    """
    cp $assembly ${assembly.baseName}_nuclear.fna
    LINES=\$(cat $matched_accessions)
    for line in \$LINES
    do
        grep -n '>' ${assembly.baseName}_nuclear.fna > row_number
        grep -A 1 \$line row_number > header_rows
        grep -oP  '.*?(?=:)' header_rows > numbers_file
        start_index=\$(head -n 1 numbers_file)
        next_index=\$(tail -n 1 numbers_file)
        if [ \$(wc -l numbers_file | awk '{print \$1}') -eq 1 ]
        then
            end_of_file=\$(wc -l ${assembly.baseName}_nuclear.fna | awk '{print \$1}')
            next_index=\$(((\$end_of_file+1)))
        fi
        row_count=\$(((\$next_index-\$start_index)))
        end_index=\$(((\$next_index-1)))
        head -n \$end_index ${assembly.baseName}_nuclear.fna | tail -n \$row_count >> ${assembly.baseName}_mitochondria.fna
        sed -e \$start_index,\$end_index\\d ${assembly.baseName}_nuclear.fna > intermediate_file.fna
        mv intermediate_file.fna ${assembly.baseName}_nuclear.fna
    done
    """
}


workflow.onComplete {
    log.info ( workflow.success ? "\nFunctional annotation input preparation complete!\n" : "Oops .. something went wrong\n" )
}