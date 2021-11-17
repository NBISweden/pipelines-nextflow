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

// filter parameters
params.bitscore = 200

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

 Filter parameters
     params.bitscore                : ${params.bitscore}
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
        filter(blastx.out.collect())
        extract(genome,
            filter.out.collect())
    emit:
        blast_result = extract.out[0]
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

process filter {

    input:
    path blast_file

    output:
    path "${blast_file.baseName}_filtered.tsv"

    script:
    """
    awk '\$12>${params.bitscore} {print}' $blast_file | awk '{print \$1}' |sort|uniq -c | awk '{print \$2}'> ${blast_file.baseName}_filtered.tsv
    """    
}

process extract {

    publishDir "${params.outdir}", mode: 'copy', pattern: "${assembly.baseName}_mitochondria.fna"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${assembly.baseName}_nuclear.fna"

    input:
    path assembly
    path result_filtered

    output:
    path "${assembly.baseName}_mitochondria.fna"
    path "${assembly.baseName}_nuclear.fna"

    script: 
    """
    cp $assembly ${assembly.baseName}_nuclear.fna
    LINES=\$(cat $result_filtered)
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
