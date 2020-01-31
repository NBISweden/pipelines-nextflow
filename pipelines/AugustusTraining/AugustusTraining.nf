#! /usr/bin/env nextflow

nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overridden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.maker_evidence_gff = "/path/to/maker/evidence.gff"
params.genome = "/path/to/genome/assembly.fasta"
params.outdir = "results"

params.gff_gene_model_filter_options = '-c -r -d 500 -a 0.3'

params.codon_table = 1

params.test_size = 100
params.flank_region_size = 500

log.info """
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Augustus training dataset workflow
 ===================================

 General Parameters
     maker_evidence_gff            : ${params.maker_evidence_gff}
     genome                        : ${params.genome}
     outdir                        : ${params.outdir}

 Gene Model Filtering parameters
     gff_gene_model_filter_options : ${params.gff_gene_model_filter_options}

 Protein Sequence extraction parameters
     codon_table                   : ${params.codon_table}

 Augustus training parameters
     test_size                     : ${params.test_size}
     flank_region_size             : ${params.flank_region_size}

 """

workflow {

    main:
        evidence = Channel.fromPath(params.maker_evidence_gff, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find gff file matching ${params.maker_evidence_gff}!\n" }
        genome = Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }

        augustus_training_dataset(evidence,genome)

}

workflow augustus_training_dataset {

    get:
        gff_annotation
        genome

    main:
        split_maker_evidence(gff_annotation)
        model_selection_by_AED(split_maker_evidence.out[0])
        retain_longest_isoform(model_selection_by_AED.out)
        remove_incomplete_gene_models(retain_longest_isoform.out,
            genome.collect())
        filter_by_locus_distance(remove_incomplete_gene_models.out)
        extract_protein_sequence(filter_by_locus_distance.out,
            genome.collect())
        blast_makeblastdb(extract_protein_sequence.out)
        blast_recursive(extract_protein_sequence.out,
            blast_makeblastdb.out.collect())
        gff_filter_by_blast(filter_by_locus_distance.out,
            blast_recursive.out.collect())
        gff2gbk(gff_filter_by_blast.out,genome.collect())
        gbk2augustus(gff2gbk.out)

}

process split_maker_evidence {

    tag "${maker_evidence.baseName}"
    publishDir "${params.outdir}/maker_results_noAbinitio_clean", mode: 'copy'
    label 'AGAT'

    input:
    path maker_evidence

    output:
    path "maker_results_noAbinitio_clean/mrna.gff"
    path "maker_results_noAbinitio_clean/*"

    script:
    """
    agat_sp_split_by_level2_feature.pl -g ${maker_evidence} -o maker_results_noAbinitio_clean
    """
    // agat_sp_split_by_level2_feature.pl is a script from AGAT
}

process model_selection_by_AED {

    tag "${mrna_gff.baseName}"
    publishDir "${params.outdir}/filter", mode: 'copy'
    label 'AGAT'

    input:
    path mrna_gff

    output:
    path "codingGeneFeatures.filter.gff"

    script:
    """
    agat_sp_filter_feature_by_attribute_value.pl --gff ${mrna_gff} --value 0.3 -a _AED -t ">=" -o codingGeneFeatures.filter.gff
    """
    // agat_sp_filter_feature_by_attribute_value.pl is a script from AGAT
}

process retain_longest_isoform {

    tag "${coding_gene_features_gff.baseName}"
    publishDir "${params.outdir}/filter", mode: 'copy'
    label 'AGAT'

    input:
    path coding_gene_features_gff

    output:
    path "codingGeneFeatures.filter.longest_cds.gff"

    script:
    """
    agat_sp_keep_longest_isoform.pl -f ${coding_gene_features_gff} -o codingGeneFeatures.filter.longest_cds.gff
    """
    // agat_sp_keep_longest_isoform.pl is a script from AGAT
}

process remove_incomplete_gene_models {

    tag "${coding_gene_features_gff.baseName}"
    publishDir "${params.outdir}/filter", mode: 'copy'
    label 'AGAT'

    input:
    path coding_gene_features_gff
    path genome_fasta

    output:
    path "codingGeneFeatures.filter.longest_cds.complete.gff"

    script:
    """
    agat_sp_filter_incomplete_gene_coding_models.pl --gff ${coding_gene_features_gff} \
        -f ${genome_fasta} -o codingGeneFeatures.filter.longest_cds.complete.gff
    """
    // agat_sp_filter_incomplete_gene_coding_models.pl is a script from AGAT
}

process filter_by_locus_distance {

    tag "${coding_gene_features_gff.baseName}"
    publishDir "${params.outdir}/filter", mode: 'copy'
    label 'AGAT'

    input:
    path coding_gene_features_gff

    output:
    path "codingGeneFeatures.filter.longest_cds.complete.good_distance.gff"

    script:
    """
    agat_sp_filter_by_locus_distance.pl --gff ${coding_gene_features_gff} -o codingGeneFeatures.filter.longest_cds.complete.good_distance.gff
    """
    // agat_sp_filter_by_locus_distance.pl is a script from AGAT
}

process extract_protein_sequence {

    tag "${gff_file.baseName}"
    label 'AGAT'

    input:
    path gff_file
    path genome_fasta

    output:
    path "${gff_file.baseName}_proteins.fasta"

    script:
    """
    agat_sp_extract_sequences.pl -o ${gff_file.baseName}_proteins.fasta -f $genome_fasta \\
        -p -cfs -cis -ct ${params.codon_table} --g $gff_file
    """
    // agat_sp_extract_sequences.pl is a script from AGAT

}

process blast_makeblastdb {

    tag "${fasta_file.baseName} type: $dbtype"
    label 'Blast'

    input:
    path fasta_file

    output:
    path "*.{phr,pin,psq}"

    script:
    """
    makeblastdb -in $fasta_file -dbtype prot
    """

}

process blast_recursive {

    tag "${fasta_file.baseName}"
    label 'Blast'

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

process gff_filter_by_blast {

    tag "${gff_file.baseName}"
    publishDir "${params.outdir}/BlastFilteredGFF", mode: 'copy'
    label 'AGAT'

    input:
    path gff_file
    path blast_file

    output:
    path "${gff_file.baseName}_blast-filtered.gff3"

    script:
    """
    agat_sp_filter_by_mrnaBlastValue.pl --gff $gff_file --blast $blast_file \\
        --outfile ${gff_file.baseName}_blast-filtered.gff3
    """
    // agat_sp_filter_by_mrnaBlastValue.pl is a script from AGAT

}

process gff2gbk {

    tag "${gff_file.baseName}"
    label 'Augustus'

    input:
    path gff_file
    path genome_fasta

    output:
    path "${gff_file.baseName}.gbk"

    script:
    """
    gff2gbSmallDNA.pl $gff_file $genome_fasta ${params.flank_region_size} ${gff_file.baseName}.gbk
    """
    // gff2gbSmallDNA.pl is a script in the Augustus package

}

process gbk2augustus {

    tag "Make Augustus training set: ${genbank_file.baseName}"
    label 'Augustus'
    publishDir "${params.outdir}/Augustus", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf(".train") > 0)        "TrainingData/$filename"
            else if (filename.indexOf(".test") > 0)    "TestingData/$filename"
            else if (filename.indexOf(".gbk") > 0)     "GenbankFile/$filename"
            else filename }

    input:
    path genbank_file

    output:
    path "${genbank_file}.train"
    path "${genbank_file}.test"
    path "${genbank_file}"

    script:
    """
    randomSplit.pl $genbank_file ${params.test_size}
    """
    // randomSplit.pl is a script in the Augustus package

}

workflow.onComplete {
    log.info ( workflow.success ? "\nAugustus training dataset complete!\n" : "Oops .. something went wrong\n" )
}
