#! /usr/bin/env nextflow

nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overridden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.maker_evidence_gff = "/path/to/maker/evidence.gff"
params.genome = "/path/to/genome/assembly.fasta"
params.outdir = "results"
params.species_label = 'test_species'  // e.g. 'asecodes_parviclava'

params.codon_table = 1

params.test_size = 100
params.flank_region_size = 500
params.maker_species_publishdir = '/path/to/shared/maker/folder/' // e.g. '/projects/references/augustus/config/species/'

log.info """
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Abintio training dataset workflow
 ===================================

 General Parameters
     maker_evidence_gff            : ${params.maker_evidence_gff}
     genome                        : ${params.genome}
     outdir                        : ${params.outdir}
     species_label                 : ${params.species_label}

 Protein Sequence extraction parameters
     codon_table                   : ${params.codon_table}

 Augustus training parameters
     test_size                     : ${params.test_size}
     flank_region_size             : ${params.flank_region_size}

 """

workflow {

    main:
        Channel.fromPath(params.maker_evidence_gff, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find gff file matching ${params.maker_evidence_gff}!\n" }
            .set {evidence}
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
            .set{genome}

        abinitio_training(evidence,genome)

}

workflow abinitio_training {

    take:
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
        augustus_training(gbk2augustus.out[0],gbk2augustus.out[1],params.species_label)
        convert_gff2zff(gff_filter_by_blast.out,genome.collect())
        snap_training(convert_gff2zff.out,params.species_label)

}

process split_maker_evidence {

    tag "${maker_evidence.baseName}"
    publishDir "${params.outdir}", mode: 'copy'
    label 'AGAT'

    input:
    path maker_evidence

    output:
    path "maker_results_noAbinitio_clean/mrna.gff"
    path "maker_results_noAbinitio_clean/*"

    script:
    """
    agat_sp_split_by_level2_feature.pl -g ${maker_evidence} -o maker_results_noAbinitio_clean
    if test -f maker_results_noAbinitio_clean/mrna.gff && test -f maker_results_noAbinitio_clean/transcript.gff; then
        agat_sp_merge_annotations.pl --gff maker_results_noAbinitio_clean/mrna.gff \\
            --gff maker_results_noAbinitio_clean/transcript.gff --out merged_transcripts.gff
        mv merged_transcripts.gff maker_results_noAbinitio_clean/mrna.gff
    elif test -f maker_results_noAbinitio_clean/transcript.gff; then
        cp maker_results_noAbinitio_clean/transcript.gff maker_results_noAbinitio_clean/mrna.gff
    fi
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

process augustus_training {

    tag "$species_label"
    label 'Augustus'
    publishDir "${params.outdir}/Augustus_training", mode: 'copy'
    publishDir "${params.maker_species_publishdir}", mode: 'copy', enabled: file(params.maker_species_publishdir).exists(), pattern: "${species_label}"

    input:
    path training_file
    path test_file
    val species_label

    output:
    path "${species_label}_run.log"
    path "${species_label}"

    script:
    """
    cp -rv \${AUGUSTUS_CONFIG_PATH}/ .
    export AUGUSTUS_CONFIG_PATH="\$PWD/config"
    new_species.pl --species=$species_label
    etraining --species=$species_label $training_file
    augustus --species=$species_label $test_file | tee ${species_label}_run.log
    mv config/species/${species_label} .
    """

}

process convert_gff2zff {

    label 'AGAT'

    input:
    path annotation
    path genome

    output:
    path "*.{ann,dna}"

    script:
    """
    agat_convert_sp_gff2zff.pl --gff $annotation \\
        --fasta $genome -o ${genome.baseName}
    """
}

process snap_training {

    input:
    path training_files
    val species_label

    output:
    path "*.hmm"

    script:
    ann_file = training_files.find { it =~ /.ann$/ }
    dna_file = training_files.find { it =~ /.dna$/ }
    """
    fathom -categorize 1000 ${ann_file} ${dna_file}
    fathom -export 1000 -plus uni.ann uni.dna
    forge export.ann export.dna
    hmm-assembler.pl "$species_label" . > "${species_label}.hmm"
    """
}

workflow.onComplete {
    log.info ( workflow.success ? "\nAugustus training dataset complete!\n" : "Oops .. something went wrong\n" )
}
