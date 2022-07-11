#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

// Abinitio training modules

// Annotation preprocessing modules
include { BUSCO } from "$projectDir/modules/nf-core/modules/busco/main"

// Functional annotation modules

// Transcript assembly modules
include { FASTQC              } from "$projectDir/modules/nf-core/modules/fastqc/main"
include { HISAT2_ALIGN        } from "$projectDir/modules/nf-core/modules/hisat2/align/main"
include { HISAT2_BUILD        } from "$projectDir/modules/nf-core/modules/hisat2/build/main"
include { FASTP               } from "$projectDir/modules/nf-core/modules/fastp/main"
include { STRINGTIE_STRINGTIE } from "$projectDir/modules/nf-core/modules/stringtie/stringtie/main"
include { MULTIQC             } from "$projectDir/modules/nf-core/modules/multiqc/main"

workflow {

    log.info """
         _  _ ___ ___ ___ 
        | \| | _ )_ _/ __|
        | .` | _ \| |\__ \
        |_|\_|___/___|___/ Annotation Service

    """

    def valid_subworkflows = [ 'abinitio_training', 'annotation_preprocessing', 'functional_annotation', 'transcript_assembly' ]
    if( ! params.subworkflow in valid_subworkflows ){
        error """
        The parameter 'subworkflow' (${params.subworkflow}) is not a valid subworkflow.
        
        Please select a valid subworkflow ([ ${valid_subworkflows.join(', ')} ]).
        """
    }

    // Ab initio training subworkflow
    if ( params.subworkflow == 'abinitio_training' ){
        log.info """
            Abintio training dataset workflow
            ===================================================
        """

        Channel.fromPath(params.maker_evidence_gff, checkIfExists: true)
            .ifEmpty { error "Cannot find gff file matching ${params.maker_evidence_gff}!\n" }
            .set { evidence }
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" }
            .set{ genome }

        split_maker_evidence(gff_annotation)
        model_selection_by_AED(split_maker_evidence.out.transcripts)
        retain_longest_isoform(model_selection_by_AED.out.selected_models)
        remove_incomplete_gene_models(retain_longest_isoform.out.longest_isoform,
            genome.collect())
        filter_by_locus_distance(remove_incomplete_gene_models.out.complete_gene_models)
        extract_protein_sequence(filter_by_locus_distance.out.distanced_models,
            genome.collect())
        blast_makeblastdb(extract_protein_sequence.out.proteins)
        blast_recursive(extract_protein_sequence.out.proteins,
            blast_makeblastdb.out.collect())
        gff_filter_by_blast(filter_by_locus_distance.out.distanced_models,
            blast_recursive.out.collect())
        gff2gbk(gff_filter_by_blast.out.blast_filtered,genome.collect())
        gbk2augustus(gff2gbk.out)
        augustus_training(gbk2augustus.out.training_data,gbk2augustus.out.testing_data,params.species_label)
        convert_gff2zff(gff_filter_by_blast.out.blast_filtered,genome.collect())
        snap_training(convert_gff2zff.out,params.species_label)
    }

    // Annotation preprocessing subworkflow
    if ( params.subworkflow == 'annotation_preprocessing' ){
        log.info """
            Annotation preprocessing workflow
            ===================================================
        """
        Channel.fromPath( params.genome, checkIfExists: true)
            .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" } \
            .set { genome_assembly }
        assembly_purify(genome_assembly)
        assembly_generate_stats(genome_assembly.mix(assembly_purify.out))
        BUSCO( assembly_purify.out, params.busco_lineage )
    }

    // Functional annotation subworkflow
    if ( params.subworkflow == 'functional_annotation' ){
        log.info """
            Functional annotation workflow
            ===================================================
        """
        Channel.fromPath(params.gff_annotation, checkIfExists: true)
            .ifEmpty { error "Cannot find gff file matching ${params.gff_annotation}!\n" }
            .set {annotation}
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" }
            .set {genome}
        Channel.fromPath("${params.blast_db_fasta}" + (params.blast_db_fasta =~ /^(ht|f)tps?:/ ? '': "*"), checkIfExists: true)
            .ifEmpty { error "Cannot find blast database files matching ${params.blast_db_fasta}?(.p*)" }
            .set {blastdb}

        makeblastdb(blastdb,blastdb.filter { it =~ /.p(hr|in|sq)$/ }.ifEmpty('DBFILES_ABSENT'))
        makeblastdb.out.mix(blastdb.filter { !(it =~ /[^.]f(ast|n|)a$/) }).unique().collect().set{blastfiles}
        gff2protein(gff_file,genome.collect())
        blastp(gff2protein.out.splitFasta(by: params.records_per_file, file: true),
            blastfiles)
        interproscan(gff2protein.out.splitFasta(by: params.records_per_file, file: true))
        merge_functional_annotation(gff_file,
            blastp.out.collectFile(name:'blast_merged.tsv').collect(),
            interproscan.out.collectFile(name:'interproscan_merged.tsv').collect(),
            blastdb.collect())
    }
    
    // Transcript assembly subworkflow
    if ( params.subworkflow == 'transcript_assembly' ){
        log.info """
            Transcript assembly using Hisat2/Stringtie workflow
            ===================================================
        """
        Channel.fromFilePairs( params.reads, size: params.single_end ? 1 : 2, checkIfExists: true)
            // .ifEmpty { error "Cannot find reads matching ${params.reads}!\n" }
            .set { reads }
        Channel.fromPath( params.genome, checkIfExists: true)
            // .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" }
            .set { genome }

        FASTQC ( reads )
        HISAT2_BUILD ( genome )
        if ( !params.skip_trimming ){
                FASTP( reads )
                FASTP.out.reads.set { trimmed_reads }
                FASTP.out.log.set { trimming_logs }
        } else {
            reads.set { trimmed_reads }
        }
        HISAT2_ALIGN ( trimmed_reads, hisat2_index.out.collect() )
        STRINGTIE_STRINGTIE ( hisat2.out[0] )
        FASTQC.out.mix( trimming_logs ).mix( hisat2.out[2] ).set {logs}
        MULTIQC( logs.collect(),params.multiqc_config )
    }
}

workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Workflow completed successfully. 

        Thank you for using our workflow.
        Results are located in the folder: $params.outdir
        """)
    } else {
        log.info("""
        The workflow completed unsuccessfully.
        Please read over the error message. If you are unable to solve it, please
        post an issue at https://github.com/NBISweden/pipelines-nextflow/issues
        where we will do our best to help.
        """)
    }
}