include { AGAT_SEPARATEBYRECORD                 as SPLIT_MAKER_EVIDENCE          } from "$projectDir/modules/local/agat/separatebyrecord"
include { AGAT_FILTERBYATTRIBUTE                as MODEL_SELECTION_BY_AED        } from "$projectDir/modules/local/agat/filterbyattribute"
include { AGAT_KEEPLONGESTISOFORM               as RETAIN_LONGEST_ISOFORM        } from "$projectDir/modules/local/agat/keeplongestisoform"
include { AGAT_FILTERINCOMPLETEGENECODINGMODELS as REMOVE_INCOMPLETE_GENE_MODELS } from "$projectDir/modules/local/agat/filterincompletegenecodingmodels"
include { AGAT_FILTERBYLOCUSDISTANCE            as FILTER_BY_LOCUS_DISTANCE      } from "$projectDir/modules/local/agat/filterbylocusdistance"
include { AGAT_EXTRACTSEQUENCES                 as EXTRACT_PROTEIN_SEQUENCE      } from "$projectDir/modules/local/agat/extractsequences"
include { BLAST_MAKEBLASTDB                                                      } from "$projectDir/modules/nf-core/blast/makeblastdb/main"
include { BLAST_BLASTP                          as BLAST_RECURSIVE               } from "$projectDir/modules/local/blast/blastp"
include { AGAT_FILTERBYMRNABLASTVALUE           as GFF_FILTER_BY_BLAST           } from "$projectDir/modules/local/agat/filterbymrnablastvalue"
include { AUGUSTUS_GFF2GBK                      as GFF2GBK                       } from "$projectDir/modules/local/augustus/gff2gbk"
include { AUGUSTUS_GBK2AUGUSTUS                 as GBK2AUGUSTUS                  } from "$projectDir/modules/local/augustus/gbk2augustus"
include { AUGUSTUS_TRAINING                                                      } from "$projectDir/modules/local/augustus/training"
include { AGAT_GFF2ZFF                          as CONVERT_GFF2ZFF               } from "$projectDir/modules/local/agat/gff2zff"
include { SNAP_TRAINING                                                          } from "$projectDir/modules/local/snap/training"

workflow ABINITIO_TRAINING {

    main:
    log.info """
        Abintio training dataset workflow
        ===================================================
    """

    Channel.fromPath( params.maker_evidence_gff, checkIfExists: true )
        .map { gff -> [ [ id: gff.baseName ], gff ] }
        .set { gff_annotation }
    Channel.fromPath( params.genome, checkIfExists: true )
        .set{ genome }

    // Make channel for sweep parameters
    Channel.fromList( params.aed_value instanceof List ? params.aed_value : [ params.aed_value ] )
        .filter( Number )
        .set{ ch_aed }
    Channel.fromList( params.locus_distance instanceof List ? params.locus_distance : [ params.locus_distance ] )
        .filter( Number )
        .combine( ch_aed )
        .map { locus_distance, aed -> [ 'aed_value': aed, 'locus_distance': locus_distance ] }
        .set { ch_sweep_parameters }

    SPLIT_MAKER_EVIDENCE( gff_annotation )
    MODEL_SELECTION_BY_AED( ch_sweep_parameters.combine( SPLIT_MAKER_EVIDENCE.out.transcripts ).map { pars, id, gff -> [ id + pars, gff ] } )
    RETAIN_LONGEST_ISOFORM( MODEL_SELECTION_BY_AED.out.selected_models )
    REMOVE_INCOMPLETE_GENE_MODELS( 
        RETAIN_LONGEST_ISOFORM.out.longest_isoform,
        genome.collect()
    )
    FILTER_BY_LOCUS_DISTANCE( REMOVE_INCOMPLETE_GENE_MODELS.out.complete_gene_models )
    EXTRACT_PROTEIN_SEQUENCE( 
        FILTER_BY_LOCUS_DISTANCE.out.distanced_models,
        genome.collect()
    )
    BLAST_MAKEBLASTDB( EXTRACT_PROTEIN_SEQUENCE.out.proteins )
    BLAST_RECURSIVE( 
        EXTRACT_PROTEIN_SEQUENCE.out.proteins
            .join( BLAST_MAKEBLASTDB.out.db )
            .multiMap { meta, proteins, proteindb -> 
                proteins: [ meta, proteins ]
                prot_db: proteindb
            }
    )
    GFF_FILTER_BY_BLAST( 
        FILTER_BY_LOCUS_DISTANCE.out.distanced_models
            .combine( BLAST_RECURSIVE.out.txt, by: 0 )
            .multiMap{ meta, dmodels, btbl ->
                dmodels  : [ meta, dmodels ] 
                blast_tbl: btbl
            }
    )
    GFF2GBK( 
        GFF_FILTER_BY_BLAST.out.blast_filtered,
        genome.collect()
    )
    GBK2AUGUSTUS( GFF2GBK.out.gbk )
    AUGUSTUS_TRAINING( 
        GBK2AUGUSTUS.out.training_data
            .combine( GBK2AUGUSTUS.out.testing_data , by: 0 )
            .multiMap{ meta, training, testing ->
                train: [ meta, training ]
                test : testing
            },
        params.species_label
    )
    CONVERT_GFF2ZFF(
        GFF_FILTER_BY_BLAST.out.blast_filtered,
        genome.collect()
    )
    SNAP_TRAINING(
        CONVERT_GFF2ZFF.out.zff,
        params.species_label
    )
}
