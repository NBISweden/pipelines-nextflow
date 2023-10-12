include { GAAS_FASTAPURIFY     as ASSEMBLY_PURIFY } from "$projectDir/modules/local/gaas/fastapurify"
include { GAAS_FASTASTATISTICS as ASSEMBLY_STATS  } from "$projectDir/modules/local/gaas/fastastatistics"
include { BUSCO                                   } from "$projectDir/modules/nf-core/busco/main"

workflow ANNOTATION_PREPROCESSING {

    main:
    log.info """
        Annotation preprocessing workflow
        ===================================================
    """
    Channel.fromPath( params.genome, checkIfExists: true)
        .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" }
        .set { genome_assembly }
    Channel.fromList( params.busco_lineage instanceof List ? params.busco_lineage : [ params.busco_lineage ] )
        .set { ch_busco_lineage }
    
    ASSEMBLY_PURIFY( genome_assembly )
    ASSEMBLY_STATS( genome_assembly.mix( ASSEMBLY_PURIFY.out.fasta ) )
    BUSCO( 
        ASSEMBLY_PURIFY.out.fasta
            .combine( ch_busco_lineage )
            .multiMap { fasta, lineage ->
                ch_fasta: [ [ id: fasta.baseName ], fasta ]
                ch_busco: lineage
            },
        params.busco_lineages_path ? file( params.busco_lineages_path, checkIfExists: true ) : [],
        []
    )
}
