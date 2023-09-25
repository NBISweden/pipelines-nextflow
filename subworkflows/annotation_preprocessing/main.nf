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
    
    ASSEMBLY_PURIFY( genome_assembly )
    ASSEMBLY_STATS( genome_assembly.mix( ASSEMBLY_PURIFY.out.fasta ) )
    BUSCO( 
        ASSEMBLY_PURIFY.out.fasta.map { fasta -> [ [ id: fasta.baseName ], fasta ] }, 
        params.busco_lineage,
        params.busco_lineages_path ? file( params.busco_lineages_path, checkIfExists: true ) : [],
        []
    )
}
