workflow ANNOTATION_PREPROCESSING {

    main:
    log.info """
        Annotation preprocessing workflow
        ===================================================
    """
    Channel.fromPath( params.genome, checkIfExists: true)
        .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" }
        .set { genome_assembly }
    assembly_purify(genome_assembly)
    assembly_generate_stats(genome_assembly.mix(assembly_purify.out))
    BUSCO( assembly_purify.out, params.busco_lineage )

}