workflow ABINITIO_TRAINING {

    main:
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