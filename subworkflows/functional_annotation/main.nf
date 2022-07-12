workflow FUNCTIONAL_ANNOTATION {

    main:
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