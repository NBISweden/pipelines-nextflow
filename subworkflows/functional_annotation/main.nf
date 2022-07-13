workflow FUNCTIONAL_ANNOTATION {

    main:
    log.info """
        Functional annotation workflow
        ===================================================
    """
    Channel.fromPath( params.gff_annotation, checkIfExists: true )
        // .ifEmpty { error "Cannot find gff file matching ${params.gff_annotation}!\n" }
        .set { gff_file }
    Channel.fromPath( params.genome, checkIfExists: true )
        // .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" }
        .set { genome }
    // Channel.fromPath("${params.blast_db_fasta}" + (params.blast_db_fasta =~ /^(ht|f)tps?:/ ? '': "*"), checkIfExists: true )
    //     .ifEmpty { error "Cannot find blast database files matching ${params.blast_db_fasta}?(.p*)" }
    //     .set { blastdb }
    Channel.fromPath( params.blast_db_fasta, checkIfExists: true )
        .branch { fasta -> 
            def db_files = [ fasta ] 
            try {
                db_files = [ fasta ] + file( fasta + ".p*", checkIfExists: true )
            } catch ( Exception e ){
                // No database files found matching the glob pattern
            }
            make_db : db_files.size() == 1
                return db_files
            with_db : db_files.size() > 1
                return db_files
        }.set { blast_fa }

    MAKEBLASTDB(
        blast_fa.make_db
    )
    blastdb_ch = MAKEBLASTDB.out.db.mix( blast_fa.with_db )
    GFF2PROTEIN( 
        gff_file, 
        genome.collect()
    )
    BLASTP(
        GFF2PROTEIN.out.splitFasta( by: params.records_per_file, file: true ),
        blastdb_ch.collect()
    )
    INTERPROSCAN( GFF2PROTEIN.out.splitFasta( by: params.records_per_file, file: true ) )
    MERGE_FUNCTIONAL_ANNOTATION( 
        gff_file,
        BLASTP.out.collectFile( name: 'blast_merged.tsv' ).collect(),
        INTERPROSCAN.out.collectFile( name: 'interproscan_merged.tsv' ).collect(),
        blastdb_ch.collect()
    )
}

process GFF2PROTEIN {

    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path gff_file
    path genome_fasta

    output:
    path "${gff_file.baseName}_proteins.fasta", emit: proteins
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 0.9.2  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_extract_sequences.pl -o ${gff_file.baseName}_proteins.fasta -f $genome_fasta \\
        -p -cfs -cis -ct ${params.codon_table} --gff $gff_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}

process MAKEBLASTDB {

    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
        'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    input:
    path fasta

    output:
    path "*.fasta*", emit: db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ "$fasta" =~ \\.f(ast|n)?a\$ ]]; then
        makeblastdb -in $fasta -dbtype prot
    else
        cp $fasta protein.fasta
        makeblastdb -in protein.fasta -dbtype prot
    fi
    """
}

process BLASTP {

    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
        'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    input:
    path fasta_file
    path blastdb

    output:
    path "${fasta_file.baseName}_blast.tsv"
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // database = blastdb[0].toString() - ~/.p\w\w$/
    database = blastdb.find { it =~ /\.f(ast|n)?a$/ }
    """
    blastp -query $fasta_file -db ${database} -num_threads ${task.cpus} \\
        -evalue ${params.blast_evalue} -outfmt 6 -out ${fasta_file.baseName}_blast.tsv
    """
}

process INTERPROSCAN {

    // FIXME:Source from path if not in container?

    input:
    path protein_fasta

    output:
    path '*.tsv'
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    applications = params.interproscan_db ? "-appl ${params.interproscan_db}" : ''
    tmpdir = task.scratch ? "-T ${task.scratch}" : ''
  
    """
    export PATH="/opt/interproscan:\$PATH"
    interproscan.sh ${applications} -i $protein_fasta -o ${protein_fasta.baseName}.tsv \\
        -f TSV --iprlookup --goterms -pa -dp -t p ${tmpdir}
    """
}

process MERGE_FUNCTIONAL_ANNOTATION {

    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path gff_annotation
    path merged_blast_results
    path merged_interproscan_results
    path blast_files

    output:
    path "${gff_annotation.baseName}_plus-functional-annotation.gff"
    path "*.tsv", includeInputs:true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    fasta = blast_files.find { it =~ /\.f(ast|n)?a$/ }
    def VERSION = 0.9.2  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_manage_functional_annotation.pl -f ${gff_annotation} \\
        -b ${merged_blast_results} -i ${merged_interproscan_results} \\
        -db ${params.blast_db_fasta} -id ${params.id_prefix} \\
        -pe ${params.protein_existence} \\
        -o ${gff_annotation.baseName}_plus-functional-annotation.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}
