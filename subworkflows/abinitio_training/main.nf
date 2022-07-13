workflow ABINITIO_TRAINING {

    main:
    log.info """
        Abintio training dataset workflow
        ===================================================
    """

    Channel.fromPath( params.maker_evidence_gff, checkIfExists: true )
        // .ifEmpty { error "Cannot find gff file matching ${params.maker_evidence_gff}!\n" }
        .set { gff_annotation }
    Channel.fromPath( params.genome, checkIfExists: true )
        // .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" }
        .set{ genome }

    SPLIT_MAKER_EVIDENCE( gff_annotation )
    MODEL_SELECTION_BY_AED( SPLIT_MAKER_EVIDENCE.out.transcripts )
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
        EXTRACT_PROTEIN_SEQUENCE.out.proteins,
        BLAST_MAKEBLASTDB.out.collect()
    )
    GFF_FILTER_BY_BLAST( 
        FILTER_BY_LOCUS_DISTANCE.out.distanced_models,
        BLAST_RECURSIVE.out.collect()
    )
    GFF2GBK( 
        GFF_FILTER_BY_BLAST.out.blast_filtered,
        genome.collect()
    )
    GBK2AUGUSTUS( GFF2GBK.out )
    AUGUSTUS_TRAINING( 
        GBK2AUGUSTUS.out.training_data,
        GBK2AUGUSTUS.out.testing_data,
        params.species_label
    )
    CONVERT_GFF2ZFF(
        GFF_FILTER_BY_BLAST.out.blast_filtered,
        genome.collect()
    )
    SNAP_TRAINING(
        CONVERT_GFF2ZFF.out,
        params.species_label
    )

}

process SPLIT_MAKER_EVIDENCE {

    tag "${maker_evidence.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path maker_evidence

    output:
    path "maker_results_noAbinitio_clean/mrna.gff", emit: transcripts
    path "maker_results_noAbinitio_clean/*", emit: all // FIXME: check
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 0.9.2  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_separate_by_record_type.pl -g ${maker_evidence} -o maker_results_noAbinitio_clean
    if test -f maker_results_noAbinitio_clean/mrna.gff && test -f maker_results_noAbinitio_clean/transcript.gff; then
        agat_sp_merge_annotations.pl --gff maker_results_noAbinitio_clean/mrna.gff \\
            --gff maker_results_noAbinitio_clean/transcript.gff --out merged_transcripts.gff
        mv merged_transcripts.gff maker_results_noAbinitio_clean/mrna.gff
    elif test -f maker_results_noAbinitio_clean/transcript.gff; then
        cp maker_results_noAbinitio_clean/transcript.gff maker_results_noAbinitio_clean/mrna.gff
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}

process MODEL_SELECTION_BY_AED {

    tag "${mrna_gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path mrna_gff

    output:
    path "codingGeneFeatures.filter.gff", emit: selected_models
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 0.9.2  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_filter_feature_by_attribute_value.pl --gff ${mrna_gff} --value ${params.model_selection_value} -a _AED -t ">" -o codingGeneFeatures.filter.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}

process RETAIN_LONGEST_ISOFORM {

    tag "${coding_gene_features_gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path coding_gene_features_gff

    output:
    path "codingGeneFeatures.filter.longest_cds.gff", emit: longest_isoform
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 0.9.2  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_keep_longest_isoform.pl -f ${coding_gene_features_gff} -o codingGeneFeatures.filter.longest_cds.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}

process REMOVE_INCOMPLETE_GENE_MODELS {

    tag "${coding_gene_features_gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path coding_gene_features_gff
    path genome_fasta

    output:
    path "codingGeneFeatures.filter.longest_cds.complete.gff", emit: complete_gene_models
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 0.9.2  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_filter_incomplete_gene_coding_models.pl --gff ${coding_gene_features_gff} \
        -f ${genome_fasta} -o codingGeneFeatures.filter.longest_cds.complete.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}

process FILTER_BY_LOCUS_DISTANCE {

    tag "${coding_gene_features_gff.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path coding_gene_features_gff

    output:
    path "codingGeneFeatures.filter.longest_cds.complete.good_distance.gff", emit: distanced_models
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 0.9.2  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_filter_by_locus_distance.pl --gff ${coding_gene_features_gff} -d ${params.locus_distance} -o codingGeneFeatures.filter.longest_cds.complete.good_distance.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}

process EXTRACT_PROTEIN_SEQUENCE {

    tag "${gff_file.baseName}"
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
        -p -cfs -cis -ct ${params.codon_table} --g $gff_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}

process BLAST_MAKEBLASTDB {

    tag "${fasta_file.baseName} type: $dbtype"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
        'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    input:
    path fasta_file

    output:
    path "*.{phr,pin,psq}", emit: db  // FIXME: Include inputs
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    makeblastdb -in $fasta_file -dbtype prot
    """

}

process BLAST_RECURSIVE {

    tag "${fasta_file.baseName}"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::blast=2.12.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
        'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    input:
    path fasta_file
    path blastdb

    output:
    path "${fasta_file.baseName}_blast.tsv", emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    database = blastdb[0].toString() - ~/.p\w\w$/
    """
    blastp -query $fasta_file -db ${database} -num_threads ${task.cpus} \\
        -outfmt 6 -out ${fasta_file.baseName}_blast.tsv
    """

}

process GFF_FILTER_BY_BLAST {

    tag "${gff_file.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path gff_file
    path blast_file

    output:
    path "${gff_file.baseName}_blast-filtered.gff3", emit: blast_filtered
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 0.9.2  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_sp_filter_by_mrnaBlastValue.pl --gff $gff_file --blast $blast_file \\
        --outfile ${gff_file.baseName}_blast-filtered.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}

process GFF2GBK {

    tag "${gff_file.baseName}"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::augustus=3.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augustus:3.4.0--pl5321h5f9f3d9_6':
        'quay.io/biocontainers/augustus:3.4.0--pl5321h5f9f3d9_6' }"

    input:
    path gff_file
    path genome_fasta

    output:
    path "${gff_file.baseName}.gbk", gbk
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gff2gbSmallDNA.pl $gff_file $genome_fasta ${params.flank_region_size} ${gff_file.baseName}.gbk
    """
    // gff2gbSmallDNA.pl is a script in the Augustus package

}

process GBK2AUGUSTUS {

    tag "Make Augustus training set: ${genbank_file.baseName}"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::augustus=3.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augustus:3.4.0--pl5321h5f9f3d9_6':
        'quay.io/biocontainers/augustus:3.4.0--pl5321h5f9f3d9_6' }"

    input:
    path genbank_file

    output:
    path "${genbank_file}.train", emit: training_data
    path "${genbank_file}.test", emit: testing_data
    path "${genbank_file}", emit: genbank_file
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    randomSplit.pl $genbank_file ${params.test_size}
    """
    // randomSplit.pl is a script in the Augustus package

}

process AUGUSTUS_TRAINING {

    tag "$species_label"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::augustus=3.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/augustus:3.4.0--pl5321h5f9f3d9_6':
        'quay.io/biocontainers/augustus:3.4.0--pl5321h5f9f3d9_6' }"

    input:
    path training_file
    path test_file
    val species_label

    output:
    path "${species_label}_run.log"
    path "${species_label}", emit: training_model
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    : \${AUGUSTUS_CONFIG_PATH:=/usr/local/config}
    cp -rv \${AUGUSTUS_CONFIG_PATH} .
    export AUGUSTUS_CONFIG_PATH="\$PWD/config"
    new_species.pl --species=$species_label
    etraining --species=$species_label $training_file
    augustus --species=$species_label $test_file | tee ${species_label}_run.log
    mv config/species/${species_label} .
    """
}

process CONVERT_GFF2ZFF {

    label 'process_single'

    conda (params.enable_conda ? "bioconda::agat=0.9.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:0.9.2--pl5321hdfd78af_1':
        'quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1' }"

    input:
    path annotation
    path genome

    output:
    path "*.{ann,dna}", emit: zff
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 0.9.2  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    agat_convert_sp_gff2zff.pl --gff $annotation \\
        --fasta $genome -o ${genome.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: $VERSION
    END_VERSIONS
    """
}

process SNAP_TRAINING {

    label 'process_single'

    conda (params.enable_conda ? "bioconda::snap=2013_11_29" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snap:2013_11_29--hec16e2b_4':
        'quay.io/biocontainers/snap:2013_11_29--hec16e2b_4' }"

    input:
    path training_files
    val species_label

    output:
    path "*.hmm", emit: training_model
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    ann_file = training_files.find { it =~ /.ann$/ }
    dna_file = training_files.find { it =~ /.dna$/ }
    """
    fathom -categorize ${params.flank_region_size} ${ann_file} ${dna_file}
    fathom -export ${params.flank_region_size} -plus uni.ann uni.dna
    forge export.ann export.dna
    hmm-assembler.pl "$species_label" . > "${species_label}.hmm"
    """
}
