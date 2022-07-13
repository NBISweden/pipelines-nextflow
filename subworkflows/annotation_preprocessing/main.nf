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
    ASSEMBLY_GENERATE_STATS( genome_assembly.mix( ASSEMBLY_PURIFY.out ) )
    BUSCO( ASSEMBLY_PURIFY.out, params.busco_lineage )
}

process ASSEMBLY_PURIFY {

    tag "${fasta_file.baseName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::gaas=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gaas:1.2.0--pl526r35_0':
        'quay.io/biocontainers/gaas:1.2.0--pl526r35_0' }"

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_purified/${fasta_file.baseName}_purified.fa", emit: fasta
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 1.2.0  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gaas_fasta_purify.pl --infile $fasta_file --size ${params.min_length} --output ${fasta_file.baseName}_purified

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gaas: $VERSION
    END_VERSIONS
    """
}

process ASSEMBLY_GENERATE_STATS {

    tag "${fasta_file.simpleName}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda (params.enable_conda ? "bioconda::gaas=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gaas:1.2.0--pl526r35_0':
        'quay.io/biocontainers/gaas:1.2.0--pl526r35_0' }"

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_report", emit: report
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 1.2.0  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gaas_fasta_statistics.pl --infile $fasta_file --output ${fasta_file.baseName}_report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gaas: $VERSION
    END_VERSIONS
    """
}

process BUSCO {

    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::busco=5.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.2--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.3.2--pyhdfd78af_0' }"

    input:
    path fasta
    each lineage

    output:
    path out, emit: busco
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out = "busco_${fasta.baseName}_${lineage}"
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        # . "/usr/local/env-activate.sh"  # Errors out because of various unbound variables
        export PATH='/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin'
        export CONDA_PREFIX='/usr/local'
        export CONDA_SHLVL='1'
        export CONDA_DEFAULT_ENV='/usr/local'
        export CONDA_PROMPT_MODIFIER=''
        . "/usr/local/etc/conda/activate.d/activate-r-base.sh"
        . "/usr/local/etc/conda/activate.d/augustus.sh"
        . "/usr/local/etc/conda/activate.d/openjdk_activate.sh"
    fi
    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi
    busco -c ${task.cpus} -i $fasta -l $lineage -m genome --out $out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
