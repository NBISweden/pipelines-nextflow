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
    label 'GAAS'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_purified/${fasta_file.baseName}_purified.fa"

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gaas_fasta_purify.pl --infile $fasta_file --size ${params.min_length} --output ${fasta_file.baseName}_purified
    """
    // gaas_fasta_purify.pl can be found in the NBIS GAAS repository
}

process ASSEMBLY_GENERATE_STATS {

    tag "${fasta_file.simpleName}"
    label 'GAAS'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_report"

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gaas_fasta_statistics.pl --infile $fasta_file --output ${fasta_file.baseName}_report
    """
    // gaas_fasta_statistics.pl can be found in the NBIS GAAS repository
}

process BUSCO {

    tag "$fasta"

    input:
    path fasta
    each lineage

    output:
    path out

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
    """
}
