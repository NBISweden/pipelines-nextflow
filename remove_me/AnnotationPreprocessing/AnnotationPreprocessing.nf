nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.genome = "/path/to/genome/assembly.fasta"
params.outdir = "results"

params.min_length = 1000

params.busco_lineage = [ 'eukaryota_odb10', 'bacteria_odb10' ]

log.info """
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Annotation preprocessing workflow
 ===================================

 General parameters
     genome          : ${params.genome}
     outdir          : ${params.outdir}

 Filtering parameters
     min_length      : ${params.min_length}

 Busco parameters
     busco_lineage      : ${params.busco_lineage}

 """

workflow {

    main:
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" } \
        | annotation_preprocessing

}

workflow annotation_preprocessing {

    take:
        genome_assembly

    main:
        assembly_purify(genome_assembly)
        assembly_generate_stats(genome_assembly.mix(assembly_purify.out))
        busco(assembly_purify.out,params.busco_lineage)

    emit:
        fasta = assembly_purify.out

}

process assembly_purify {

    tag "${fasta_file.baseName} ; min length = ${params.min_length}"
    publishDir "${params.outdir}/assembly", mode: 'copy'
    label 'GAAS'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_purified/${fasta_file.baseName}_purified.fa"

    script:
    """
    gaas_fasta_purify.pl --infile $fasta_file --size ${params.min_length} --output ${fasta_file.baseName}_purified
    """
    // gaas_fasta_purify.pl can be found in the NBIS GAAS repository

}

process assembly_generate_stats {

    tag "${fasta_file.simpleName}"
    publishDir "${params.outdir}/stats", mode: 'copy'
    label 'GAAS'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_report"

    script:
    """
    gaas_fasta_statistics.pl --infile $fasta_file --output ${fasta_file.baseName}_report
    """
    // gaas_fasta_statistics.pl can be found in the NBIS GAAS repository
}

process busco {

    tag "$fasta"
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    path fasta
    each lineage

    output:
    path out

    script:
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

workflow.onComplete {
    log.info ( workflow.success ? "\nAnnotation preprocessing complete!\n" : "Oops .. something went wrong\n" )
}
