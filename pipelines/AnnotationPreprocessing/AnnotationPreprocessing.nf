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

    get:
        genome_assembly

    main:
        fasta_filter_size(genome_assembly)
        assembly_generate_stats(fasta_filter_size.out)
        busco(fasta_filter_size.out,params.busco_lineage)

}

process fasta_filter_size {

    tag "${fasta_file.baseName} ; min length = ${params.min_length}"
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_min${params.min_length}.fasta"

    script:
    """
    seqtk seq -A $fasta_file -L ${params.min_length} > ${fasta_file.baseName}_min${params.min_length}.fasta
    """

}

process assembly_generate_stats {

    // FIXME: Options:
    // a) Replace with agat script,
    // b) Include script in bin directory,
    // c) Include GAAS as subtree

    tag "${fasta_file.simpleName}"
    publishDir "${params.outdir}/stats", mode: 'copy'
    label 'AGAT'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_assembly_report"

    script:
    """
    fasta_statisticsAndPlot.pl --infile $fasta_file --output ${fasta_file.baseName}_assembly_report
    """
    // fasta_statisticsAndPlot.pl can be found in the NBIS GAAS repository
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
    busco -c ${task.cpus} -i $fasta -l $lineage -m genome --out $out
    """
}

workflow.onComplete {
    log.info ( workflow.success ? "\nAnnotation preprocessing complete!\n" : "Oops .. something went wrong\n" )
}
