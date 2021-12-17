process MAPPING {

    conda "${task.ext.enable_conda ? 'bioconda::minimap2=2.18 bioconda::samtools=1.12' : '' }"
    container "${workflow.containerEngine == 'singularity' &&
                  !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:fd08f3e36cad38244f4a20cb4aebe7d795774ee9-0' :
              'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:fd08f3e36cad38244f4a20cb4aebe7d795774ee9-0' }"

    scratch = true

    input:
    path assembly
    path reads_file

    output:
    path "depth_full_alignment.tsv", emit: depth_file

    script:
    """
    minimap2 -d assembly.mmi $assembly
    minimap2 -a assembly.mmi $reads_file | samtools sort --threads $task.cpus -o sorted_assembly_full_alignment.bam
    samtools depth -o depth_full_alignment.tsv sorted_assembly_full_alignment.bam
    """
}