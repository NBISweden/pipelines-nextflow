process TBLASTN {

    conda "${task.ext.enable_conda ? 'bioconda::blast:2.12.0' : '' }"
    container "${workflow.containerEngine == 'singularity' &&
                  !task.ext.singularity_pull_docker_container ?
              'https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0' :
              'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0' }"

    label 'blast'

    input:
    path reference_organelle
    path blastdb
    val blast_evalue

    output:
    path "output_blast.tsv", emit: output_blast

    script:
    database = blastdb.find { it =~ /\.f(ast|n)?a$/ }
    """
    tblastn -query $reference_organelle -db ${database} -evalue $blast_evalue -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -out output_blast.tsv
    """

}
