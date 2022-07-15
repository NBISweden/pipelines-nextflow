# Functional annotation pipeline

## Quick start

Run workflow using the singularity profile

`params.yml`:

```yml
subworkflow: 'functional_annotation'
genome: '/path/to/genome/assembly.fasta'
gff_annotation: '/path/to/annotation.gff3'
blast_db_fasta: '/path/to/protein/database.fasta'
```

Command line:

```bash
nextflow run main.nf \
    -profile singularity \
    -params-file params.yml
```

> **note**
>
> The Interproscan conda package does not appear to be fully functional. Please use a local installation
> by overriding the workflow configuration.
>
> `nextflow.config`:
>
> ```nextflow
> process {
>     withName: 'INTERPROSCAN' {
>         conda = null
>         container = null
>         module = 'bioinfo-tools:InterProScan/5.30-69.0'
>     }    
> }
> ```

## Parameters

- General:
  - `gff_annotation`:  Path to GFF genome annotation.
  - `genome`: Path to the genome assembly.
  - `outdir`: Path to the results folder.
  - `records_per_file`: Number of fasta records per file to distribute to blast and interproscan (default: 1000).
  - `codon_table`: (default: 1).

## Workflow Stages

1. Extract protein sequences based on GFF coordinates.
2. Blast protein sequences against protein database.
3. Query protein sequences against interproscan databases.
4. Merge functional annotations.
