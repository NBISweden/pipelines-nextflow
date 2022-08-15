# Functional annotation pipeline

The functional annotation workflow takes a draft assembly (parameter: `genome`) and
predicted gene coordinates (e.g., from Maker; parameter: `gff_annotation`), and assigns functional
annotation based on similarity to existing protein databases (parameter: `blast_db_fasta`).

## Quick start

Run workflow using the singularity profile:

`params.yml`:

```yml
subworkflow: 'functional_annotation'
genome: '/path/to/genome/assembly.fasta'
gff_annotation: '/path/to/annotation.gff3'
blast_db_fasta: '/path/to/protein/database.fasta'
outdir: '/path/to/save/results'
```

Command line:

```bash
nextflow run NBISweden/pipelines-nextflow \
    -profile singularity \
    -params-file params.yml
```

> **note**
>
> The Interproscan conda package is temperamental. Please use a local installation
> by overriding the workflow configuration.
>
> `nextflow.config`:
>
> ```nextflow
> process {
>     withName: 'INTERPROSCAN' {
>         conda = null
>         container = null
>         module = 'bioinfo-tools:InterProScan/5.30-69.0'  // Load Uppmax modules `bioinfo-tools` and `InterProScan/5.30-69.0`
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

### Tool specific parameters

Process specific options are passed by overriding the `ext.args` variable using a process selector in a configuration file.

`nextflow.config`:

```nextflow
process {
    withName: 'INTERPROSCAN' {
        ext.args = '-f TSV --iprlookup --goterms -pa -dp -t p'
    }
}
```

See [Functional annotation modules config](../../config/functional_annotation_modules.config) for the default tool configuration.

## Workflow Stages

1. Extract protein sequences based on GFF coordinates.
2. Blast protein sequences against protein database.
3. Query protein sequences against interproscan databases.
4. Merge functional annotations.
