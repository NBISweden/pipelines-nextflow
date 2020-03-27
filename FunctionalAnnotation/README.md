# Functional annotation pipeline

## Quickstart

```bash
nextflow run -profile nbis,singularity FunctionalAnnotation.nf \
  --genome '/path/to/genome_assembly.fasta' \
  --gff_annotation 'path/to/annotation.gff3'
```

## Usage

### Parameters

- General:
    * `gff_annotation`:  Path to GFF genome annotation.
    * `genome`: Path to the genome assembly.
    * `outdir`: Path to the results folder.
    * `records_per_file`: Number of fasta records per file to distribute to blast and interproscan (default: 1000).
- GFF2Protein:
    * `codon_table`: (default: 1).
- Blastp:
    * `blast_db_fasta`: Path to protein database.
    * `blast_evalue`: Maximum expected value (E) for saving hits (default: '1e-6').
- Interproscan:
    * `interproscan_db`: Names of interproscan database to check against.
- Merge Functional Annotation:
    * `id_prefix`: Prefix set in the new identifider (default: 'NBIS').
    * `protein_existence`: Filter blast result by type of evidence that supports the existence of the protein in the uniprot DB (from 1 to 5, default: '5').

Parameters to the workflow can be provided either using `--parameter` notation or via a config file as follows:

`params.config`:
```
// Workflow parameters
params.gff_annotation = "/path/to/annotation.gff"
params.genome = "/path/to/genome.fasta"
params.outdir = "results"
params.codon_table = 1
params.records_per_file = 1000
params.blast_db_fasta = '/path/to/protein/database.fasta'
params.blast_evalue = '1e-6'
params.interproscan_db = ''
params.id_prefix = 'NBIS'
params.protein_existence = '5'

// Nextflow parameters
resume = true
workDir = '/path/to/temporary/workspace'
conda.cacheDir = '$HOME/.nextflow/conda'
singularity.cacheDir = '$HOME/.nextflow/singularity'
```

Run nextflow with config file:
```bash
# Open screen terminal
screen -S my_nextflow_analysis
# Load Nextflow
conda activate nextflow-env
# Run Nextflow analysis
nextflow run -c params.config -profile nbis,singularity FunctionalAnnotation.nf
```

## Workflow Stages

1. Extract protein sequences based on GFF coordinates.
2. Blast protein sequences against protein database.
3. Query protein sequences against interproscan databases.
4. Merge functional annotations.
