# Annotation preprocessing pipeline

## Quickstart

```
nextflow run -profile nbis,conda AnnotationPreprocessing.nf --genome_assembly '/path/to/genome_assembly.fasta'
```

Parameters can also be stored in a config file:
```bash
# Put all available parameter settings in a file.
grep "^params." AnnotationPreprocessing.nf > params.config
# Edit config file parameter values.
vim params.config
# Run workflow with config file.
nextflow run -c params.config -profile nbis,conda AnnotationPreprocessing.nf
```

Use `-resume` to restart failed jobs.
```bash
nextflow run -resume -c params.config -profile nbis,conda AnnotationPreprocessing.nf
```

## Parameters

* `genome_assembly`: The path to the genome assembly in quotes.
* `outdir`: The name of the results folder.
* `min_length`: The minimum_length for fasta sequences in the assembly to be. 

## Stages

* Filter: Remove fasta sequences less than `min_length` bases.
* Get Assembly Metrics: Calculate and plot summary metrics on the assembly.

