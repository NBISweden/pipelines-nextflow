# Functional annotation preparation pipeline

## Quickstart

```
nextflow run -profile nbis,conda FunctionalAnnotationPreparation.nf \
  --genome '/path/to/genome_assembly.fasta' \
  --gff_annotation 'path/to/annotation.gff3'
```

Parameters can also be stored in a config file:
```bash
# Put all available parameter settings in a file.
grep "^params\." FunctionalAnnotationPreparation.nf > params.config
# Edit config file parameter values.
vim params.config
# Run workflow with config file.
nextflow run -c params.config -profile nbis,conda FunctionalAnnotationPreparation.nf
```

Use `-resume` to restart failed jobs.
```bash
nextflow run -resume -c params.config -profile nbis,conda FunctionalAnnotationPreparation.nf
```

## Parameters

* `genome`: The path to the genome assembly in quotes.
* `gff_annotation`: The path to the gff annotation in quotes.
* `outdir`: The name of the results folder

* `records_per_file`: The number of records per file. A parallelisation option 
to improve blast and interproscan speed.

* `blast_db`: The path to protein database files in quotes.

## Stages


