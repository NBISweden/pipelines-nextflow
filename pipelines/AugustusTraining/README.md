# Augustus training pipeline

## Quickstart

```
nextflow run -profile nbis,conda AugustusTraining.nf \
  --genome '/path/to/genome_assembly.fasta' \
  --maker_evidence_gff 'path/to/annotation.gff3'
```

Parameters can also be stored in a config file:
```bash
# Put all available parameter settings in a file.
grep "^params." AugustusTraining.nf > params.config
# Edit config file parameter values.
vim params.config
# Run workflow with config file.
nextflow run -c params.config -profile nbis,conda AugustusTraining.nf
```

Use `-resume` to restart failed jobs.
```bash
nextflow run -resume -c params.config -profile nbis,conda AugustusTraining.nf
```

## Parameters

* `genome`: The path to the genome assembly in quotes.
* `maker_evidence_gff`: The path to the gff annotation in quotes.
* `outdir`: The name of the results folder

* `gff_gene_model_filter_options`: Options to be passed to the filter by gene model script (default:'-c -r -d 500 -a 0.3').

* `codon_table`: The number of the codon table to use for translation.

* `test_size`: The size of the test data set
* `flank_region_size`: The fize of the flank region to include. 

## Stages


