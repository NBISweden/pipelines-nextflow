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

## Workflow description

1. Separate maker evidence into .
2. Select model by AED.
3. Keep the longest isoform.
4. Remove incomplete gene models.
5. Filter by locus distance.
6. Extract the protein sequence.
7. Blast sequences against themselves.
8. Filter sequences.
9. Create a training and test dataset.
10. Train augustus using the species and training data.

## Parameters

| **General** | Description |
| :------- | :--- |
| `genome` | The path to the genome assembly in quotes. |
| `maker_evidence_gff` | The path to the gff annotation in quotes. |
| `outdir` | The name of the results folder |
| **Gene Model Filter parameters** | |
| `gff_gene_model_filter_options` | Options to be passed to the filter by gene model script (default:'-c -r -d 500 -a 0.3'). |
| **Protein Sequence Extraction parameters** | |
| `codon_table` | The number of the codon table to use for translation. |
| **Augustus parameters** | |
| `test_size` | The size of the test data set. |
| `flank_region_size` | The size of the flank region to include. |
| `augustus_training_species` | The name of the species (a folder with the name of the species and containing the augustus `species` profile will be created) e.g. `[ 'species1' ]` |
| `maker_species_publishdir` | The shared directory where a copy of the augustus `species` profile is saved. | 
