# Abinitio training pipeline

## Quickstart

```bash
nextflow run -profile nbis,singularity AbinitioTraining.nf \
  --genome '/path/to/genome_assembly.fasta' \
  --maker_evidence_gff 'path/to/annotation.gff3'
```

## Usage

### Parameters

- General:
    * `maker_evidence_gff`: Path to the GFF annotation.
    * `genome`: Path to the genome assembly.
    * `outdir`: Path to the results folder.
- Extract protein sequence:
    * `codon_table`: The number of the codon table to use for translation (default: 1).
- Augustus:
    * `flank_size`: The size of the flank region to include (default: 500).
    * `test_size`: The size of the test data set (default: 100).
    * `augustus_training_species`: A label for the Augustus `species` profile to use for training e.g. `'my_species_label'`.
    * `maker_species_publishdir`: A shared directory where a copy of the augustus `my_species_label` profile is saved.

Parameters to the workflow can be provided either using `--parameter` notation or via a config file as follows:

`params.config`:
```
// Workflow parameters
params.maker_evidence_gff = "/path/to/maker/evidence.gff"
params.genome = "/path/to/genome/assembly.fasta"
params.outdir = "results"
params.codon_table = 1
params.test_size = 100
params.flank_region_size = 500
params.augustus_training_species = ''  // e.g. 'asecodes_parviclava'
params.maker_species_publishdir = '/path/to/shared/maker/folder/' // e.g. '/projects/references/augustus/config/species/'

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
nextflow run -c params.config -profile nbis,singularity AbinitioTraining.nf
```

## Workflow Stages

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
