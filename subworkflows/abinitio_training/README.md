# Abinitio training pipeline

The abinitio training workflow takes an assembly (parameter:`genome`) and
an evidence file from Maker (parameter:`maker_evidence_gff`) to filter
gene models and create training and test data sets for abinitio evidence
driven prediction in Maker.

## Quick start

Run workflow using the singularity profile:

`params.yml`:

```yml
subworkflow: 'abinitio_training'
genome: '/path/to/genome/assembly.fasta'
maker_evidence_gff: '/path/to/evidence/annotation.gff'
species_label: 'species_name'
codon_table: 1
aed_value:
  - 0.2
  - 0.3
locus_distance:
  - 3000
  - 4000
outdir: '/path/to/save/results'
```

To add the result folder to the augustus folder (to be used in maker for instance), add to the yml file :

```yml
maker_species_publishdir : '/PATH/augustus/config/species/'
```

Command line:

```bash
nextflow run NBISweden/pipelines-nextflow \
    -profile singularity \
    -params-file params.yml
```

## Parameters

- General:
  - `maker_evidence_gff`: Path to the GFF annotation.
  - `genome`: Path to the genome assembly.
  - `outdir`: Path to the results folder.
  - `species_label`: A species label for the training data.
  - `maker_species_publishdir`: A shared directory where a copy of the augustus `species_label` profile is saved.
  - `codon_table`: The number of the codon table to use for translation (default: 1).
  - `aed_value`: A list of model selection values to explore (smaller values mean higher stringency).
  - `locus_distance`: A list of locus distances (average distance between genes) to explore.
  - `flank_region_size`: The size of the flank region to include (default: 1000).

### Tool specific parameters

In these workflows, the Nextflow process directive `ext.args` is used to inject command line tool parameters directly to the shell script.
These command line tool parameters can be changed by overriding the `ext.args` variable for the respective process in a configuration file.

`nextflow.config`:

```nextflow
process {
    withName: 'MODEL_SELECTION_BY_AED' {
        ext.args = '--value 0.3 -a _AED -t ">"'
    }
}
```

See [Abinitio training modules config](../../config/abinitio_training_modules.config) for the default tool configuration.

## Workflow Stages

1. Separate maker evidence by record type.
2. Select model by AED.
3. Keep the longest isoform.
4. Remove incomplete gene models.
5. Filter by locus distance.
6. Extract the protein sequence.
7. Blast sequences against themselves.
8. Filter sequences.
9. Create a training and test dataset.
10. Train augustus.
11. Train snap.
