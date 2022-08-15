# Annotation preprocessing pipeline

The annotation preprocessing workflow cleans contigs in an assembly (parameter: `genome`),
calculates assembly statistics pre and post cleaning, along with Busco scores post cleaning
(parameter: `busco_lineage`).

## Quick start

Run workflow using the singularity profile:

`params.yml`:

```yml
subworkflow: 'annotation_preprocessing'
genome: '/path/to/genome/assembly.fasta'
busco_lineage:
  - 'eukaryota_odb10'
  - 'bacteria_odb10'
outdir: '/path/to/save/results'
```

Command line:

```bash
nextflow run NBISweden/pipelines-nextflow \
    -profile singularity \
    -params-file params.yml
```

## Parameters

- General:
  - `genome`: The path to the genome assembly in quotes.
  - `outdir`: The name of the results folder.
- Busco:
  - `busco_lineage`: The busco lineages to compare against (default: '[ 'eukaryota_odb10', 'bacteria_odb10' ]').
  - `busco_lineages_path`: The folder where busco lineages have been downloaded for shared use (default: unset -
    the selected busco lineage is downloaded by each process).

### Tool specific parameters

Process specific options are passed by overriding the `ext.args` variable using a process selector in a configuration file.

`nextflow.config`:

```nextflow
process {
    withName: 'ASSEMBLY_PURIFY' {
        ext.args   = '--size 1000'
    }
}
```

See [Annotation preprocessing modules config](../../config/annotation_preprocessing_modules.config) for the default tool configuration.

## Workflow Stages

1. Filter: Remove fasta sequences less than `min_length` bases.
2. Summarise and plot assembly metrics.
3. Run BUSCO on filtered assembly.

## Known issues

1. The Busco conda package does not resolve dependencies when `channel_priority: strict` is used.
