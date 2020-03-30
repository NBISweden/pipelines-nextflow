# Annotation preprocessing pipeline

## Quickstart (NBIS)

```bash
module load Singularity
nextflow run -profile nbis,singularity AnnotationPreprocessing.nf \
	--genome '/path/to/genome_assembly.fasta' \
	--outdir 'results' \
	--min_length '1000'
```

Or:
```bash
nextflow run -profile nbis,conda AnnotationPreprocessing.nf \
	--genome '/path/to/genome_assembly.fasta' \
	--outdir 'results' \
	--min_length '1000'
```


## Usage

### Parameters

- General:
	* `genome`: The path to the genome assembly in quotes.
	* `outdir`: The name of the results folder.
	* `min_length`: The minimum_length for fasta sequences in the assembly to be (default: 1000).
- Busco:
	* `busco_lineage`: The busco lineages to compare against (default: '[ 'eukaryota_odb10', 'bacteria_odb10' ]').

Parameters to the workflow can be provided either using `--parameter` notation or via a config file as follows:

`params.config`:
```
// Workflow parameters
params.genome = '/path/to/genome'
params.outdir = '/path/to/results'
params.min_length = 1000
	// Use `busco --list-datasets` for full list of available lineage sets
params.busco_lineage = [ 'eukaryota_odb10', 'bacteria_odb10' ]

// Nextflow parameters
resume = true
workDir = '/path/to/temporary/workspace'
conda.cacheDir = "$HOME/.nextflow/conda"
singularity.cacheDir = "$HOME/.nextflow/singularity"
```

Run nextflow with config file:
```bash
# Open screen terminal
screen -S my_nextflow_analysis
# Load Nextflow
conda activate nextflow-env
# Run Nextflow analysis
nextflow run -c params.config -profile nbis,singularity AnnotationPreprocessing.nf
```

### Workflow Stages

1. Filter: Remove fasta sequences less than `min_length` bases.
2. Summarise and plot assembly metrics.
3. Run BUSCO on filtered assembly.
