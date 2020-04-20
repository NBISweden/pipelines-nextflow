# Transcript assembly pipeline

## Quickstart (NBIS Staff)

```bash
module load Singularity
nextflow run -profile nbis,singularity /path/to/TranscriptAssembly.nf \
  --reads '/path/to/reads*_{R1,R2}.fastq.gz' \
  --genome 'path/to/genome.fasta'
```

Or:
```bash
nextflow run -profile nbis,conda /path/to/TranscriptAssembly.nf \
  --reads '/path/to/reads*_{R1,R2}.fastq.gz' \
  --genome 'path/to/genome.fasta'
```

## Usage

### Parameters

- General:
    * `reads`: Path to reads.
    * `genome`: Path to genome.
    * `single_end`: True if reads are single end reads, false if paired end (default: false).
    * `outdir`: Path to results folder.
    * `skip_trimming`: True if trimming should be skipped (default: false).
    * `trimmer`: Read trimming tool to use ( available: 'fastp' - default, 'trimmomatic').
- Fastp:
    * `fastp_options`: Command line options for fastp (default: ' -Q -L' - disable quality trimming; disable length trimming).
- Trimmomatic:
    * `trimmomatic_adapter_path`: Path to trimmomatic adapter sequences.
    * `trimmomatic_clip_options`: Trimmomatic clipping options ( default: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36' ).
- Hisat2:
    * `hisat2_options`: Command line options for hisat2, e.g. strandedness (`--hisat2_options ' --fr'`). **Note:** Quote the options and precede `--` with a space, otherwise nextflow interprets it as a workflow parameter. See the [Hisat2 Manual](https://ccb.jhu.edu/software/hisat2/manual.shtml) for the full range of options. (Default: `''`).
- Stringtie:
    * `stringtie_options`: Command line options for Stringtie, e.g. strandedness (`--stringtie_options ' --fr'`). **Note:** Quote the options and precede `--` with a space, otherwise nextflow interprets it as a workflow parameter. See the [StringTie Manual](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) for the full range of options. (Default: `''`).
- MultiQC:
    * `multiqc_config`: Path to MultiQC config (default: "$baseDir/config/multiqc_conf.yml").


Parameters to the workflow can be provided either using `--parameter` notation or via a config file as follows:

`params.config`:
```
// Workflow parameters
params.reads = '/path/to/reads*_{R1,R2}.fastq.gz'
params.genome = '/path/to/genome.fasta'
params.single_end = false
params.outdir = '/path/to/results'
params.skip_trimming = false
params.trimmer = 'fastp'
params.fastp_options = ' -Q -L'
// params.trimmomatic_adapter_path = '/path/to/trimmomatic/adapters.fasta'
// params.trimmomatic_clip_options = 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
params.hisat2_options = ''
params.stringtie_options = ''
params.multiqc_config = "$baseDir/config/multiqc_conf.yml"

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
nextflow run -c params.config -profile nbis,singularity /path/to/TranscriptAssembly.nf
```

## Workflow Stages

1. Read QC.
    * **FastQC**: Reads properties are summarised and checked for common issues relating to sequence content and quality.
2. Read trimming (optional; select one).
    * **Fastp**: (*default*) Reads are trimmed for adapter read-through and a QC summary is then provided.
    * **Trimmomatic**: Reads are trimmed for adapter read-through and low quality regions.
3. Guided assembly.
    * **Hisat2 Build**: Builds an index database for the input genome.
    * **Hisat2**: Align trimmed reads to the genome.
    * **Stringtie**: Assemble transcripts from the aligned reads.
4. Summary report.
    * **MultiQC**: Provide a consolidated report of the statistics from each previous tool.
