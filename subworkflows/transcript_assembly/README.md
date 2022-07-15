# Transcript assembly pipeline

## Quickstart

Run workflow using the singularity profile

`params.yml`:

```yml
subworkflow: 'transcript_assembly'
reads: '/path/to/reads*_{R1,R2}.fastq.gz'
genome: '/path/to/genome.fasta'
single_end: false
outdir: '/path/to/save/results'
```

Command line:

```bash
nextflow run main.nf \
    -profile singularity \
    -params-file params.yml
```

## Parameters

- `reads`: Path to reads.
- `genome`: Path to genome.
- `single_end`: True if reads are single end reads, false if paired end (default: false).
- `outdir`: Path to save results.
- `skip_trimming`: True if trimming should be skipped (default: false).

### Tool specific parameters

Process specific options are passed by overriding the `ext.args` variable using a process selector in a configuration file.

`nextflow.config`:

```nextflow
process {
    // Override FASTP command line options
    withName: 'FASTP' {
        ext.args   = '-Q -L'
    }
}
```

See [Transcript Assembly modules config](../../config/transcript_assembly_modules.config) for the default tool configuration.
## Workflow Stages

1. Read QC.
    - **FastQC**: Reads properties are summarised and checked for common issues relating to sequence content and quality.
2. Read trimming (optional).
    - **Fastp**: Reads are trimmed for adapter read-through and a QC summary is then provided.
3. Guided assembly.
    - **Hisat2 Build**: Builds an index database for the input genome.
    - **Hisat2**: Align trimmed reads to the genome.
    - **Stringtie**: Assemble transcripts from the aligned reads.
4. Summary report.
    - **MultiQC**: Provide a consolidated report of the statistics from each previous tool.
