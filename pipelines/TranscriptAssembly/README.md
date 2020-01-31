# Transcript assembly pipeline

## Quickstart

```bash
nextflow run -profile nbis,conda TranscriptAssembly.nf \
  --reads '/path/to/reads*_{R1,R2}.fastq.gz' \
  --genome 'path/to/genome.fasta'
```

Parameters can also be stored in a config file:
```bash
# Put all available parameter settings in a file.
grep "^params." TranscriptAssembly.nf > params.config
# Edit config file parameter values.
vim params.config
# Run workflow with config file.
nextflow run -c params.config -profile nbis,conda TranscriptAssembly.nf
```

## Workflow description

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


## Parameters

| **General** | Description |
| :------- | :--- |
| `reads` | The path to the reads in quotes. The read pairs path must use the `{}` notation to define what a read pair is. |
| `genome` | The path to the genome assembly in quotes. |
| `single_end` | `true` if the reads are single-end, or `false` if reads are paired-end (default: `false`). |
| `outdir` | The name of the results folder. |
| `skip_trimming` | Skips trimming if true (default: `false`). |
| `trimmer` | The trimming tool to use (avail:  [`fastp` (default), `trimmomatic`]).
| **Fastp parameters** | |
| `fastp_options` | Default is ` -Q -L` which disables quality and length trimming by default. |
| **Trimmomatic parameters** | |
| `trimmomatic_adapter_path` | The path to the trimmomatic adapter file (must be specified when `trimmomatic` is selected as the trimmer).  |
| `trimmomatic_clip_options` | Read clipping options (Default:`'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'`). |
| **Hisat2 parameters** | |
| `hisat2_options` | Additional options for hisat2, e.g. strandedness (`--hisat2_options ' --fr'`). **Note:** Quote the options and preceed `--` with a space, otherwise nextflow interprets it as a workflow parameter. See the [Hisat2 Manual](https://ccb.jhu.edu/software/hisat2/manual.shtml) for the full range of options. (Default: `''`). |
| **StringTie parameters** | |
| `stringtie_options` | Additional options for stringtie, e.g. strandedness (`--stringtie_options ' --fr'`). **Note:** Quote the options and preceed `--` with a space, otherwise nextflow interprets it as a workflow parameter. See the [StringTie Manual](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) for the full range of options. (Default: `''`). |
