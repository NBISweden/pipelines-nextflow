# NBIS Annotation service nextflow pipelines

### Table of Contents

* [Disclaimer](#disclaimer)
* [Maintainers](#maintainers)
* [Installation and Usage](#installation-and-usage)
* [Available pipelines](#available-pipelines)

## Disclaimer

If you use these pipelines in your work, please acknowledge NBIS within your communication according to this example: "Support by NBIS (National Bioinformatics Infrastructure Sweden) is gratefully acknowledged."

## Maintainers

* Mahesh Binzer-Panchal (Nextflow, Reproducibility)
* Jacques Dainat (Annotation expert, Nextflow, Conda)
* Lucile Soler (Annotation expert)

## Installation and Usage

Requirements:

* Nextflow
* A container platform (recommended for reproducibility)
    * Singularity
    * Docker
* The conda package manager if a container platform is not available.
* If conda, singularity, or docker is unavailable, all tool dependencies
must be installed in your PATH.

### Nextflow:

Installation using conda:

```bash
# Install nextflow from bioconda
conda create -c conda-forge -c bioconda -n nextflow-env nextflow
conda activate nextflow-env
```

Or:

```bash
# Install nextflow without using conda:
curl -fsSL get.nextflow.io | bash
mv ./nextflow ~/bin
```

### General Usage.

A workflow is run in the following way:
```bash
nextflow run [-profile <profile_name1>[,<profile_name2>,...] ] workflow.nf [--workflow_parameters]
```

If running on a compute cluster infrastructure, `nextflow` must be able to communicate
with the scheduler at all times, otherwise tasks will be cancelled.
The best way to do this is to run `nextflow` using a `screen` or `tmux`
terminal.

E.g. Screen
```bash
# Open a named screen terminal session
screen -S my_nextflow_run
# load nextflow with conda
conda activate nextflow-env
# run nextflow
nextflow run -c <config> -profile <profile> <nextflow_script>
# "Detach" screen terminal
<ctrl + a> <ctrl + d>
# list screen sessions
screen -ls
# "Attach" screen session
screen -r my_nextflow_run
```

#### Nextflow on Uppmax

Nextflow is available under the module system on Uppmax.

Nextflow scripts can be generally be run in the following way.

```bash
module load bioinfo-tools Nextflow
# Uses Nextflow version 20.01.0 - check Nextflow website for current release version
NXF_VER=20.01.0 nextflow run [ -c <config> ] <nextflow_script> [ --script_parameters ]
```

This downloads the specific version of Nextflow locally for you to use before running the script. This version is cached in your `$HOME/.nextflow/` folder.

In general, one should write their own configuration for nextflow scripts for their specific system.
The annotation pipelines come with prebuilt profiles that are useful on Uppmax systems. More specifically,
the profile `uppmax` is suitable to use as configuration for utilising the SLURM queuing system and the
singularity application to load necessary software. Alternatively the profile `nbis` and either
the profile `singularity` (recommended - the application `singularity` is available by default) or
`conda` (very slow - the application `conda` is loaded using the module system) to load the necessary
software needed for the workflows can be used. Additional configuration should also be added to utilise the Uppmax
clusters efficiently. In your own configuration file we suggest adding the following additional settings:

`workflow.config`:
```
// Workflow parameters
params.outdir = '/proj/<snic_storage_project>/results'

// Nextflow configuration options
workDir = '/proj/<snic_storage_project>/nobackup/work'
resume = true
process {
    clusterOptions = '-A <snic_compute_project>'
    // You can also override existing process cpu or time settings here too.
}
```
used like so:
```
NXF_VER=<version> nextflow run -c workflow.config -profile uppmax <nextflow_script>
```

Note: The FunctionalAnnotation pipeline needs one tool installed in the `PATH` along with its databases.
See [FunctionalAnnotation](./FunctionalAnnotation) for details.

## Available pipelines

### Pipelines

See their respective README for operation instructions.

* [AbinitioTraining.nf](./AbinitioTraining):
A pipeline for creating a training and testing data set for Augustus and Snap.

* [AnnotationPreprocessing.nf](AnnotationPreprocessing):
A pipeline for preprocessing genome assemblies in preparation for genome annotation.

* [FunctionalAnnotation.nf](./FunctionalAnnotation):
A pipeline for functional annotation.

* [TranscriptAssembly.nf](./TranscriptAssembly):
A transcript assembly pipeline using hisat2 and stringtie.

### General workflow profiles

* `uppmax`: A slurm and singularity profile for Uppmax clusters.
	- Usage: `nextflow run -c <config> -profile uppmax <nextflow_script>`.
* `nbis`: A slurm profile for the NBIS annotation cluster.
	- Usage: `nextflow run -c <config> -profile nbis,singularity <nextflow_script>`.
* `bils`: An LSF profile for the NBIS annotation cluster.
	- Usage: `nextflow run -c <config> -profile bils,conda <nextflow_script>`.
* `conda`: A conda software profile for use with compute infrastructures without `singularity` or `docker`.
	- Usage: `nextflow run -c <config> -profile nbis,conda <nextflow_script>`.
* `singularity`: A singularity software profile for compute infrastructures with `singularity` installed.
	- Usage: `nextflow run -c <config> -profile nbis,singularity <nextflow_script>`.
* `docker`: A docker software profile for compute infrastructures with `docker` installed.
	- Usage: `nextflow run -c <config> -profile docker <nextflow_script>`.
* `test`: Test profiles for each pipeline.
	- Usage: `nextflow run -profile nbis,singularity,test <nextflow_script>`.
