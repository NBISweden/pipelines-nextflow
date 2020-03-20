# NBIS Annotation service nextflow pipelines

### Table of Contents

* [Disclaimer](#disclaimer)
* [Setup](#setup)
* [Available pipelines](#available-pipelines)

## Disclaimer

If you use these pipelines in your work, please acknowledge NBIS within your communication according to this example: "Support by NBIS (National Bioinformatics Infrastructure Sweden) is gratefully acknowledged."

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
# Install both nextflow and nf-core tools using conda
conda create -n nextflow-env nf-core nextflow
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

If running on a grid infrastructure, `nextflow` must be able to communicate
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

Nextflow is available under the module system on Uppmax, but could be outdated.

Nextflow scripts can be run in the following way.

```bash
module load bioinfo-tools Nextflow
NXF_VER=20.01.0 nextflow run [ -c <config> ] <nextflow_script> [ --script_parameters ]
```

This downloads the specific version of Nextflow locally for you to use before running the script. This version is cached in your `$HOME/.nextflow/` folder.

## Available pipelines

### Pipelines

See their respective README for operation instructions.

* [AnnotationPreprocessing.nf](AnnotationPreprocessing):
A pipeline for preprocessing genome assemblies in preparation for genome annotation.

* [AugustusTraining.nf](./AugustusTraining):
A pipeline for creating a training and testing data set for Augustus.

* [FunctionalAnnotationPreparation.nf](./FunctionalAnnotationPreparation):
A pipeline for functional annotation preparation.

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
