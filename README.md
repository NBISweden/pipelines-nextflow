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

### Nextflow:

Installation using conda:

```
# Install both nextflow and nf-core tools using conda
conda create -n nextflow-env nf-core nextflow
conda activate nextflow-env
```

Or:

```
# Install nextflow without using conda:
curl -fsSL get.nextflow.io | bash
mv ./nextflow ~/bin
```

#### Nextflow on Uppmax

Nextflow is available under the module system on Uppmax, but is outdated.

Nextflow scripts can be run in the following way.
```
module load bioinfo-tools Nextflow
NXF_VER=20.01.0 nextflow run [ -c <config> ] <nextflow_script> [ --script_parameters ]
```

### General Usage.

A workflow is run in the following way:
```
nextflow run [-profile <profile_name1>[,<profile_name2>,...] ] workflow.nf [--workflow_parameters]
```

## Available pipelines

See their respective README for operation instructions.

### [AnnotationPreprocessing.nf](AnnotationPreprocessing)

A pipeline for preprocessing genome assemblies in preparation for genome annotation.

### [AugustusTraining.nf](./AugustusTraining)

A pipeline for creating a training and testing data set for Augustus.

### [FunctionalAnnotationPreparation.nf](./FunctionalAnnotationPreparation)

A pipeline for functional annotation preparation.

### [TranscriptAssembly.nf](./TranscriptAssembly)

A transcript assembly pipeline using hisat2 and stringtie.
