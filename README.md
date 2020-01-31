# NBIS Annotation service nextflow pipelines

### Table of Contents

* [Disclaimer](#disclaimer)
* [Setup](#setup)
* [Available pipelines](#available-pipelines)

## Disclaimer

In the case where this work is used for communications (presentation, publication, etc.), we invite you to acknowledge NBIS within your communication according to this example: "Support by NBIS (National Bioinformatics Infrastructure Sweden) is gratefully acknowledged."


## Setup

It's expected that Nextflow is available on your system, otherwise it can
be installed using:
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

It might be good practice to set the `NXF_WORK` directory, otherwise nextflow creates
a work folder from where you run the nextflow pipeline.
```
# Set work dir to no-backup
export NXF_WORK=$SNIC_NOBACKUP/work
```

The general idea to run a workflow is in the following way:
```
nextflow run [-profile <profile_name1>[,<profile_name2>,...] ] workflow.nf [--workflow_parameters]
```

## Available pipelines

For a list of available pipelines see the [pipelines](./pipelines) folder.
