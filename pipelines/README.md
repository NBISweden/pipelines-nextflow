# NBIS Annotation service nextflow pipelines

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

These pipelines are translated from the bpipe pipelines in this repo along with some updates.

See their respective README for operation instructions.

### AnnotationPreprocessing.nf

A pipeline for preprocessing genome assemblies in preparation for genome annotation.

### AugustusTraining.nf

A pipeline for creating a training and testing data set for Augustus.

### FunctionalAnnotationPreparation.nf

A pipeline for functional annotation preparation.

### TranscriptAssemblyHisat2Stringtie.nf

A transcript assembly pipeline using hisat2 and stringtie.

