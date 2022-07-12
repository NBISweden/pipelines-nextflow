#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ABINITIO_TRAINING        } from "$projectDir/subworkflows/abinitio_training"
include { ANNOTATION_PREPROCESSING } from "$projectDir/subworkflows/annotation_preprocessing"
include { FUNCTIONAL_ANNOTATION    } from "$projectDir/subworkflows/functional_annotation"
include { TRANSCRIPT_ASSEMBLY      } from "$projectDir/subworkflows/transcript_assembly"

workflow {

    log.info """
         _  _ ___ ___ ___ 
        | \| | _ )_ _/ __|
        | .` | _ \| |\__ \
        |_|\_|___/___|___/ Annotation Service

    """

    def valid_subworkflows = [ 'abinitio_training', 'annotation_preprocessing', 'functional_annotation', 'transcript_assembly' ]
    if( ! params.subworkflow in valid_subworkflows ){
        error """
        The parameter 'subworkflow' (value: ${params.subworkflow}) is not a valid subworkflow.
        
        Please select a valid subworkflow ([ ${valid_subworkflows.join(', ')} ]).
        """
    }

    if ( params.subworkflow == 'abinitio_training' ){
        ABINITIO_TRAINING()
    }

    if ( params.subworkflow == 'annotation_preprocessing' ){
        ANNOTATION_PREPROCESSING()
    }

    if ( params.subworkflow == 'functional_annotation' ){
        FUNCTIONAL_ANNOTATION()
    }
    
    if ( params.subworkflow == 'transcript_assembly' ){
        TRANSCRIPT_ASSEMBLY()
    }
}

workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Workflow completed successfully. 

        Thank you for using our workflow.
        Results are located in the folder: $params.outdir
        """)
    } else {
        log.info("""
        The workflow completed unsuccessfully.
        Please read over the error message. If you are unable to solve it, please
        post an issue at https://github.com/NBISweden/pipelines-nextflow/issues
        where we will do our best to help.
        """)
    }
}