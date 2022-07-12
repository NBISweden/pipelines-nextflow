include { FASTQC              } from "$projectDir/modules/nf-core/modules/fastqc/main"
include { HISAT2_ALIGN        } from "$projectDir/modules/nf-core/modules/hisat2/align/main"
include { HISAT2_BUILD        } from "$projectDir/modules/nf-core/modules/hisat2/build/main"
include { FASTP               } from "$projectDir/modules/nf-core/modules/fastp/main"
include { STRINGTIE_STRINGTIE } from "$projectDir/modules/nf-core/modules/stringtie/stringtie/main"
include { MULTIQC             } from "$projectDir/modules/nf-core/modules/multiqc/main"

workflow TRANSCRIPT_ASSEMBLY {

    main:
    log.info """
        Transcript assembly using Hisat2/Stringtie workflow
        ===================================================
    """
    Channel.fromFilePairs( params.reads, size: params.single_end ? 1 : 2, checkIfExists: true)
        // .ifEmpty { error "Cannot find reads matching ${params.reads}!\n" }
        .set { reads }
    Channel.fromPath( params.genome, checkIfExists: true)
        // .ifEmpty { error "Cannot find genome matching ${params.genome}!\n" }
        .set { genome }

    FASTQC ( reads )
    HISAT2_BUILD ( genome )
    if ( !params.skip_trimming ){
            FASTP( reads )
            FASTP.out.reads.set { trimmed_reads }
            FASTP.out.log.set { trimming_logs }
    } else {
        reads.set { trimmed_reads }
    }
    HISAT2_ALIGN ( trimmed_reads, hisat2_index.out.collect() )
    STRINGTIE_STRINGTIE ( hisat2.out[0] )
    FASTQC.out.mix( trimming_logs ).mix( hisat2.out[2] ).set {logs}
    MULTIQC( logs.collect(),params.multiqc_config )

}