include { FASTQC              } from "$projectDir/modules/nf-core/fastqc/main"
include { HISAT2_ALIGN        } from "$projectDir/modules/local/hisat2/align"
include { HISAT2_BUILD        } from "$projectDir/modules/local/hisat2/build"
include { FASTP               } from "$projectDir/modules/nf-core/fastp/main"
include { STRINGTIE_STRINGTIE } from "$projectDir/modules/local/stringtie/stringtie"
include { MULTIQC             } from "$projectDir/modules/nf-core/multiqc/main"

workflow TRANSCRIPT_ASSEMBLY {

    main:
    log.info """
        Transcript assembly using Hisat2/Stringtie workflow
        ===================================================
    """
    Channel.fromFilePairs( params.reads, size: params.single_end ? 1 : 2, checkIfExists: true )
        .map { filestem, files -> [ [ id: filestem, single_end: params.single_end ], files ] }
        .set { reads }
    Channel.fromPath( params.genome, checkIfExists: true )
        .set { genome }

    FASTQC ( reads )
    HISAT2_BUILD ( genome )
    FASTP(  // Disabled using `when:` when params.skip_trimming
        reads, 
        [],     // Adapter file
        false,  // save trimmed fail
        false   // save merged
    ) // Disabled when params.skip_trimming
    HISAT2_ALIGN ( 
        params.skip_trimming ? reads : FASTP.out.reads,
        HISAT2_BUILD.out.index.collect()
    )
    STRINGTIE_STRINGTIE ( HISAT2_ALIGN.out.bam, [] )
    MULTIQC(
        FASTQC.out.zip.map{ meta, log -> log }.mix( 
            FASTP.out.log.map{ meta, log -> log }, 
            HISAT2_ALIGN.out.summary.map{ meta, log -> log } 
        ).collect(),
        file( params.multiqc_config, checkIfExists: true ),
        [], // extra MQC config
        []  // MQC logo
    )
}
