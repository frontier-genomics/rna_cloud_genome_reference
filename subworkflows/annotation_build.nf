nextflow.enable.dsl=2

include { CONVERT_ANNOTATION_REFSEQ_TO_UCSC } from '../modules/annotation_build.nf'
include { GET_TARGET_CONTIGS } from '../modules/genome_build.nf'
include { SUBSET_GTF } from '../modules/annotation_build.nf'
include { APPEND_GTF as APPEND_NC_000021 } from '../modules/annotation_build.nf'
include { APPEND_GTF as APPEND_NT_167214 } from '../modules/annotation_build.nf'
include { APPEND_GTF as APPEND_NT_187388 } from '../modules/annotation_build.nf'

workflow BUILD_ANNOTATION_REFERENCE {
    take:
    gtf
    assembly_report
    grc_fixes_assessment

    main:
    CONVERT_ANNOTATION_REFSEQ_TO_UCSC(
        gtf,
        assembly_report
    )

    APPEND_NC_000021(
        CONVERT_ANNOTATION_REFSEQ_TO_UCSC.out.gtf,
        Channel.fromPath(params.rRNA.NC_000021, checkIfExists: true)
    )

    APPEND_NT_167214(
        APPEND_NC_000021.out.gtf,
        Channel.fromPath(params.rRNA.NT_167214, checkIfExists: true)
    )

    APPEND_NT_187388(
        APPEND_NT_167214.out.gtf,
        Channel.fromPath(params.rRNA.NT_187388, checkIfExists: true)
    )

    GET_TARGET_CONTIGS(
        assembly_report,
        grc_fixes_assessment
    )

    def target_contigs = GET_TARGET_CONTIGS.out

    SUBSET_GTF(
        APPEND_NT_187388.out.gtf,
        target_contigs
    )
}