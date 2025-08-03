nextflow.enable.dsl=2

include { CONVERT_ANNOTATION_REFSEQ_TO_UCSC } from '../modules/annotation_build.nf'
include { GET_TARGET_CONTIGS } from '../modules/genome_build.nf'
include { SUBSET_GTF } from '../modules/annotation_build.nf'
include { APPEND_GTFS } from '../modules/annotation_build.nf'
include { DECOMPRESS_GTF } from '../modules/annotation_build.nf'
include { REMOVE_SECTIONS } from '../modules/annotation_build.nf'

workflow BUILD_ANNOTATION_REFERENCE {
    take:
    gtf
    assembly_report
    grc_fixes_assessment

    main:
    DECOMPRESS_GTF(
        gtf
    )

    REMOVE_SECTIONS(
        DECOMPRESS_GTF.out.gtf,
        Channel.value([["NC_000021.9", "rRNA"],
                       ["NT_187388.1", "rRNA"],
                       ["NT_167214.1", "rRNA"]])
    )

    APPEND_GTFS(
        REMOVE_SECTIONS.out.gtf,
        Channel.value([params.rRNA.NC_000021,
                       params.rRNA.NT_187388,
                       params.rRNA.NT_167214])
    )

    CONVERT_ANNOTATION_REFSEQ_TO_UCSC(
        APPEND_GTFS.out.gtf,
        assembly_report
    )

    GET_TARGET_CONTIGS(
        assembly_report,
        grc_fixes_assessment
    )

    def target_contigs = GET_TARGET_CONTIGS.out

    // Obtain the final output prefix from the GTF filename
    def gtf_filename_from_url = "${params.genome.annotation_url.tokenize('/')[-1]}"
    def (full, final_output_prefix, suffix) = (gtf_filename_from_url =~ /(G.+p\d+)(_.+)/)[0]

    SUBSET_GTF(
        CONVERT_ANNOTATION_REFSEQ_TO_UCSC.out.gtf,
        target_contigs,
        final_output_prefix // Prefix for output files
    )

    emit:
    gtf                   = SUBSET_GTF.out.gtf
    gtf_index             = SUBSET_GTF.out.gtf_index
}