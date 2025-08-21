nextflow.enable.dsl=2

include { DOWNLOAD_EBV_GTF } from '../modules/download.nf'
include { CONVERT_ANNOTATION_REFSEQ_TO_UCSC } from '../modules/annotation_build.nf'
include { GET_TARGET_CONTIGS } from '../modules/genome_build.nf'
include { SUBSET_GTF } from '../modules/annotation_build.nf'
include { APPEND_GTFS as APPEND_RIBOSOMAL_RNA_GTFS } from '../modules/annotation_build.nf'
include { APPEND_GTFS as APPEND_EBV_GTF } from '../modules/annotation_build.nf'
include { DECOMPRESS_GTF } from '../modules/annotation_build.nf'
include { REMOVE_SECTIONS } from '../modules/annotation_build.nf'
include { SORT_GTF } from '../modules/annotation_build.nf'

workflow BUILD_ANNOTATION_REFERENCE {
    take:
    gtf
    assembly_report
    grc_fixes_assessment

    main:
    DOWNLOAD_EBV_GTF(
        params.genome.ebv_annotation_url
    )

    REMOVE_SECTIONS(
        gtf,
        Channel.value([["NC_000021.9", "rRNA"],
                       ["NT_187388.1", "rRNA"],
                       ["NT_167214.1", "rRNA"]])
    )

    APPEND_RIBOSOMAL_RNA_GTFS(
        "rRNA",
        REMOVE_SECTIONS.out.gtf,
        Channel.value([params.rRNA.NC_000021,
                       params.rRNA.NT_187388,
                       params.rRNA.NT_167214])
    )

    CONVERT_ANNOTATION_REFSEQ_TO_UCSC(
        APPEND_RIBOSOMAL_RNA_GTFS.out.gtf,
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
        CONVERT_ANNOTATION_REFSEQ_TO_UCSC.out.gtf_index,
        target_contigs
    )

    APPEND_EBV_GTF(
        "EBV",
        SUBSET_GTF.out.gtf,
        DOWNLOAD_EBV_GTF.out.gtf.toList()
    )

    SORT_GTF(
        final_output_prefix, // Prefix for output files
        APPEND_EBV_GTF.out.gtf
    )

    emit:
    gtf                   = SORT_GTF.out.gtf
    gtf_index             = SORT_GTF.out.gtf_index
}