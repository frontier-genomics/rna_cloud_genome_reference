nextflow.enable.dsl=2

include { CONVERT_GENOME_ANNOT_REFSEQ_TO_UCSC } from '../modules/genome_build.nf'
include { GET_TARGET_CONTIGS } from '../modules/genome_build.nf'
include { SUBSET_FASTQ } from '../modules/genome_build.nf'
include { ADD_EBV } from '../modules/genome_build.nf'
include { REDUNDANT_5S_MASK_REGIONS } from '../modules/genome_build.nf'
include { GRC_FIX_AND_ASSEMBLY_MASK_REGIONS } from '../modules/genome_build.nf'
include { MASK_FASTA } from '../modules/genome_build.nf'

workflow BUILD_GENOME_REFERENCE {
    take:
    fasta
    fasta_fai_index
    fasta_gzi_index
    gtf
    gtf_index
    assembly_report
    grc_fixes_assessment
    cen_par_mask_regions
    ebv_fasta

    main:
    CONVERT_GENOME_ANNOT_REFSEQ_TO_UCSC(
        fasta,
        fasta_fai_index,
        fasta_gzi_index,
        assembly_report
    )

    GET_TARGET_CONTIGS(
        assembly_report,
        grc_fixes_assessment
    )

    def target_contigs = GET_TARGET_CONTIGS.out

    SUBSET_FASTQ(
        CONVERT_GENOME_ANNOT_REFSEQ_TO_UCSC.out.fasta,
        target_contigs
    )

    ADD_EBV(
        SUBSET_FASTQ.out.fasta,
        ebv_fasta
    )

    REDUNDANT_5S_MASK_REGIONS(gtf)

    GRC_FIX_AND_ASSEMBLY_MASK_REGIONS(
        grc_fixes_assessment,
        gtf,
        gtf_index,
        cen_par_mask_regions
    )

    MASK_FASTA(
        GRC_FIX_AND_ASSEMBLY_MASK_REGIONS.out.bed,
        REDUNDANT_5S_MASK_REGIONS.out.bed,
        ADD_EBV.out.fasta,
        ADD_EBV.out.fasta_fai_index,
        ADD_EBV.out.fasta_gzi_index
    )
}