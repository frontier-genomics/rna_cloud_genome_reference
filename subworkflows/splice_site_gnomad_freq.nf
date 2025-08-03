nextflow.enable.dsl = 2

include { GET_CLINICALLY_SIGNIFICANT_PROTEIN_CODING_GENES } from '../modules/splice_site.nf'
include { EXTRACT_SJ_POSITIONS_FROM_CLINICALLY_SIGNIFICANT_GENES } from '../modules/splice_site.nf'
include { OET_SPLICE_SITE_GNOMAD_FREQ } from '../modules/splice_site.nf'

workflow SPLICE_SITE_GNOMAD_FREQ {
    take:
    gtf
    gtf_index
    genes_path
    clinically_significant_genes_path
    genome_regions_report_path

    main:
    GET_CLINICALLY_SIGNIFICANT_PROTEIN_CODING_GENES(
        genes_path,
        clinically_significant_genes_path,
        genome_regions_report_path
    )

    EXTRACT_SJ_POSITIONS_FROM_CLINICALLY_SIGNIFICANT_GENES(
        GET_CLINICALLY_SIGNIFICANT_PROTEIN_CODING_GENES.out.clinically_significant_protein_coding_genes,
        gtf,
        gtf_index
    )

    OET_SPLICE_SITE_GNOMAD_FREQ(
        params.gnomad.reference,
        params.gnomad.freq,
        params.gnomad.hemizygote_count,
        params.gnomad.homozygote_count,
        EXTRACT_SJ_POSITIONS_FROM_CLINICALLY_SIGNIFICANT_GENES.out.sj_positions
    )

    emit:
    splice_site_pop_freq = OET_SPLICE_SITE_GNOMAD_FREQ.out.splice_site_pop_freq
}