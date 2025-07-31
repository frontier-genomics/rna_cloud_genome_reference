nextflow.enable.dsl=2

include { DOWNLOAD_GENOME_AND_REFERENCES } from './subworkflows/download_genome_and_references.nf'
include { EXTRACT_PROTEIN_CODING_GENES } from './modules/grc_fixes.nf'
include { GRC_FIXES_ASSESSMENT } from './subworkflows/grc_fixes_assessment.nf'
include { SPLICE_SITE_GNOMAD_FREQ } from './subworkflows/splice_site_gnomad_freq.nf'
include { BUILD_GENOME_REFERENCE } from './subworkflows/genome_build.nf'
include { BUILD_ANNOTATION_REFERENCE } from './subworkflows/annotation_build.nf'

workflow {
    DOWNLOAD_GENOME_AND_REFERENCES()

    EXTRACT_PROTEIN_CODING_GENES(
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf
    )

    GRC_FIXES_ASSESSMENT(
        DOWNLOAD_GENOME_AND_REFERENCES.out.grc_fixes,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta_fai_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta_gzi_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report,
        EXTRACT_PROTEIN_CODING_GENES.out.protein_coding_genes,
        DOWNLOAD_GENOME_AND_REFERENCES.out.clinically_relevant
    )

    SPLICE_SITE_GNOMAD_FREQ(
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf_index,
        EXTRACT_PROTEIN_CODING_GENES.out.protein_coding_genes,
        DOWNLOAD_GENOME_AND_REFERENCES.out.clinically_relevant,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report
    )
    
    BUILD_GENOME_REFERENCE(
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta_fai_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta_gzi_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report,
        GRC_FIXES_ASSESSMENT.out.grc_fixes_assessment,
        DOWNLOAD_GENOME_AND_REFERENCES.out.cen_par_mask_regions,
        DOWNLOAD_GENOME_AND_REFERENCES.out.ebv_fasta
    )

    BUILD_ANNOTATION_REFERENCE(
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report,
        GRC_FIXES_ASSESSMENT.out.grc_fixes_assessment
    )
}