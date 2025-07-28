nextflow.enable.dsl=2

include { DOWNLOAD_AND_INDEX_GENOME } from './modules/download.nf'
include { DOWNLOAD_AND_INDEX_GTF } from './modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_ASSEMBLY_REPORT } from './modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_GRC_FIXES } from './modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_CLINICALLY_RELEVANT_GENES } from './modules/download.nf'
include { EXTRACT_PROTEIN_CODING_GENES } from './modules/grc_fixes.nf'
include { SIMPLIFY_AND_ANNOTATE_GRC_FIXES } from './modules/grc_fixes.nf'
include { COMBINE_GRC_FIXES_AND_PROTEIN_CODING_GENES } from './modules/grc_fixes.nf'
include { COMPARE_FEATURES } from './modules/grc_fixes.nf'
include { FLAG_CLINICALLY_RELEVANT_GENES } from './modules/grc_fixes.nf'
include { DOWNLOAD_GENOME_AND_REFERENCES } from './subworkflows/download_genome_and_references.nf'
include { GRC_FIXES_ASSESSMENT } from './subworkflows/grc_fixes_assessment.nf'


workflow {
    println "🎬 Starting GRC Fixes Assessment Workflow"

    println "⬇️ Downloading genome and references"
    DOWNLOAD_GENOME_AND_REFERENCES()

    println "🏃‍♂️ Extracting protein coding genes"
    EXTRACT_PROTEIN_CODING_GENES(
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf
    )

    println "🏃‍♂️ Assessing GRC fixes"
    GRC_FIXES_ASSESSMENT(
        DOWNLOAD_GENOME_AND_REFERENCES.out.grc_fixes,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report,
        EXTRACT_PROTEIN_CODING_GENES.out.protein_coding_genes,
        DOWNLOAD_GENOME_AND_REFERENCES.out.clinically_relevant
    )

    println "🏁 GRC Fixes Assessment Workflow completed successfully"
}