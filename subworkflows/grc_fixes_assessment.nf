nextflow.enable.dsl=2

include { SIMPLIFY_AND_ANNOTATE_GRC_FIXES } from '../modules/grc_fixes.nf'
include { COMBINE_GRC_FIXES_AND_PROTEIN_CODING_GENES } from '../modules/grc_fixes.nf'
include { COMPARE_FEATURES } from '../modules/grc_fixes.nf'
include { FLAG_CLINICALLY_RELEVANT_GENES } from '../modules/grc_fixes.nf'


workflow GRC_FIXES_ASSESSMENT {
    take:
    grc_fixes
    gtf
    gtf_index
    fasta
    fasta_index
    assembly_report
    protein_coding_genes
    clinically_relevant

    main:
    SIMPLIFY_AND_ANNOTATE_GRC_FIXES(
        grc_fixes,
        assembly_report
    )

    COMBINE_GRC_FIXES_AND_PROTEIN_CODING_GENES(
        SIMPLIFY_AND_ANNOTATE_GRC_FIXES.out.simplified_grc_fixes,
        protein_coding_genes
    )

    COMPARE_FEATURES(
        gtf,
        fasta,
        COMBINE_GRC_FIXES_AND_PROTEIN_CODING_GENES.out.combined_grc_fixes,
        gtf_index,
        fasta_index
    )

    FLAG_CLINICALLY_RELEVANT_GENES(
        COMPARE_FEATURES.out.comparison_results,
        clinically_relevant
    )
}