nextflow.enable.dsl=2

include { SIMPLIFY_AND_ANNOTATE_GRC_FIXES } from '../modules/grc_fixes.nf'
include { COMBINE_GRC_FIXES_AND_GENES } from '../modules/grc_fixes.nf'
include { COMPARE_FEATURES } from '../modules/grc_fixes.nf'
include { FLAG_CLINICALLY_RELEVANT_GENES } from '../modules/grc_fixes.nf'


workflow GRC_FIXES_ASSESSMENT {
    take:
    grc_fixes
    gtf
    gtf_index
    fasta
    fasta_fai_index
    fasta_gzi_index
    assembly_report
    genes
    clinically_relevant

    main:
    SIMPLIFY_AND_ANNOTATE_GRC_FIXES(
        grc_fixes,
        assembly_report
    )

    COMBINE_GRC_FIXES_AND_GENES(
        SIMPLIFY_AND_ANNOTATE_GRC_FIXES.out.simplified_grc_fixes,
        genes
    )

    COMPARE_FEATURES(
        gtf,
        fasta,
        COMBINE_GRC_FIXES_AND_GENES.out.combined_grc_fixes,
        gtf_index,
        fasta_fai_index,
        fasta_gzi_index
    )

    FLAG_CLINICALLY_RELEVANT_GENES(
        COMPARE_FEATURES.out.comparison_results,
        clinically_relevant
    )

    emit:
    grc_fixes_assessment = FLAG_CLINICALLY_RELEVANT_GENES.out.grc_fixes_assessment
}