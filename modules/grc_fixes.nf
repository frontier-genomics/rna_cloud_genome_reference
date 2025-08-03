nextflow.enable.dsl=2

process EXTRACT_GENES {
    tag "EXTRACT_GENES"
    label "python"

    input:
    path gtf_file

    output:
    path "genes.tsv", emit: genes

    script:
    """
    set -euo pipefail
    python3 -m rnacloud_genome_reference.grc_fixes.extract_genes ${gtf_file} genes.tsv
    """
}

process SIMPLIFY_AND_ANNOTATE_GRC_FIXES {
    tag "SIMPLIFY_AND_ANNOTATE_GRC_FIXES"
    tag "python"

    input:
    path grc_fixes_file
    path assembly_report_file

    output:
    path "simplified_grc_fixes.tsv", emit: simplified_grc_fixes

    script:
    """
    set -euo pipefail
    python3 -m rnacloud_genome_reference.grc_fixes.simplify_and_annotate_grc_fixes ${grc_fixes_file} ${assembly_report_file} simplified_grc_fixes.tsv
    """
}

process COMBINE_GRC_FIXES_AND_GENES {
    tag "COMBINE_GRC_FIXES_AND_GENES"
    label "python"

    input:
    path genes_file
    path simplified_grc_fixes_file

    output:
    path "combined_grc_fixes.tsv", emit: combined_grc_fixes

    script:
    """
    set -euo pipefail
    python3 -m rnacloud_genome_reference.grc_fixes.combine_grc_fixes_and_genes ${genes_file} ${simplified_grc_fixes_file} combined_grc_fixes.tsv
    """
}

process COMPARE_FEATURES {
    tag "COMPARE_FEATURES"
    label "python"

    input:
    path gtf
    path fasta
    path combined_grc_fixes_file

    path gtf_index
    path fasta_fai_index
    path fasta_gzi_index

    output:
    path "comparison_results.tsv", emit: comparison_results

    script:
    """
    set -euo pipefail
    python3 -m rnacloud_genome_reference.grc_fixes.compare_features ${gtf} ${fasta} ${combined_grc_fixes_file} comparison_results.tsv
    """
}

process FLAG_CLINICALLY_RELEVANT_GENES {
    tag "FLAG_CLINICALLY_RELEVANT_GENES"
    label "python"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path comparison_results_file
    path clinically_relevant_genes_file

    output:
    path "grc_fixes_assessment.tsv", emit: grc_fixes_assessment

    script:
    """
    set -euo pipefail
    python3 -m rnacloud_genome_reference.grc_fixes.flag_clinically_relevant_genes ${comparison_results_file} ${clinically_relevant_genes_file} grc_fixes_assessment.tsv
    """
}