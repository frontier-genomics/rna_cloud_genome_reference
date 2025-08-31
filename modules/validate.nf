nextflow.enable.dsl=2

process VALIDATE_GENOME_ANNOTATION {
    tag "VALIDATE_GENOME_ANNOTATION"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path fasta
    path fasta_index
    path gtf
    path gtf_index
    path masked_regions_bed
    path unmasked_regions_bed

    output:
    path "validation_report.txt", emit: report

    script:
    """
    set -euo pipefail

    /app/rnacloud_genome_reference/validation/scripts/validation.sh \
      ${fasta} \
      ${fasta_index} \
      ${gtf} \
      ${gtf_index} \
      ${masked_regions_bed} \
      ${unmasked_regions_bed} | tee validation_report.txt
    """
}

process GENOME_AND_ANNOTATION_REPORT {
    tag "GENOME_AND_ANNOTATION_REPORT"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    // Required genome + indices
    path fasta                                   // .fasta.gz
    path fasta_gzi                               // .fasta.gz.gzi
    path fasta_fai                               // .fasta.gz.fai

    // Annotation + index
    path gtf                                     // .gtf.gz
    path gtf_index                               // .gtf.gz.tbi

    // References & metadata
    path assembly_report
    path cen_par_regions
    path grc_fixes_summary
    path config_json                             // conf/sources.json

    // Genome and annotation version
    val version

    output:
    path 'genome_and_annotation_report.md', emit: report

    script:
    def genome_and_annotation_version = version ?: "0.0.0"
    """
    set -euo pipefail

    python -m rnacloud_genome_reference.validation.report \
      --fasta ${fasta} \
      --gtf ${gtf} \
      --assembly-report ${assembly_report} \
      --cen-par-regions ${cen_par_regions} \
      --grc-fixes-summary ${grc_fixes_summary} \
      --config ${config_json} \
      --version ${genome_and_annotation_version} \
      --out genome_and_annotation_report.md
    """
}