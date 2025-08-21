nextflow.enable.dsl=2

process VALIDATE_GENOME_ANNOTATION {
    tag "VALIDATE_GENOME_ANNOTATION"

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