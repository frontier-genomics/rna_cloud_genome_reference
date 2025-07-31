nextflow.enable.dsl=2

process CONVERT_ANNOTATION_REFSEQ_TO_UCSC {
    tag "CONVERT_GENOME_ANNOT_REFSEQ_TO_UCSC"

    input:
    path gtf
    path assembly_report

    output:
    path "${gtf.simpleName}_ucsc.gtf", emit: gtf

    script:
    """
    set -euo pipefail

    echo "Obtaining UCSC name substitutions from assembly report"
    awk -F '\t' '
    !/^#/ {
        refseq = \$7
        ucsc   = \$NF
        sub(/\\r\$/, "", refseq)
        sub(/\\r\$/, "", ucsc)
        printf "s/^%s/%s/\\n", refseq, ucsc
    }' ${assembly_report} > chrom_subst.sed

    echo "Substituting RefSeq names with UCSC names in GTF file"
    gunzip -c ${gtf} | sed -f chrom_subst.sed > ${gtf.simpleName}_ucsc.gtf
    """
}

process APPEND_GTF {
    tag "APPEND_GTF + ${additional_gtf.simpleName}"

    input:
    path gtf
    path additional_gtf

    output:
    path "${gtf.baseName}.${additional_gtf.simpleName}.gtf", emit: gtf

    script:
    """
    set -euo pipefail

    echo "Appending additional GTF to main GTF"
    cat ${gtf} ${additional_gtf} > "${gtf.baseName}.${additional_gtf.simpleName}.gtf"
    """
}

process SUBSET_GTF {
    tag "SUBSET_GTF"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path gtf
    val target_contigs

    output:
    path "${gtf.simpleName}_subset.gtf.gz", emit: gtf
    path "${gtf.simpleName}_subset.gtf.gz.tbi", emit: gtf_index

    script:
    """
    set -euo pipefail

    echo "Compressing GTF file"
    bgzip ${gtf}

    echo "Indexing GTF file"
    tabix ${gtf}.gz

    echo "Filtering GTF for target contigs"
    tabix ${gtf}.gz ${target_contigs} | bgzip -c > ${gtf.simpleName}_subset.gtf.gz
    tabix ${gtf.simpleName}_subset.gtf.gz
    """
}