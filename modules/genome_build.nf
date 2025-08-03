nextflow.enable.dsl=2

process CONVERT_GENOME_ANNOT_REFSEQ_TO_UCSC {
    tag "CONVERT_GENOME_ANNOT_REFSEQ_TO_UCSC"

    input:
    path fasta
    path fasta_fai_index
    path fasta_gzi_index
    path assembly_report

    output:
    path "${fasta.simpleName}_ucsc.fasta.gz", emit: fasta
    path "${fasta.SimpleName}_ucsc.fasta.gz.fai", emit: fasta_fai_index
    path "${fasta.SimpleName}_ucsc.fasta.gz.gzi", emit: fasta_gzi_index

    script:
    def base = fasta.simpleName

    """
    set -euo pipefail

    # Build a sed script to swap RefSeq → UCSC names
    awk -F '\\t' '
    BEGIN { print "/^>/!b" }
    !/^#/ {
        refseq = \$7
        ucsc   = \$NF
        sub(/\\r\$/, "", refseq)
        sub(/\\r\$/, "", ucsc)
        printf "s/^>%s.*/>%s/g\\n", refseq, ucsc
    }' ${assembly_report} > chrom_subst.sed

    # Apply name mapping, recompress, and index
    gunzip -c ${fasta} \
      | sed -f chrom_subst.sed \
      | bgzip -c > ${base}_ucsc.fasta.gz

    samtools faidx ${base}_ucsc.fasta.gz
    """
}

process GET_TARGET_CONTIGS {
    tag "GET_TARGET_CONTIGS"
    label "python"

    input:
    path assembly_report
    path grc_fixes_assessment

    output:
    stdout

    script:
    """
    set -euo pipefail

    python -m rnacloud_genome_reference.genome_build.get_target_contigs ${assembly_report} ${grc_fixes_assessment}
    """
}

process SUBSET_FASTA {
    tag "SUBSET_FASTA"

    input:
    path fasta
    val target_contigs

    output:
    path "${fasta.baseName}_subset.fasta", emit: fasta

    script:
    """
    set -euo pipefail

    samtools faidx ${fasta} ${target_contigs} > ${fasta.baseName}_subset.fasta
    """
}

process ADD_EBV {
    tag "ADD_EBV"

    input:
    path fasta
    path chrEBV_fasta

    output:
    path "${fasta.baseName}_with_ebv.fasta.gz", emit: fasta
    path "${fasta.baseName}_with_ebv.fasta.gz.fai", emit: fasta_fai_index
    path "${fasta.baseName}_with_ebv.fasta.gz.gzi", emit: fasta_gzi_index

    script:
    """
    set -euo pipefail

    echo "Adding EBV sequence to FASTA..."
    cat ${fasta} ${chrEBV_fasta} | bgzip -c > ${fasta.baseName}_with_ebv.fasta.gz

    echo "Indexing the combined FASTA file..."
    samtools faidx ${fasta.baseName}_with_ebv.fasta.gz
    """
}

process GRC_FIX_AND_ASSEMBLY_MASK_REGIONS {
    tag "GRC_FIX_AND_ASSEMBLY_MASK_REGIONS"
    label "python"

    input:
    path grc_fixes_assessment
    path gtf
    path gtf_index
    path cen_par_mask_regions

    output:
    path "grc_fixes_and_assembly_mask_regions.bed", emit: bed

    script:
    """
    set -euo pipefail
    python3 -m rnacloud_genome_reference.genome_build.generate_mask_bed ${grc_fixes_assessment} ${gtf} ${cen_par_mask_regions} grc_fixes_and_assembly_mask_regions.bed
    """
}

process REDUNDANT_5S_MASK_REGIONS {
    tag "REDUNDANT_5S_REGIONS"

    input:
    path gtf

    output:
    path "redundant_5s_regions.bed", emit: bed

    script:
    """
    set -euo pipefail
    echo "Extracting redundant 5S regions from GTF file..."
    gunzip -c ${gtf} | awk -F"\\t" '\$1=="NC_000001.11" && \$9~/gene_id "RNA5S([2-9]|1[0-7])"/ && \$3=="gene"' > redundant_5s_regions.gtf
    
    echo "Converting redundant 5S regions to BED format..."
    cat redundant_5s_regions.gtf | awk -F"\\t" 'OFS="\\t" { print "chr1", \$4-1, \$5, ".", 0, "." }' > redundant_5s_regions.bed
    """
}

process MASK_FASTA {
    tag "MASK_FASTA"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path grc_fixes_and_assembly_mask_regions_bed
    path redundant_5s_regions_bed
    path fasta
    path fasta_fai_index
    path fasta_gzi_index
    val final_output_prefix  // Prefix for output files

    output:
    path "${final_output_prefix}_rna_cloud.fasta.gz", emit: fasta
    path "${final_output_prefix}_rna_cloud.fasta.gz.fai", emit: fasta_fai_index
    path "${final_output_prefix}_rna_cloud.fasta.gz.gzi", emit: fasta_gzi_index

    script:
    """
    set -euo pipefail

    echo "Combining bed files for masking..."
    cat ${grc_fixes_and_assembly_mask_regions_bed} ${redundant_5s_regions_bed} | sort -u > combined_mask_regions.bed

    echo "Masking FASTA file with GRC fixes and redundant 5S regions..."
    bedtools maskfasta -fi <(gunzip -c ${fasta}) -bed combined_mask_regions.bed -fo ${final_output_prefix}_rna_cloud.fasta

    echo "Compressing the masked FASTA file..."
    bgzip ${final_output_prefix}_rna_cloud.fasta

    echo "Indexing the masked FASTA file..."
    samtools faidx ${final_output_prefix}_rna_cloud.fasta.gz
    """
}