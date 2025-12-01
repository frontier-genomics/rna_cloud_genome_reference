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
    path assembly_report
    path grc_fixes_assessment
    path gtf
    path gtf_index
    path cen_par_mask_regions

    output:
    path "grc_fixes_and_assembly_mask_regions.bed", emit: bed

    script:
    """
    set -euo pipefail
    python3 -m rnacloud_genome_reference.genome_build.generate_mask_bed ${assembly_report} ${grc_fixes_assessment} ${gtf} ${cen_par_mask_regions} grc_fixes_and_assembly_mask_regions.bed
    """
}

process GRC_FIX_UNMASK_REGIONS {
    tag "GRC_FIX_UNMASK_REGIONS"
    label "python"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path grc_fixes_assessment
    path gtf
    path gtf_index

    output:
    path "unmasked_regions.bed", emit: bed

    script:
    """
    set -euo pipefail
    python3 -m rnacloud_genome_reference.genome_build.generate_unmask_bed ${grc_fixes_assessment} ${gtf} unmasked_regions.bed
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
    gunzip -c ${gtf} | awk -F"\t" '
    {
        chrom = \$1;                        # Chromosome name
        feature = \$3;                      # Feature type (e.g., gene, exon)
        attributes = \$9;                   # Attributes column
    }
    chrom == "NC_000001.11" && attributes ~ /gene_id "RNA5S([2-9]|1[0-7])"/ && feature == "gene"
    ' > redundant_5s_regions.gtf
    
    echo "Converting redundant 5S regions to BED format..."
    awk -F"\t" '
    {
        # Assign meaningful names to columns
        chrom = "chr1";                 # NOTE: This is hardcoded for chromosome 1 since the GTF extracted is for chromosome 1 only
        start = \$4 - 1;                # BED is 0-based, GTF is 1-based
        end = \$5;
        score = 0;
        strand = ".";

        # Extract gene_id from attributes (column 9)
        attributes = \$9;
        name = ".";
        n = split(attributes, attribute_array, ";");
        for (i = 1; i <= n; i++) {
            if (attribute_array[i] ~ /gene_id/) {
                # Remove label and quotes to get only the value
                sub(/gene_id "/, "", attribute_array[i]);
                sub(/"/, "", attribute_array[i]);
                name = attribute_array[i];
            }
        }

        # Set tab as output separator
        OFS = "\t";

        # Print in BED format: chrom, start, end, name, score, strand
        print chrom, start, end, name"-Redundant5S", score, strand
    }' redundant_5s_regions.gtf > redundant_5s_regions.bed
    """
}

process MASK_FASTA {
    tag "MASK_FASTA"
    publishDir "${params.output_dir}", mode: 'copy', pattern: "masked_regions.bed"

    input:
    path grc_fixes_and_assembly_mask_regions_bed
    path redundant_5s_regions_bed
    path fasta
    path fasta_fai_index
    path fasta_gzi_index

    output:
    path "masked.fasta.gz", emit: fasta
    path "masked.fasta.gz.fai", emit: fasta_fai_index
    path "masked.fasta.gz.gzi", emit: fasta_gzi_index
    path "masked_regions.bed", emit: mask_regions_bed

    script:
    """
    set -euo pipefail

    echo "Combining bed files for masking..."
    cat ${grc_fixes_and_assembly_mask_regions_bed} ${redundant_5s_regions_bed} | sort -u > masked_regions.bed

    echo "Masking FASTA file with GRC fixes and redundant 5S regions..."
    bedtools maskfasta -fi <(gunzip -c ${fasta}) -bed masked_regions.bed -fo masked.fasta

    echo "Compressing the masked FASTA file..."
    bgzip masked.fasta

    echo "Indexing the masked FASTA file..."
    samtools faidx masked.fasta.gz
    """
}

process SORT_FASTA {
    tag "SORT_FASTA"
    publishDir "${params.output_dir}", mode: 'copy'

    // Set default CPUs (2), can be overridden in nextflow.config or command line
    cpus 2

    input:
    val final_output_prefix  // Prefix for output files
    path fasta
    path fasta_fai_index
    path fasta_gzi_index

    output:
    path "${final_output_prefix}_rna_cloud.fasta.gz", emit: compressed_fasta
    path "${final_output_prefix}_rna_cloud.fasta.gz.fai", emit: compressed_fasta_fai_index
    path "${final_output_prefix}_rna_cloud.fasta.gz.gzi", emit: compressed_fasta_gzi_index
    path "${final_output_prefix}_rna_cloud.fasta", emit: uncompressed_fasta
    path "${final_output_prefix}_rna_cloud.fasta.fai", emit: uncompressed_fasta_fai_index

    script:
    """
    set -euo pipefail

    echo "Sorting and compressing RNA cloud FASTA file..."
    seqkit sort -j ${task.cpus} -N -o ${final_output_prefix}_rna_cloud.fasta ${fasta}

    echo "Compressing RNA cloud FASTA file..."
    bgzip -@ ${task.cpus} -c ${final_output_prefix}_rna_cloud.fasta > ${final_output_prefix}_rna_cloud.fasta.gz

    echo "Indexing uncompressed RNA cloud FASTA file..."
    samtools faidx ${final_output_prefix}_rna_cloud.fasta

    echo "Indexing compressed RNA cloud FASTA file..."
    samtools faidx ${final_output_prefix}_rna_cloud.fasta.gz
    """
}