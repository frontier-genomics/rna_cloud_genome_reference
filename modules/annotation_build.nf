nextflow.enable.dsl=2

process DECOMPRESS_GTF {
    tag "DECOMPRESS_GTF"

    input:
    path gtf

    output:
    path "${gtf.simpleName}.gtf", emit: gtf

    script:
    """
    set -euo pipefail

    echo "Decompressing GTF file"
    gunzip -c ${gtf} > ${gtf.simpleName}.gtf
    """
}

process REMOVE_SECTIONS {
    tag "REMOVE_SECTIONS"

    input:
    path gtf
    val pairs  // list of [contig, biotype] pairs

    output:
    path "${gtf.simpleName}.filtered.gtf", emit: gtf

    script:
    """
    set -euo pipefail

    # Print the filtering criteria
    echo "Removing multiple contig-biotype pairs from GTF: ${pairs}"

    # Copy the input GTF to a temporary file for iterative filtering
    cp ${gtf} temp.gtf

    # For each [contig, biotype] pair, filter out lines matching both criteria
    ${pairs.collect { pair -> 
        // Construct an awk command that excludes lines where:
        // - the first column matches the contig
        // - and the attribute field contains the specified biotype
        // Logic: keep line if it's NOT the matching biotype OR NOT the matching contig
        "awk '!(/\\_biotype \\\"${pair[1]}\\\"/) || \$1 != \"${pair[0]}\"' temp.gtf > temp2.gtf && mv temp2.gtf temp.gtf"
    }.join('\n')}

    # Rename the filtered GTF file to the final output
    mv temp.gtf ${gtf.simpleName}.filtered.gtf
    """
}

process APPEND_GTFS {
    tag "APPEND_GTFS"

    input:
    path gtf
    val additional_gtfs  // List of additional GTF paths as strings

    output:
    path "${gtf.baseName}.appended.gtf", emit: gtf

    script:
    """
    set -euo pipefail

    echo "Appending ${additional_gtfs.size()} additional GTF files to ${gtf}"

    # Check that all additional GTF files exist
    for f in ${additional_gtfs.join(' ')}; do
        if [ ! -f "\$f" ]; then
            echo "ERROR: Additional GTF file '\$f' not found." >&2
            exit 1
        fi
    done

    # Concatenate main GTF and all additional GTFs
    cat ${gtf} ${additional_gtfs.join(' ')} > "${gtf.baseName}.appended.gtf"
    """
}

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
    sed -f chrom_subst.sed ${gtf} > ${gtf.simpleName}_ucsc.gtf
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
    bedtools sort -i ${gtf} | bgzip -c > ${gtf}.gz

    echo "Indexing GTF file"
    tabix ${gtf}.gz

    echo "Filtering GTF for target contigs"
    tabix ${gtf}.gz ${target_contigs} | bgzip -c > ${gtf.simpleName}_subset.gtf.gz
    tabix ${gtf.simpleName}_subset.gtf.gz
    """
}