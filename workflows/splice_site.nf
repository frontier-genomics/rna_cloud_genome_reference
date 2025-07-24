// Splice site population frequency workflow

nextflow.enable.dsl = 2

process EXTRACT_SPLICE_JUNCTIONS {
    tag "extract_sj"
    publishDir "${params.folders.temp_dir}", mode: 'copy'
    
    input:
    path gtf_file
    path clinically_relevant_genes
    
    output:
    path "clinically_significant_protein_coding_genes_sj_positions.tsv", emit: sj_positions
    
    script:
    """
    # Create necessary directories and copy input files
    mkdir -p ${params.folders.data_dir}
    mkdir -p ${params.folders.temp_dir}
    
    cp ${gtf_file} ${params.folders.data_dir}/
    cp ${clinically_relevant_genes} ${params.folders.data_dir}/clinically_relevant_genes/
    
    # Extract splice junction positions
    python -m rnacloud_genome_reference.splice_site_population_freq.extract_sj_pos
    """
}

process DOWNLOAD_GNOMAD_FREQ {
    tag "download_gnomad"
    publishDir "${params.folders.data_dir}/gnomad", mode: 'copy'
    
    output:
    path "**/*.tsv.gz", emit: gnomad_freq
    
    script:
    """
    # Create necessary directories
    mkdir -p ${params.folders.data_dir}
    mkdir -p ${params.folders.temp_dir}
    
    # Download gnomAD frequency data
    python -m rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq \\
        --chunk-size ${params.gnomad.chunk_size}
    """
}

process COMBINE_GNOMAD_FREQ {
    tag "combine_gnomad"
    publishDir "${params.folders.data_dir}/gnomad", mode: 'copy'
    
    input:
    path gnomad_files
    
    output:
    path "*.tsv.gz", emit: combined_freq
    
    script:
    """
    # Get variables from Python modules
    GNOMAD_DATA_PATH=\$(python -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_DATA_PATH; print(GNOMAD_DATA_PATH)")
    GNOMAD_VERSION=\$(python -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_VERSION; print(GNOMAD_VERSION)")
    GNOMAD_REFERENCE_GENOME=\$(python -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_REFERENCE_GENOME; print(GNOMAD_REFERENCE_GENOME)")
    GNOMAD_COMBINED_FILE="${GNOMAD_REFERENCE_GENOME}/${GNOMAD_VERSION}_freq.tsv.gz"
    
    # Run the combine script
    bash ${projectDir}/rnacloud_genome_reference/splice_site_population_freq/scripts/combine_gnomad_freq.sh \\
        \${GNOMAD_DATA_PATH} \\
        \${GNOMAD_COMBINED_FILE}
    """
}

process COMBINE_SPLICE_FREQ {
    tag "combine_splice_freq"
    publishDir "${params.folders.output_dir}", mode: 'copy'
    
    input:
    path sj_positions
    path gnomad_combined
    
    output:
    path "*_splice_site_pop_freq.tsv", emit: splice_freq
    
    script:
    """
    # Create necessary directories
    mkdir -p ${params.folders.data_dir}/gnomad
    mkdir -p ${params.folders.temp_dir}
    mkdir -p ${params.folders.output_dir}
    
    # Copy input files to expected locations
    cp ${sj_positions} ${params.folders.temp_dir}/
    cp ${gnomad_combined} ${params.folders.data_dir}/gnomad/
    
    # Get variables from Python modules
    GNOMAD_VERSION=\$(python -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_VERSION; print(GNOMAD_VERSION)")
    GNOMAD_REFERENCE_GENOME=\$(python -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_REFERENCE_GENOME; print(GNOMAD_REFERENCE_GENOME)")
    GNOMAD_COMBINED_FILE="${params.folders.data_dir}/gnomad/${GNOMAD_REFERENCE_GENOME}/${GNOMAD_VERSION}_freq.tsv.gz"
    GNOMAD_COMBINED_DUCKDB="${params.folders.data_dir}/gnomad/${GNOMAD_REFERENCE_GENOME}/${GNOMAD_VERSION}_freq.duckdb"
    OUTPUT="${GNOMAD_VERSION}_${GNOMAD_REFERENCE_GENOME}_splice_site_pop_freq.tsv"
    
    # Create DuckDB database and query
    duckdb "\$GNOMAD_COMBINED_DUCKDB" <<EOF
-- drop tables if they exist (so you can re-run without errors)
DROP TABLE IF EXISTS gnomad_freq;
DROP TABLE IF EXISTS splice_junctions;

CREATE TABLE gnomad_freq AS
SELECT *
FROM read_csv(
  "\${GNOMAD_COMBINED_FILE}",
  delim='\\t',
  header=TRUE,
  types={
    'chrom':'VARCHAR',
    'pos':'INTEGER',
    'ref':'VARCHAR',
    'alt':'VARCHAR',
    'lof_filter':'VARCHAR',
    'ac':'INTEGER',
    'an':'INTEGER',
    'hemizygote_count':'INTEGER',
    'homozygote_count':'INTEGER',
    'clinvar_variation_id':'INTEGER',
    'clinical_significance':'VARCHAR',
    'review_status':'VARCHAR'
  },
  compression='gzip'
);

CREATE TABLE splice_junctions AS
SELECT *
FROM read_csv(
  '${params.folders.temp_dir}/clinically_significant_protein_coding_genes_sj_positions.tsv',
  delim='\\t',
  header=TRUE,
  types={
    'chrom':'VARCHAR',
    'chrom_refseq':'VARCHAR',
    'pos':'INTEGER',
    'entrez_gene_id':'INTEGER',
    'gene_name':'VARCHAR',
    'transcript':'VARCHAR',
    'transcript_is_mane_select':'BOOLEAN',
    'exon_no':'INTEGER',
    'dist_from_annot':'INTEGER',
    'category':'VARCHAR',
  }
);
EOF

    # Query splice sites with high population frequency
    duckdb "\$GNOMAD_COMBINED_DUCKDB" <<EOF
COPY (
    WITH splice_site_pop_freq AS
            (SELECT sj.chrom,
                    sj.chrom_refseq,
                    sj.pos,
                    sj.entrez_gene_id,
                    sj.gene_name,
                    sj.transcript,
                    sj.transcript_is_mane_select,
                    sj.exon_no,
                    sj.dist_from_annot,
                    sj.category,
                    gf.ref,
                    gf.alt,
                    gf.lof_filter,
                    gf.ac,
                    gf.an,
                    gf.ac / gf.an AS af,
                    gf.hemizygote_count,
                    gf.homozygote_count,
                    gf.clinvar_variation_id,
                    gf.clinical_significance,
                    gf.review_status
            FROM splice_junctions sj
                    JOIN gnomad_freq gf
                            ON sj.chrom = gf.chrom
                                AND sj.pos = gf.pos)
    SELECT *
    FROM splice_site_pop_freq
    WHERE af > 0.1 OR homozygote_count > 100 OR hemizygote_count > 100
)
TO "\$OUTPUT"
  (FORMAT CSV, DELIMITER '\\t', HEADER true);
EOF

    # Clean up DuckDB file
    rm -f "\$GNOMAD_COMBINED_DUCKDB"
    """
}

workflow SPLICE_SITE_POPULATION_FREQ {
    take:
    gtf_file
    clinically_relevant_genes
    
    main:
    // Extract splice junction positions
    EXTRACT_SPLICE_JUNCTIONS(gtf_file, clinically_relevant_genes)
    
    // Download and combine gnomAD frequency data
    DOWNLOAD_GNOMAD_FREQ()
    COMBINE_GNOMAD_FREQ(DOWNLOAD_GNOMAD_FREQ.out.gnomad_freq.collect())
    
    // Combine splice junctions with gnomAD frequency data
    COMBINE_SPLICE_FREQ(
        EXTRACT_SPLICE_JUNCTIONS.out.sj_positions,
        COMBINE_GNOMAD_FREQ.out.combined_freq
    )
    
    emit:
    splice_freq = COMBINE_SPLICE_FREQ.out.splice_freq
}