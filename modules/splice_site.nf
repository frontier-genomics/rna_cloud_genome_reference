nextflow.enable.dsl=2

process GET_CLINICALLY_SIGNIFICANT_PROTEIN_CODING_GENES {
    tag "GET_CLINICALLY_SIGNIFICANT_PROTEIN_CODING_GENES"
    label "python"
    publishDir "${params.temp_dir}", mode: 'copy'

    input:
    path protein_coding_genes_path
    path clinically_significant_genes_path
    path genome_regions_report_path

    output:
    path "clinically_significant_protein_coding_genes.tsv", emit: clinically_significant_protein_coding_genes

    script:
    """
    python -m rnacloud_genome_reference.splice_site_population_freq.get_clinically_significant_protein_coding_genes \
        --protein-coding-genes ${protein_coding_genes_path} \
        --clinically-significant-genes ${clinically_significant_genes_path} \
        --genome-regions-report ${genome_regions_report_path} \
        --output clinically_significant_protein_coding_genes.tsv
    """
}

process EXTRACT_SJ_POSITIONS_FROM_CLINICALLY_SIGNIFICANT_GENES {
    tag "EXTRACT_SJ_POSITIONS_FROM_CLINICALLY_SIGNIFICANT_GENES"
    label "python"
    publishDir "${params.temp_dir}", mode: 'copy'

    input:
    path clinical_genes_path
    path gtf
    path gtf_index

    output:
    path "sj_positions_from_clinically_significant_genes.tsv", emit: sj_positions

    script:
    """
    python -m rnacloud_genome_reference.splice_site_population_freq.extract_sj_pos \
        --clinical_genes_path ${clinical_genes_path} \
        --gtf_file_path ${gtf} \
        --output_path sj_positions_from_clinically_significant_genes.tsv
    """
}

process OET_SPLICE_SITE_GNOMAD_FREQ {
    tag "SPLICE_SITE_GNOMAD_FREQ"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path gnomad_freq
    val gnomad_freq_threshold
    val gnomad_hemizygote_count_threshold
    val gnomad_homozygote_count_threshold
    path splice_junctions

    output:
    path "splice_site_pop_freq.tsv", emit: splice_site_pop_freq

    script:
    """
    echo "🏃‍♂️ Creating DuckDB database from combined gnomAD frequency data and splice junctions"

    duckdb "gnomad.duckdb" <<EOF
    -- drop tables if they exist (so you can re-run without errors)
    DROP TABLE IF EXISTS gnomad_freq;
    DROP TABLE IF EXISTS splice_junctions;

    CREATE TABLE gnomad_freq AS
    SELECT *
    FROM read_csv(
        "${gnomad_freq}",
        delim='\t',
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
        "${splice_junctions}",
        delim='\t',
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

    echo "🏃‍♂️ Querying splice sites with high population frequency"

    duckdb "gnomad.duckdb" <<EOF
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
    TO "splice_site_pop_freq.tsv"
    (FORMAT CSV, DELIMITER '\t', HEADER true);
    EOF
    echo "✅ Splice sites with high population frequency written to splice_site_pop_freq.tsv"
    """
}