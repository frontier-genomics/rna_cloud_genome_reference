#!/bin/bash
set -e # exit on first error

GNOMAD_DATA_PATH=$(python3 -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_DATA_PATH; print(GNOMAD_DATA_PATH)")
GNOMAD_VERSION=$(python3 -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_VERSION; print(GNOMAD_VERSION)")
GNOMAD_REFERENCE_GENOME=$(python3 -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_REFERENCE_GENOME; print(GNOMAD_REFERENCE_GENOME)")
GNOMAD_COMBINED_FILE="data/gnomad/${GNOMAD_REFERENCE_GENOME}/${GNOMAD_VERSION}_freq.tsv.gz"
GNOMAD_COMBINED_DUCKDB="data/gnomad/${GNOMAD_REFERENCE_GENOME}/${GNOMAD_VERSION}_freq.duckdb"
OUTPUT="output/${GNOMAD_VERSION}_${GNOMAD_REFERENCE_GENOME}_splice_site_pop_freq.tsv"

echo "🏃‍♂️ Starting script to combine gnomAD frequency data and splice junctions"
echo "  GNOMAD_DATA_PATH: $GNOMAD_DATA_PATH"
echo "  GNOMAD_VERSION: $GNOMAD_VERSION"
echo "  GNOMAD_REFERENCE_GENOME: $GNOMAD_REFERENCE_GENOME"
echo "  GNOMAD_COMBINED_FILE: $GNOMAD_COMBINED_FILE"
echo "  GNOMAD_COMBINED_DUCKDB: $GNOMAD_COMBINED_DUCKDB"
echo "  OUTPUT: $OUTPUT"

echo "🏃‍♂️ Obtaining splice site positions"
python3 -m rnacloud_genome_reference.splice_site_population_freq.extract_sj_pos
echo "✅ Splice site positions extracted to temp/clinically_significant_protein_coding_genes_sj_positions.tsv"

if [ ! -f ${GNOMAD_COMBINED_FILE} ]; then
    echo "⛔️ Combined gnomAD frequency file does not exist: ${GNOMAD_COMBINED_FILE}"
    exit 1
fi

echo "🏃‍♂️ Creating DuckDB database from combined gnomAD frequency data and splice junctions"
duckdb "$GNOMAD_COMBINED_DUCKDB" <<EOF
-- drop tables if they exist (so you can re-run without errors)
DROP TABLE IF EXISTS gnomad_freq;
DROP TABLE IF EXISTS splice_junctions;

CREATE TABLE gnomad_freq AS
SELECT *
FROM read_csv(
  "${GNOMAD_COMBINED_FILE}",
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
  'temp/clinically_significant_protein_coding_genes_sj_positions.tsv',
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
echo "✅ Finished creating DuckDB database from combined gnomAD frequency data."

echo "🏃‍♂️ Querying splice sites with high population frequency"
duckdb "$GNOMAD_COMBINED_DUCKDB" <<EOF
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
TO "$OUTPUT"
  (FORMAT CSV, DELIMITER '\t', HEADER true);
EOF
echo "✅ Splice sites with high population frequency written to $OUTPUT"

# Remove duckdb file if it exists
if [ -f "$GNOMAD_COMBINED_DUCKDB" ]; then
    echo "ℹ️ Removing DuckDB file: $GNOMAD_COMBINED_DUCKDB"
    rm "$GNOMAD_COMBINED_DUCKDB"
else
    echo "DuckDB file does not exist: $GNOMAD_COMBINED_DUCKDB"
fi

echo "🏁 Finished combining gnomAD frequency data and splice junctions"