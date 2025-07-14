#!/bin/bash
set -e # exit on first error

echo "🏃‍♂️ Obtaining splice site positions"
python -m rnacloud_genome_reference.splice_site_population_freq.population_frequency

echo "🏃‍♂️ Running gnomAD frequency download and combination script"

COMBINE_GNOMAD_FREQ_SCRIPT="rnacloud_genome_reference/splice_site_population_freq/scripts/combine_gnomad_freq.sh"
GNOMAD_DATA_PATH=$(python -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_DATA_PATH; print(GNOMAD_DATA_PATH)")
GNOMAD_VERSION=$(python -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_VERSION; print(GNOMAD_VERSION)")
GNOMAD_REFERENCE_GENOME=$(python -c "from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import GNOMAD_REFERENCE_GENOME; print(GNOMAD_REFERENCE_GENOME)")
GNOMAD_COMBINED_FILE="data/gnomad/${GNOMAD_REFERENCE_GENOME}/${GNOMAD_VERSION}_freq.tsv.gz"
GNOMAD_COMBINED_DUCKDB="data/gnomad/${GNOMAD_REFERENCE_GENOME}/${GNOMAD_VERSION}_freq.duckdb"
OUTPUT="output/${GNOMAD_VERSION}_${GNOMAD_REFERENCE_GENOME}_splice_site_pop_freq.tsv"

if [ -f ${GNOMAD_COMBINED_FILE} ]; then
    echo "ℹ️ Combined gnomAD frequency file already exists. Skipping download and combination step."
else
    echo " Downloading gnomAD frequency data"
    python -m rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq

    if [ -z "$GNOMAD_DATA_PATH" ]; then
        echo "⛔️ GNOMAD_DATA_PATH is not set. Please check the download_gnomad_freq script."
        exit 1
    fi

    echo "---------⬇️----------------------------"
    echo "Combining gnomAD frequency data"
    echo "-------------------------------------"

    echo "Combining gnomAD frequency data into a single file."
    ${COMBINE_GNOMAD_FREQ_SCRIPT} ${GNOMAD_DATA_PATH} ${GNOMAD_COMBINED_FILE}
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
    'alt':'VARCHAR',
    'ac':'INTEGER',
    'an':'INTEGER',
    'hemizygote_count':'INTEGER',
    'homozygote_count':'INTEGER'
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
            (SELECT sj.*,
                    gf.alt,
                    gf.ac,
                    gf.an,
                    gf.ac / gf.an AS af,
                    gf.hemizygote_count,
                    gf.homozygote_count
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

echo "✅ Output written to $OUTPUT"

# Remove duckdb file if it exists
if [ -f "$GNOMAD_COMBINED_DUCKDB" ]; then
    echo "ℹ️ Removing DuckDB file: $GNOMAD_COMBINED_DUCKDB"
    rm "$GNOMAD_COMBINED_DUCKDB"
else
    echo "DuckDB file does not exist: $GNOMAD_COMBINED_DUCKDB"
fi
echo "ℹ️ Script completed successfully."

echo "You can now use the output file for further analysis: $OUTPUT"