import logging
import sqlite3

import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def combine_grc_fixes_and_protein_coding_genes(grc_fixes_file: str, protein_coding_genes_file: str, output_file: str) -> None:
    # Read GRC fixes and protein coding genes files
    grc_fixes_df = pd.read_csv(grc_fixes_file, sep='\t')
    protein_coding_genes_df = pd.read_csv(protein_coding_genes_file, sep='\t')

    # Write both to sqlite database
    conn = sqlite3.connect(":memory:")
    grc_fixes_df.to_sql('grc_fixes', conn, if_exists='replace', index=False)
    protein_coding_genes_df.to_sql('protein_coding_genes', conn, if_exists='replace', index=False)

    sql = """
        SELECT gf.chr_refseq,
        gf.chr_ucsc,
        pcg.start,
        pcg.end,
        pcg.strand,
        pcg.gene_name,
        pcg.entrez_gene_id,
        gf.issue_id,
        gf.type,
        gf.summary,
        gf.description,
        gf.alt_chr_refseq,
        gf.alt_chr_ucsc,
        gf.alt_scaf_start,
        gf.alt_scaf_stop
    FROM grc_fixes gf
    JOIN protein_coding_genes pcg
    ON pcg.end >= gf.parent_start
    AND pcg.start <= gf.parent_stop
    AND pcg.chr = gf.chr_refseq;
    """

    combined_df = pd.read_sql_query(sql, conn)
    combined_df.to_csv(output_file, sep='\t', index=False, header=True)
    logger.info(f"Combined GRC fixes and protein coding genes saved to {output_file}")

    conn.close()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Combine GRC fixes and protein coding genes into a single file.")
    parser.add_argument("grc_fixes_file", help="Path to the GRC fixes TSV file.")
    parser.add_argument("protein_coding_genes_file", help="Path to the protein coding genes TSV file.")
    parser.add_argument("output_file", help="Path to save the combined output TSV file.")

    args = parser.parse_args()
    combine_grc_fixes_and_protein_coding_genes(args.grc_fixes_file, args.protein_coding_genes_file, args.output_file)