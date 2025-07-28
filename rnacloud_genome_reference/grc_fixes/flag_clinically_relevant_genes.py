
import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def flag_clinically_relevant_genes(gene_alt_contigs_comparison_file: str, clinically_relevant_genes_file: str, output_file: str) -> None:
    logger.info(f"Flagging clinically relevant genes in {gene_alt_contigs_comparison_file}")

    # Load gene-alt contigs comparison results
    comparison_results = pd.read_csv(gene_alt_contigs_comparison_file, sep="\t", low_memory=False)
    logger.info(f"Loaded {len(comparison_results)} gene-alt contigs comparison results.")

    # Load clinically relevant genes
    clinically_relevant_genes = pd.read_csv(clinically_relevant_genes_file, sep="\t", low_memory=False)
    logger.info(f"Loaded {len(clinically_relevant_genes)} clinically relevant genes.")

    comparison_results['clinically_relevant_gene'] = np.where(comparison_results['entrez_gene_id'].isin(clinically_relevant_genes['entrez_id']), True, False)
    logger.info(f"Writing output to {output_file}")
    comparison_results.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Flag clinically relevant genes in gene-alt contigs comparison results.')
    parser.add_argument('gene_alt_contigs_comparison_file', type=str, help='Path to the gene-alt contigs comparison file')
    parser.add_argument('clinically_relevant_genes_file', type=str, help='Path to the clinically relevant genes file')
    parser.add_argument('output_file', type=str, help='Path to the output file for flagged genes')

    args = parser.parse_args()

    flag_clinically_relevant_genes(args.gene_alt_contigs_comparison_file, args.clinically_relevant_genes_file, args.output_file)