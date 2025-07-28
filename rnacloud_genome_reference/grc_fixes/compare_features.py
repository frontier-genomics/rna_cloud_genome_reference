import logging

import pandas as pd

from rnacloud_genome_reference.grc_fixes.comparator import FeatureComparator

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def compare_features(gtf_file_path: str,
                     fasta_file_path: str,
                     gene_alt_contigs_mapping_file: str,
                     output_file: str) -> None:
    logger.info(f"Comparing features using GTF file: {gtf_file_path} and FASTA file: {fasta_file_path}")
    
    logger.info(f"Loading gene-alt contigs mapping from {gene_alt_contigs_mapping_file}")
    mappings = pd.read_csv(gene_alt_contigs_mapping_file, sep="\t", low_memory=False)

    comparator = FeatureComparator(gtf_file_path=gtf_file_path,
                                   fasta_file_path=fasta_file_path)

    logger.info(f"Comparing features for {len(mappings)} mappings.")
    mappings.join(
        mappings.apply(lambda x: pd.Series(comparator.compare_features(x['chr_refseq'],
                                                                       x['start'],
                                                                       x['end'],
                                                                       x['alt_chr_refseq'],
                                                                       x['alt_scaf_start'],
                                                                       x['alt_scaf_stop'],
                                                                       x['entrez_gene_id'])), axis=1)
    ).to_csv(output_file, sep="\t", index=False)

    logger.info(f"Feature comparison results saved to {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compare features in GTF and FASTA files.")
    parser.add_argument("gtf_file", help="Path to the GTF file.")
    parser.add_argument("fasta_file", help="Path to the FASTA file.")
    parser.add_argument("gene_alt_contigs_mapping_file", help="Path to the gene-alt contigs mapping file.")
    parser.add_argument("output_file", help="Path to save the comparison results.")

    args = parser.parse_args()

    compare_features(args.gtf_file, args.fasta_file, args.gene_alt_contigs_mapping_file, args.output_file)