
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

def get_clinically_significant_protein_coding_genes(protein_coding_genes_path: str,
                                                    clinically_significant_genes_path: str,
                                                    genome_regions_report_path: str,
                                                    output_path: str) -> None:
    logger.info("Loading protein-coding genes and clinically significant genes...")
    protein_coding_genes = pd.read_csv(protein_coding_genes_path, sep="\t", low_memory=False)
    logger.info(f"Loaded {len(protein_coding_genes)} protein-coding genes.")
    
    logger.info("Loading clinically significant genes...")
    clinical_genes = pd.read_csv(clinically_significant_genes_path, sep="\t", low_memory=False)
    logger.info(f"Loaded {len(clinical_genes)} clinically significant genes.")

    logger.info("Flagging which protein-coding genes are clinically significant...")
    protein_coding_genes['clinically_relevant'] = np.where(
        protein_coding_genes['entrez_gene_id'].isin(clinical_genes.loc[clinical_genes.locus_group == "protein-coding gene", 'entrez_id']),
        True,
        False
    )

    logger.info("Loading genome regions report...")
    genome_regions = pd.read_csv(genome_regions_report_path, 
                      sep='\t',
                      comment='#',
                      low_memory=False,
                      header=None,
                      names=['Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn','Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name'])

    logger.info("Merging protein-coding genes with genome regions...")
    protein_coding_genes = protein_coding_genes.merge(
        genome_regions[['Sequence-Name','RefSeq-Accn','Sequence-Role']].rename(columns={'Sequence-Name': 'chrom', 
                                                                                'RefSeq-Accn': 'chr', 
                                                                                'Sequence-Role': 'role'}),
        how='left',
        on='chr'
    ).rename(columns={'chr': 'chrom_refseq'})

    logger.info("Filtering protein coding genes to those on primary contig and are clinically relevant...")
    pcg_of_interest = protein_coding_genes.query('role == "assembled-molecule" and chrom != "MT" and clinically_relevant == True')
    logger.info(f"Filtered down to {len(pcg_of_interest)} clinically relevant protein-coding genes on primary contig.")

    logger.info(f"Found {len(pcg_of_interest)} clinically significant protein-coding genes.")
    pcg_of_interest[['chrom', 'chrom_refseq','start', 'end', 'gene_name', 'entrez_gene_id']].to_csv(
        output_path,
        sep='\t',
        index=False
    )
    logger.info(f"Clinically significant protein-coding genes saved to {output_path}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Get clinically significant protein-coding genes.")
    parser.add_argument('--protein-coding-genes', type=str, required=True, help='Path to protein-coding genes file')
    parser.add_argument('--clinically-significant-genes', type=str, required=True, help='Path to clinically significant genes file')
    parser.add_argument('--genome-regions-report', type=str, required=True, help='Path to genome regions report file')
    parser.add_argument('--output', type=str, required=True, help='Output path for clinically significant protein-coding genes')

    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    get_clinically_significant_protein_coding_genes(args.protein_coding_genes, args.clinically_significant_genes, args.genome_regions_report, args.output)