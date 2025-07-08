import argparse
import logging
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from rnacloud_genome_reference.config import Config
from rnacloud_genome_reference.common.gtf import extract_protein_coding_genes
from rnacloud_genome_reference.grc_fixes.common import combine_grc_fixes_and_protein_coding_genes, download_file, index_fasta_file, sort_index_gtf_file
from rnacloud_genome_reference.grc_fixes.comparator import FeatureComparator

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)

DATA_DIR = Config.get_str('folders', 'data_dir')
TEMP_DIR = Config.get_str('folders', 'temp_dir')
OUTPUT_DIR = Config.get_str('folders', 'output_dir')

ANNOTATION_URL = Config.get_str('genome','annotation')
ANNOTATION_DESTINATION_FOLDER = os.path.join(DATA_DIR, os.path.basename(os.path.dirname(ANNOTATION_URL)))
ANNOTATION_DESTINATION_FILE = os.path.basename(ANNOTATION_URL)
ANNOTATION_DESTINATION_SORTED_FILE = (Path(ANNOTATION_DESTINATION_FILE).with_suffix('')).with_suffix('.sorted.gtf.gz').name

GENOME_URL = Config.get_str('genome', 'fasta')
GENOME_DESTINATION_FOLDER = os.path.join(DATA_DIR, os.path.basename(os.path.dirname(GENOME_URL)))
GENOME_DESTINATION_FILE = os.path.basename(GENOME_URL)

GENOME_REPORT_URL = Config.get_str('genome', 'assembly_report')
GENOME_REPORT_DESTINATION_FOLDER = os.path.join(DATA_DIR, os.path.basename(os.path.dirname(GENOME_REPORT_URL)))
GENOME_REPORT_DESTINATION_FILE = os.path.basename(GENOME_REPORT_URL)

GRC_FIXES_URL = Config.get_str('references', 'grc_fixes')
GRC_FIXES_DESTINATION_FOLDER = os.path.join(DATA_DIR, 'grc_fixes', os.path.basename(os.path.dirname(GRC_FIXES_URL)))
GRC_FIXES_DESTINATION_FILE = os.path.basename(GRC_FIXES_URL)

CLINICALLY_RELEVANT_GENES_URL = Config.get_str('references', 'clinically_relevant_genes')
CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER = os.path.join(DATA_DIR, 'clinically_relevant_genes', os.path.basename(os.path.dirname(CLINICALLY_RELEVANT_GENES_URL)))
CLINICALLY_RELEVANT_GENES_DESTINATION_FILE = os.path.basename(CLINICALLY_RELEVANT_GENES_URL)

PROTEIN_CODING_GENES = os.path.join(TEMP_DIR, 'protein_coding_genes.tsv')
SIMPLIFIED_GRC_FIXES = os.path.join(TEMP_DIR, 'simplified_grc_fixes.tsv')

GENE_ALT_CONTIGS_MAPPING = os.path.join(TEMP_DIR, 'gene_alt_contigs_mapping.tsv')
GENE_ALT_CONTIGS_COMPARISON = os.path.join(TEMP_DIR, 'gene_alt_contigs_mapping_comparison.tsv')

def check_and_create_folder(folder_path: str) -> None:
    """Check if a folder exists, if not create it."""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        logger.info(f"Created folder: {folder_path}")
    else:
        logger.info(f"Folder already exists: {folder_path}")

def simplify_and_annotate_grc_fixes(grc_fixes_path: str, assembly_report_path: str, output_path: str) -> None:
    grc_fixes_temp = pd.read_csv(grc_fixes_path, sep='\t', low_memory=False)
    logger.info(f"Loaded GRC fixes from {grc_fixes_path} with {len(grc_fixes_temp)} entries.")

    regions = pd.read_csv(assembly_report_path, 
                      sep='\t',
                      comment='#',
                      low_memory=False,
                      header=None,
                      names=['Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn','Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name'])
    logger.info("Loaded assembly report with {} entries.".format(len(regions)))

    grc_fixes = grc_fixes_temp.groupby(['parent_name', 'parent_start', 'parent_stop', 'ori', 'alt_scaf_acc', 'alt_scaf_start', 'alt_scaf_stop']).agg({
                    'issue_id': lambda x: ';'.join(x.astype(str)),
                    'type': lambda x: ';'.join(x.astype(str)),
                    'summary': lambda x: ';'.join(x.astype(str)),
                    'description': lambda x: ';'.join(x.astype(str))
                }).reset_index()
    logger.info(f"Grouped GRC fixes to {len(grc_fixes)} entries.")
    
    # Merge with assembly report to get UCSC-style names
    regions['parent_name'] = regions['Sequence-Name'].where(regions['Sequence-Name'] != regions['Assigned-Molecule'], regions['Assigned-Molecule'])
    
    # Merge GRC fixes with regions to get UCSC-style names
    grc_fixes = grc_fixes.merge(regions[['parent_name','UCSC-style-name','RefSeq-Accn']], on='parent_name', how='inner')

    # Rename columns for clarity
    grc_fixes.rename(columns={'UCSC-style-name': 'chr_ucsc',
                              'RefSeq-Accn': 'chr_refseq'}, inplace=True)
    
    # Rename alt_scaf_acc to GenBank-Accn for consistency
    grc_fixes.rename(columns={'alt_scaf_acc': 'GenBank-Accn'}, inplace=True)

    # Merge with regions again to get alternative chromosome names
    grc_fixes = grc_fixes.merge(regions[['GenBank-Accn','UCSC-style-name','RefSeq-Accn']], on='GenBank-Accn', how='inner')

    # Rename columns for alternative chromosome names
    grc_fixes.rename(columns={'UCSC-style-name': 'alt_chr_ucsc',
                              'RefSeq-Accn': 'alt_chr_refseq'}, inplace=True)
    
    # Check if any columns are NA
    if grc_fixes.isna().any().any():
        logger.error("NA values found in the grc_fixes DataFrame.")
        raise ValueError("There are NA values in the grc_fixes DataFrame. Please check the data.")

    grc_fixes.to_csv(output_path, sep='\t', index=False, header=True)
    logger.info(f"Simplified GRC fixes saved to {output_path}")

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
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="GRC fixes assessment pipeline")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file")
    
    args = parser.parse_args()

    # Validate output file path
    output_file = args.output
    if not output_file:
        logger.error("Output file path cannot be empty")
        sys.exit(1)
    
    logger.info("Starting GRC fixes assessment...")

    check_and_create_folder(DATA_DIR)
    check_and_create_folder(TEMP_DIR)
    check_and_create_folder(OUTPUT_DIR)

    # Download necessary files
    download_file(ANNOTATION_URL, ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_FILE)
    sort_index_gtf_file(os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_FILE),
                        os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_SORTED_FILE))
    download_file(GENOME_REPORT_URL, GENOME_REPORT_DESTINATION_FOLDER, GENOME_REPORT_DESTINATION_FILE)
    download_file(GENOME_URL, GENOME_DESTINATION_FOLDER, GENOME_DESTINATION_FILE)
    index_fasta_file(os.path.join(GENOME_DESTINATION_FOLDER, GENOME_DESTINATION_FILE))
    download_file(GRC_FIXES_URL, GRC_FIXES_DESTINATION_FOLDER, GRC_FIXES_DESTINATION_FILE)
    download_file(CLINICALLY_RELEVANT_GENES_URL, CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER, CLINICALLY_RELEVANT_GENES_DESTINATION_FILE)
    
    logger.info("Extracting protein coding genes from annotation file...")
    extract_protein_coding_genes(os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_FILE), PROTEIN_CODING_GENES)

    logger.info("Simplifying and annotating GRC fixes...")
    simplify_and_annotate_grc_fixes(grc_fixes_path=os.path.join(GRC_FIXES_DESTINATION_FOLDER, GRC_FIXES_DESTINATION_FILE),
                                    assembly_report_path=os.path.join(GENOME_REPORT_DESTINATION_FOLDER, GENOME_REPORT_DESTINATION_FILE),
                                    output_path=SIMPLIFIED_GRC_FIXES)

    logger.info("Combining GRC fixes with protein coding genes...")
    combine_grc_fixes_and_protein_coding_genes(grc_fixes_file=SIMPLIFIED_GRC_FIXES,
                                               protein_coding_genes_file=PROTEIN_CODING_GENES,
                                               output_file=GENE_ALT_CONTIGS_MAPPING)
    
    logger.info("Comparing primary and alternative contigs for genes...")
    compare_features(gtf_file_path=os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_SORTED_FILE),
                     fasta_file_path=os.path.join(GENOME_DESTINATION_FOLDER, GENOME_DESTINATION_FILE),
                     gene_alt_contigs_mapping_file=GENE_ALT_CONTIGS_MAPPING,
                     output_file=GENE_ALT_CONTIGS_COMPARISON)
    
    logger.info("Flag clinically relevant genes...")
    flag_clinically_relevant_genes(gene_alt_contigs_comparison_file=GENE_ALT_CONTIGS_COMPARISON,
                                   clinically_relevant_genes_file=os.path.join(CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER, CLINICALLY_RELEVANT_GENES_DESTINATION_FILE),
                                   output_file=output_file)
    
    logger.info(f"Pipeline completed successfully. Results written to: {output_file}")