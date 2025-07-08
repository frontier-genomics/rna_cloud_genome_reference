from dataclasses import dataclass
from typing import Literal
import logging
import os
from pathlib import Path

import pandas as pd
import numpy as np

from rnacloud_genome_reference.config import Config
from rnacloud_genome_reference.common.gtf import Feature, Exon, Intron, GTFHandler, SpliceJunctionPosition
from rnacloud_genome_reference.common.gnomad import GenomicRegionQuerier, query_genomic_regions, GenomicRegion

logger = logging.getLogger(__name__)

logging.basicConfig(
    level=logging.DEBUG)

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

# Output format:
# Chrom
# gene
# entrez_gene_id
# transcript
# exon_no
# Category (D/A)
# pos
# alt
# ac
# an
# homozygous_count
# hemizygous_count

@dataclass
class SpliceJunctionPositionPopulationFrequency(SpliceJunctionPosition):
    alt: str
    ac: int
    an: int
    homozygous_count: int = 0
    hemizygous_count: int = 0
                
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
    print("ABC")
    get_clinically_significant_protein_coding_genes(
        protein_coding_genes_path=PROTEIN_CODING_GENES,
        clinically_significant_genes_path=os.path.join(CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER, CLINICALLY_RELEVANT_GENES_DESTINATION_FILE),
        genome_regions_report_path=os.path.join(GENOME_REPORT_DESTINATION_FOLDER, GENOME_REPORT_DESTINATION_FILE),
        output_path=os.path.join(TEMP_DIR, 'clinically_significant_protein_coding_genes.tsv')
    )