
from dataclasses import dataclass
from typing import Literal
import logging
import os
from pathlib import Path

import pandas as pd
import numpy as np

from rnacloud_genome_reference.config import Config
from rnacloud_genome_reference.grc_fixes.common import download_file, sort_index_gtf_file
from rnacloud_genome_reference.splice_site_population_freq.helper import get_clinically_significant_protein_coding_genes
from rnacloud_genome_reference.grc_fixes.assess_grc_fixes import ANNOTATION_DESTINATION_FILE, ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_SORTED_FILE, ANNOTATION_URL, CLINICALLY_RELEVANT_GENES_DESTINATION_FILE, CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER, CLINICALLY_RELEVANT_GENES_URL, GENOME_REPORT_DESTINATION_FILE, GENOME_REPORT_DESTINATION_FOLDER, GENOME_REPORT_URL, PROTEIN_CODING_GENES, TEMP_DIR
from rnacloud_genome_reference.gtf import Feature, Exon, Intron, GTFHandler, SpliceJunctionPosition, extract_protein_coding_genes
logger = logging.getLogger(__name__)

def extract_sj_positions_from_clinically_significant_genes(clinical_genes_path: str, gtf_file_path: str, output_path: str) -> None:
    logger.info("Extracting splice junction positions from clinically significant genes...")
    
    clinical_genes = pd.read_csv(clinical_genes_path, sep='\t', low_memory=False)
    
    gtf_file = GTFHandler(gtf_file_path=gtf_file_path)

    with open(output_path, 'w') as f:
        f.write("chrom\tchrom_refseq\tpos\tentrez_gene_id\tgene_name\ttranscript\ttranscript_is_mane_select\texon_no\tdist_from_annot\tcategory\n")

        for _, row in clinical_genes.iterrows():
            sj_positions = gtf_file.obtain_sj_positions(
                chrom=row['chrom_refseq'],
                start=row['start'],
                end=row['end'],
                entrez_gene_id=row['entrez_gene_id']
            )

            for sj_pos in sj_positions:
                f.write("{chrom}\t{chrom_refseq}\t{pos}\t{entrez_gene_id}\t{gene_name}\t{transcript}\t{transcript_is_mane_select}\t{exon_no}\t{dist_from_annot}\t{category}\n".format(
                    chrom=row['chrom'],
                    chrom_refseq=row['chrom_refseq'],
                    pos=sj_pos.pos,
                    entrez_gene_id=row['entrez_gene_id'],
                    gene_name=row['gene_name'],
                    transcript=sj_pos.transcript,
                    transcript_is_mane_select=sj_pos.transcript_is_mane_select,
                    dist_from_annot=sj_pos.dist_from_exon,
                    exon_no=sj_pos.exon_no,
                    category=sj_pos.category
                ))

if __name__ == "__main__":
    logger.info("Starting process to extract clinically significant protein-coding genes and their splice junction positions")

    logger.info("Downloading necessary files...")
    download_file(GENOME_REPORT_URL, GENOME_REPORT_DESTINATION_FOLDER, GENOME_REPORT_DESTINATION_FILE)
    
    download_file(ANNOTATION_URL, ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_FILE)
    sort_index_gtf_file(os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_FILE),
                        os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_SORTED_FILE))
    
    download_file(CLINICALLY_RELEVANT_GENES_URL, CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER, CLINICALLY_RELEVANT_GENES_DESTINATION_FILE)
    
    extract_protein_coding_genes(os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_FILE), PROTEIN_CODING_GENES)
    
    logger.info("Get clinically significant protein-coding genes...")
    get_clinically_significant_protein_coding_genes(
        protein_coding_genes_path=PROTEIN_CODING_GENES,
        clinically_significant_genes_path=os.path.join(CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER, CLINICALLY_RELEVANT_GENES_DESTINATION_FILE),
        genome_regions_report_path=os.path.join(GENOME_REPORT_DESTINATION_FOLDER, GENOME_REPORT_DESTINATION_FILE),
        output_path=os.path.join(TEMP_DIR, 'clinically_significant_protein_coding_genes.tsv')
    )

    logger.info("Extracting splice junction positions from clinically significant genes...")
    extract_sj_positions_from_clinically_significant_genes(
        clinical_genes_path=os.path.join(TEMP_DIR, 'clinically_significant_protein_coding_genes.tsv'),
        gtf_file_path=os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_SORTED_FILE),
        output_path=os.path.join(TEMP_DIR, 'clinically_significant_protein_coding_genes_sj_positions.tsv')
    )

    logger.info("Finished processing clinically significant protein-coding genes and their splice junction positions.")