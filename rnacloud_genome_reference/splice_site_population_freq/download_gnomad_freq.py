import logging
import os
from typing import Optional

import pandas as pd

from rnacloud_genome_reference.grc_fixes.assess_grc_fixes import ANNOTATION_DESTINATION_FILE, ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_SORTED_FILE, ANNOTATION_URL, CLINICALLY_RELEVANT_GENES_DESTINATION_FILE, CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER, CLINICALLY_RELEVANT_GENES_URL, GENOME_REPORT_DESTINATION_FILE, GENOME_REPORT_DESTINATION_FOLDER, GENOME_REPORT_URL, PROTEIN_CODING_GENES, TEMP_DIR
from rnacloud_genome_reference.grc_fixes.common import download_file, sort_index_gtf_file
from rnacloud_genome_reference.gtf import extract_protein_coding_genes
from rnacloud_genome_reference.splice_site_population_freq.helper import get_clinically_significant_protein_coding_genes
from rnacloud_genome_reference.config import Config
from rnacloud_genome_reference.common.gnomad import GnomadFrequency, GnomadProvider, GNOMAD_VERSION, GNOMAD_REFERENCE_GENOME

logger = logging.getLogger(__name__)

DATA_DIR = Config.get_str('folders', 'data_dir')
TEMP_DIR = Config.get_str('folders', 'temp_dir')
OUTPUT_DIR = Config.get_str('folders', 'output_dir')

GNOMAD_DATA_PATH = os.path.join(DATA_DIR, 'gnomad', GNOMAD_REFERENCE_GENOME, GNOMAD_VERSION)

def download_gnomad_frequency(clinically_significant_protein_coding_genes: str,
                              gnomad_data_path: str = GNOMAD_DATA_PATH) -> None:
    logger.info("Loading clinically significant protein-coding genes...")
    data = pd.read_csv(clinically_significant_protein_coding_genes, sep='\t', low_memory=False)

    gnomad_provider = GnomadProvider(reference_genome=GNOMAD_REFERENCE_GENOME, gnomad_version=GNOMAD_VERSION)
    os.makedirs(gnomad_data_path, exist_ok=True)

    for _, row in data.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])
        entrez_gene_id = row['entrez_gene_id']

        output_filename = os.path.join(
            gnomad_data_path,
            f"gnomad_frequencies_{chrom}_{entrez_gene_id}.tsv.gz"
        )

        logger.info(
            f"Processing chromosome: {chrom}, gene: {row['gene_name']}, entrez_gene_id: {entrez_gene_id}"
        )

        if os.path.exists(output_filename):
            logger.info(f"File {output_filename} already exists. Skipping.")
            continue

        df = gnomad_provider.fetch_gnomad_stats_for_region(chrom, start, end, entrez_gene_id)
        if df is None:
            continue

        filtered = df.query('filters_count == 0')
        logger.info(f"Filtered {len(df)} rows to {len(filtered)} rows with no filters.")

        filtered[['chrom', 'pos', 'ref', 'alt', 'lof_filter', 'ac', 'an', 'hemizygote_count', 'homozygote_count', 'clinvar_variation_id', 'clinical_significance', 'review_status']].to_csv(
            output_filename,
            sep='\t',
            index=False,
            compression='gzip'
        )

if __name__ == "__main__":
    logger.info("Starting the gnomAD frequency download process...")

    logger.info("Downloading GRC annotation file...")
    download_file(ANNOTATION_URL, ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_FILE)
    sort_index_gtf_file(os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_FILE),
                        os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_SORTED_FILE))
    
    logger.info("Downloading GRC genome report...")
    download_file(GENOME_REPORT_URL, GENOME_REPORT_DESTINATION_FOLDER, GENOME_REPORT_DESTINATION_FILE)

    logger.info("Downloading clinically relevant genes file...")
    download_file(CLINICALLY_RELEVANT_GENES_URL, CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER, CLINICALLY_RELEVANT_GENES_DESTINATION_FILE)

    logger.info("Starting to extract protein-coding genes from GRC...")
    extract_protein_coding_genes(os.path.join(ANNOTATION_DESTINATION_FOLDER, ANNOTATION_DESTINATION_FILE), PROTEIN_CODING_GENES)

    logger.info("Starting to get clinically significant protein-coding genes...")
    get_clinically_significant_protein_coding_genes(
        protein_coding_genes_path=PROTEIN_CODING_GENES,
        clinically_significant_genes_path=os.path.join(CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER, CLINICALLY_RELEVANT_GENES_DESTINATION_FILE),
        genome_regions_report_path=os.path.join(GENOME_REPORT_DESTINATION_FOLDER, GENOME_REPORT_DESTINATION_FILE),
        output_path=os.path.join(TEMP_DIR, 'clinically_significant_protein_coding_genes.tsv')
    )

    logger.info("Downloading gnomAD frequency data for clinically significant protein-coding genes...")
    download_gnomad_frequency(
        clinically_significant_protein_coding_genes=os.path.join(TEMP_DIR, 'clinically_significant_protein_coding_genes.tsv'),
        gnomad_data_path=GNOMAD_DATA_PATH
    )

    logger.info("gnomAD frequency download process completed.")

