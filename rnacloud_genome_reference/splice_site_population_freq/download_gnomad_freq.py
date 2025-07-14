import logging
import os

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

def split_ranges(start: int, stop: int, max_range: int = 50000) -> list[tuple[int, int]]:
    """
    Split a genomic range [start, stop] into subranges no larger than max_range.

    Args:
        start: The starting position (inclusive).
        stop: The ending position (inclusive).
        max_range: Maximum width of each subrange.

    Returns:
        A list of (sub_start, sub_stop) tuples.
    """
    ranges: list[tuple[int, int]] = []
    current_start = start
    while current_start <= stop:
        current_stop = min(current_start + max_range - 1, stop)
        ranges.append((current_start, current_stop))
        current_start = current_stop + 1
    return ranges

def download_gnomad_frequency(clinically_significant_protein_coding_genes: str,
                              gnomad_data_path: str = GNOMAD_DATA_PATH) -> None:
    logger.info("Loading clinically significant protein-coding genes...")
    data = pd.read_csv(clinically_significant_protein_coding_genes,
                       sep='\t', low_memory=False)

    gnomad_provider = GnomadProvider(reference_genome=GNOMAD_REFERENCE_GENOME,
                                     gnomad_version=GNOMAD_VERSION)

    os.makedirs(gnomad_data_path, exist_ok=True)

    for _, row in data.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        stop = int(row['end'])
        gene_id = row['entrez_gene_id']
        output_filename = os.path.join(
            gnomad_data_path,
            f"gnomad_frequencies_{chrom}_{gene_id}.tsv.gz"
        )

        logger.info(
            f"Processing chromosome: {chrom}, gene: {row['gene_name']}, entrez_gene_id: {gene_id}"
        )

        if os.path.exists(output_filename):
            logger.info(f"File {output_filename} already exists. Skipping.")
            continue

        try:
            total_range = stop - start + 1
            if total_range > 10000:
                sub_ranges = split_ranges(start, stop, 10000)
                logger.info(
                    f"Requested range {chrom}:{start}-{stop} (size={total_range}) "
                    f"exceeds 10000. Splitting into {len(sub_ranges)} sub-queries."
                )
            else:
                sub_ranges = [(start, stop)]

            all_variants: list[GnomadFrequency] = []
            for idx, (sub_start, sub_stop) in enumerate(sub_ranges, start=1):
                logger.info(
                    f"Querying gnomAD ({idx}/{len(sub_ranges)}) for region "
                    f"{chrom}:{sub_start}-{sub_stop}"
                )
                variants = gnomad_provider.query_gnomad(chrom, sub_start, sub_stop)
                if variants:
                    all_variants.extend(variants)
                else:
                    logger.warning(f"No gnomAD data returned for sub-range {sub_start}-{sub_stop}")

            if not all_variants:
                logger.warning(f"No variants found in gnomAD for {chrom}:{start}-{stop}")
                continue

            logger.info(f"Total variants found for {chrom}:{start}-{stop}: {len(all_variants)}")
            df = pd.DataFrame(all_variants)

            filtered = df.query('lof_filter.isna() and filters_count == 0')
            filtered[['chrom', 'pos', 'alt', 'ac', 'an', 'hemizygote_count', 'homozygote_count']] \
                .to_csv(
                    output_filename,
                    sep='\t',
                    index=False,
                    compression='gzip'
                )

        except Exception as e:
            logger.error(f"Error querying gnomAD for {chrom}:{start}-{stop}: {e}")


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

