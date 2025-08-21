import argparse
from dataclasses import dataclass
import logging

import pandas as pd

from rnacloud_genome_reference.common.gtf import GTFHandler
from rnacloud_genome_reference.genome_build.common import Region, write_bed_file, GRC_FIXES_QUERY

logger = logging.getLogger(__name__)

def get_fix_unmasked_regions(grc_fixes_assessment: str,
                         gtf: str,
                         query: str) -> list[Region]:
    unmasked_fix_regions = []

    logger.info(f'Loading GRC fixes assessment from {grc_fixes_assessment}')
    grc = pd.read_csv(grc_fixes_assessment, sep='\t', low_memory=False)

    logger.info('Filtering GRC fixes assessment for clinically relevant genes with specific comparison statuses')
    grc_filtered = grc.query(query)
    logger.info(f'Found {len(grc)} clinically relevant genes with discrepancies')

    logger.info(f"Retrieving fix contig regions that are to be kept")
    gtf_handler = GTFHandler(gtf)

    for _, row in grc_filtered.iterrows():
        fix_contig_ucsc = row['alt_chr_ucsc']
        fix_contig_refseq = row['alt_chr_refseq']
        entrez_gene_id = row['entrez_gene_id']
        gene_name = f"{row['gene_name']}-FIX"

        fixed_gene = gtf_handler.get_gene_by_entrez_id(fix_contig_refseq, entrez_gene_id)
        
        if fixed_gene is not None:
            unmasked_fix_regions.append(
                Region(
                    chrom=fix_contig_ucsc,
                    start=fixed_gene.start,
                    end=fixed_gene.end,
                    name=gene_name
                )
            )

    return unmasked_fix_regions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate GRC mask regions from GRC fixes assessment and GTF file.")
    parser.add_argument("grc_fixes_assessment", help="Path to the GRC fixes assessment TSV file.")
    parser.add_argument("gtf", help="Path to the GTF file.")
    parser.add_argument("output_bed", default="unmask_regions.bed", help="Output BED file containing regions that should be unmasked.")

    args = parser.parse_args()

    grc_fix_unmasked_regions = get_fix_unmasked_regions(args.grc_fixes_assessment, args.gtf, GRC_FIXES_QUERY)

    write_bed_file(grc_fix_unmasked_regions, output_file=args.output_bed)
