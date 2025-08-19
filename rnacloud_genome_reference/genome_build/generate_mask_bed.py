import argparse
from dataclasses import dataclass
import logging

import pandas as pd

from rnacloud_genome_reference.common.gtf import GTFHandler
from rnacloud_genome_reference.genome_build.common import GRC_FIXES_QUERY, Region, subtract_ranges, write_bed_file

logger = logging.getLogger(__name__)

def range_diff(start1: int, end1: int, start2: int, end2: int) -> list[tuple[int, int]] | None:
    if start1 > end1:
        raise ValueError(f"Invalid first range: start1 ({start1}) > end1 ({end1})")
    if start2 > end2:
        raise ValueError(f"Invalid second range: start2 ({start2}) > end2 ({end2})")

    # No overlap (strictly before or after)
    if end1 < start2 or end2 < start1:
        logger.debug(f"No overlap between [{start1}, {end1}] and [{start2}, {end2}]")
        return None

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    diffs = []
    is_int = all(isinstance(v, int) for v in (start1, end1, start2, end2))

    if is_int:
        left_end = overlap_start - 1
        if left_end >= start1:
            diffs.append((start1, left_end))
        right_start = overlap_end + 1
        if right_start <= end1:
            diffs.append((right_start, end1))
    else:
        if start1 < overlap_start:
            diffs.append((start1, overlap_start))
        if overlap_end < end1:
            diffs.append((overlap_end, end1))

    return diffs

def get_grc_mask_regions(grc_fixes_assessment: str,
                         gtf: str,
                         query: str) -> list[Region]:
    mask_regions = []

    logger.info(f'Loading GRC fixes assessment from {grc_fixes_assessment}')
    grc = pd.read_csv(grc_fixes_assessment, sep='\t', low_memory=False)

    logger.info('Filtering GRC fixes assessment for clinically relevant genes with specific comparison statuses')
    grc_filtered = grc.query(query)
    logger.info(f'Found {len(grc)} clinically relevant genes with discrepancies')

    logger.info(f"Retrieving primary contig regions that are to be masked")
    for _, row in grc_filtered.iterrows():
        region = Region(
            chrom=row['chr_ucsc'],
            start=row['start'],
            end=row['end'],
            name=f"{row['gene_name']}-PRIMARY",
            strand=row['strand']
        )
        mask_regions.append(region)

    gtf_handler = GTFHandler(gtf)

    for _, contig in grc_filtered[['alt_chr_ucsc','alt_scaf_start','alt_scaf_stop']].drop_duplicates().iterrows():
        fix_contig_range = (contig['alt_scaf_start'], 
                            contig['alt_scaf_stop'])
        
        grc_fixes_for_contig = grc_filtered.query(f'alt_chr_ucsc == "{contig["alt_chr_ucsc"]}"')

        fix_genes_ranges = []

        for _, row in grc_fixes_for_contig.iterrows():
            fix_contig_ucsc = row['alt_chr_ucsc']
            fix_contig_refseq = row['alt_chr_refseq']
            entrez_gene_id = row['entrez_gene_id']

            fixed_gene = gtf_handler.get_gene_by_entrez_id(fix_contig_refseq, entrez_gene_id)

            if fixed_gene is not None:
                fix_genes_ranges.append((fixed_gene.start, fixed_gene.end))

        mask_ranges = subtract_ranges(fix_contig_range, fix_genes_ranges)

        gene_names = '-'.join(grc_fixes_for_contig['gene_name'].to_list())

        for gene_range in mask_ranges:
            mask_regions.append(Region(
                chrom=contig['alt_chr_ucsc'],
                start=gene_range[0],
                end=gene_range[1],
                name=f"{gene_names}-FIX"
            ))

    return mask_regions

def get_cen_par_regions(cen_par_regions: str) -> list[Region]:
    logger.info(f'Loading centromere and PAR regions from {cen_par_regions}')
    df = pd.read_csv(cen_par_regions, sep='\t', low_memory=False)
    
    regions = []
    for _, row in df.iterrows():
        region = Region(
            chrom=row['masked_copy_chr_name'],
            start=row['masked_copy_start'],
            end=row['masked_copy_stop'],
            name=f"{row['masked_scaf_accn']}-CEN_PAR"
        )
        regions.append(region)

    return regions    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate GRC mask regions from GRC fixes assessment and GTF file.")
    parser.add_argument("grc_fixes_assessment", help="Path to the GRC fixes assessment TSV file.")
    parser.add_argument("gtf", help="Path to the GTF file.")
    parser.add_argument("cen_par_regions", help="Path to the centromere and PAR regions file.")
    parser.add_argument("output_bed", default="mask_regions.bed", help="Output BED file containing regions that should be masked.")

    args = parser.parse_args()

    grc_fix_mask_regions = get_grc_mask_regions(args.grc_fixes_assessment, args.gtf, GRC_FIXES_QUERY)
    cen_par_mask_regions = get_cen_par_regions(args.cen_par_regions)

    write_bed_file(grc_fix_mask_regions, cen_par_mask_regions, output_file=args.output_bed)