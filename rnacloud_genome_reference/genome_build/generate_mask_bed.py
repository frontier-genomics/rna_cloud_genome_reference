import argparse
from dataclasses import dataclass
import logging

import pandas as pd

from rnacloud_genome_reference.common.gtf import GTFHandler
from rnacloud_genome_reference.genome_build.get_target_contigs import GRC_FIXES_QUERY

logger = logging.getLogger(__name__)


@dataclass
class Region:
    chrom: str
    start: int # 1-based index
    end: int
    name: str = ''
    score: int | None = 0
    strand: str | None = '.'

    def __repr__(self) -> str:
        return (f"BedRegion(chrom={self.chrom}, start={self.start}, "
                f"end={self.end}, name={self.name}, score={self.score}, "
                f"strand={self.strand})")

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
            name=row['gene_name'],
            strand=row['strand']
        )
        mask_regions.append(region)

    gtf_handler = GTFHandler(gtf)

    for _, row in grc_filtered.iterrows():
        fix_contig_ucsc = row['alt_chr_ucsc']
        fix_contig_refseq = row['alt_chr_refseq']
        entrez_gene_id = row['entrez_gene_id']
        start = row['alt_scaf_start']
        end = row['alt_scaf_stop']
        gene_name = row['gene_name']

        fixed_gene = gtf_handler.get_gene_by_entrez_id(fix_contig_refseq, entrez_gene_id)

        if fixed_gene is None:
            logger.error(f"Could not find gene with Entrez ID {entrez_gene_id} in contig {fix_contig_refseq}")
            raise ValueError(f"Gene with Entrez ID {entrez_gene_id} not found in contig {fix_contig_refseq}")
        else:
            logger.info(f"Found gene {fixed_gene} in contig {fix_contig_refseq}")

        fix_contig_mask_regions = range_diff(
            start1=start,
            end1=end,
            start2=fixed_gene.start,
            end2=fixed_gene.end
        )

        if fix_contig_mask_regions is None:
            logger.warning(f"No range difference found for gene {gene_name} in contig {fix_contig_refseq}")
        else:
            for diff in fix_contig_mask_regions:
                mask_region = Region(
                    chrom=fix_contig_ucsc,
                    start=diff[0],
                    end=diff[1],
                    name=gene_name
                )
                mask_regions.append(mask_region)
                logger.info(f"Added mask region: {mask_region}")

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
            name=row['masked_scaf_accn']
        )
        regions.append(region)

    return regions    

def write_bed_file(*region_lists: list[Region], output_file: str):
    with open(output_file, 'w') as bed_file:
        for region_list in region_lists:
            for region in region_list:
                adjusted_start = region.start - 1  # Convert to 0-based index for BED format
                bed_file.write(f"{region.chrom}\t{adjusted_start}\t{region.end}\t{region.name}\t{region.score}\t{region.strand}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate GRC mask regions from GRC fixes assessment and GTF file.")
    parser.add_argument("grc_fixes_assessment", help="Path to the GRC fixes assessment TSV file.")
    parser.add_argument("gtf", help="Path to the GTF file.")
    parser.add_argument("cen_par_regions", help="Path to the centromere and PAR regions file.")
    parser.add_argument("output_file", default="mask_regions.bed", help="Output BED file name.")

    args = parser.parse_args()

    grc_fix_mask_regions = get_grc_mask_regions(args.grc_fixes_assessment, args.gtf, GRC_FIXES_QUERY)
    cen_par_mask_regions = get_cen_par_regions(args.cen_par_regions)

    write_bed_file(grc_fix_mask_regions, cen_par_mask_regions, output_file=args.output_file)
