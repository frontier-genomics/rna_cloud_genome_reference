from dataclasses import dataclass

GRC_FIXES_QUERY = '''
    (comparison_status == 'Different - Exons or introns sequences differ') or \
    (comparison_status == 'Different - No. of exons and introns differ') or \
    (comparison_status == 'Different - Discordant/flipped exon numbering') or \
    (comparison_status == 'Different - Splice sites differ')
'''

ASSEMBLY_REPORT_QUERY = '''
    `Sequence-Role` == 'assembled-molecule' or \
    `RefSeq-Accn` == 'NT_187388.1' or \
    `RefSeq-Accn` == 'NT_167214.1' or \
    `RefSeq-Accn` == 'NT_187633.1'
'''

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
    
def subtract_ranges(contig: tuple[int, int], genes: list[tuple[int, int]]) -> list[tuple[int, int]]:
    start, end = contig

    for g_start, g_end in genes:
        if g_start < start or g_end > end:
            raise ValueError(f"Gene range ({g_start}, {g_end}) is out of contig bounds ({start}, {end})")

    genes.sort()
    merged_genes = []
    for gene in genes:
        if not merged_genes or merged_genes[-1][1] < gene[0] - 1:
            merged_genes.append(gene)
        else:
            merged_genes[-1] = (merged_genes[-1][0], max(merged_genes[-1][1], gene[1]))

    result = []
    current = start
    for gene_start, gene_end in merged_genes:
        if current < gene_start:
            result.append((current, gene_start - 1))
        current = max(current, gene_end + 1)

    if current <= end:
        result.append((current, end))

    return result

def write_bed_file(*region_lists: list[Region], output_file: str):
    with open(output_file, 'w') as bed_file:
        for region_list in region_lists:
            for region in region_list:
                adjusted_start = region.start - 1  # Convert to 0-based index for BED format
                bed_file.write(f"{region.chrom}\t{adjusted_start}\t{region.end}\t{region.name}\t{region.score}\t{region.strand}\n")