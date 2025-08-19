from dataclasses import dataclass

GRC_FIXES_QUERY = '''
    (comparison_status == 'Different - Sequences differ' and clinically_relevant_gene == True) or \
    (comparison_status == 'Different - Exon numbering is discordant' and clinically_relevant_gene == True) or \
    (comparison_status == 'Not comparable - Partial transcript annotation in GTF file' and clinically_relevant_gene == True and fix_contig_transcript_partial == False) or \
    (comparison_status == 'Different - No. of exons or introns differ' and clinically_relevant_gene == True)
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
    
def write_bed_file(*region_lists: list[Region], output_file: str):
    with open(output_file, 'w') as bed_file:
        for region_list in region_lists:
            for region in region_list:
                adjusted_start = region.start - 1  # Convert to 0-based index for BED format
                bed_file.write(f"{region.chrom}\t{adjusted_start}\t{region.end}\t{region.name}\t{region.score}\t{region.strand}\n")