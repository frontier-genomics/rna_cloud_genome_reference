import argparse
import logging

import pandas as pd

from rnacloud_genome_reference.common.gtf import GTFHandler
from rnacloud_genome_reference.common.utils import ChromosomeConverter

logger = logging.getLogger(__name__)

GRC_FIXES_QUERY = '''
    (comparison_status == 'Different - Sequences differ' and clinically_relevant_gene == True) or \
    (comparison_status == 'Different - Exon numbering is discordant' and clinically_relevant_gene == True) or \
    (comparison_status == 'Not comparable - Partial transcript annotation in GTF file' and clinically_relevant_gene == True and fix_contig_transcript_partial == False) or \
    (comparison_status == 'Different - No. of exons or introns differ' and clinically_relevant_gene == True)
'''

ADDITIONAL_CONTIGS = [
    "NT_187388.1",  # unlocalized rRNA
    "NT_167214.1",  # unplaced rRNA
    "NT_187633.1"   # GSTT1
]

def get_ucsc_contigs(assembly_report: str, refseq_contigs: list[str]) -> list[str]:
    df = pd.read_csv(assembly_report, 
                      sep='\t',
                      comment='#',
                      low_memory=False,
                      header=None,
                      names=['Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn','Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name'])
    df = df[df['RefSeq-Accn'].isin(refseq_contigs)]

    if df.empty:
        logger.warning(f"No contigs found in assembly report for {refseq_contigs}")
        return []

    return df['UCSC-style-name'].tolist()

def get_grc_fixes_contigs(grc_fixes_assessment: str, query: str) -> list[str]:
    """
    Retrieve contigs from the GRC fixes assessment file.
    """
    grc = pd.read_csv(grc_fixes_assessment, sep='\t', low_memory=False)
    grc_filtered = grc.query(query)

    if grc_filtered.empty:
        logger.warning("No clinically relevant genes with discrepancies found in GRC fixes assessment.")
        return []

    contigs = grc_filtered['alt_chr_ucsc'].unique().tolist()
    
    return sorted(set(contigs))

def get_target_contigs(assembly_report: str,
                       grc_fixes_assessment: str,
                       custom_contigs: list[str],
                       query: str):
    contigs_set_1 = get_ucsc_contigs(assembly_report, custom_contigs)
    logger.debug(f"No. of contigs contigs_set_1: {len(contigs_set_1)}")

    contigs_set_2 = get_grc_fixes_contigs(grc_fixes_assessment, query)
    logger.debug(f"No. of contigs contigs_set_2: {len(contigs_set_2)}")

    contigs_set_1.extend(contigs_set_2)
    logger.debug(f"No. of contigs after merging: {len(contigs_set_1)}")

    print(' '.join(contigs_set_1))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate contigs to be extracts from FASTA and GTF.")
    parser.add_argument("assembly_report", help="Path to the GRC assembly report.")
    parser.add_argument("grc_fixes_assessment", help="Path to the GRC fixes assessment TSV file.")
    parser.add_argument("--custom_contigs", nargs='*', default=ADDITIONAL_CONTIGS, help="Additional contigs to include.")
    parser.add_argument("--query", default=GRC_FIXES_QUERY, help="Query to filter GRC fixes assessment.")

    args = parser.parse_args()

    get_target_contigs(args.assembly_report, args.grc_fixes_assessment, args.custom_contigs, args.query)