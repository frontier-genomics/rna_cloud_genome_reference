import argparse
import logging

import pandas as pd

from rnacloud_genome_reference.common.gtf import GTFHandler
from rnacloud_genome_reference.common.utils import AssemblyReportParser
from rnacloud_genome_reference.genome_build.common import ASSEMBLY_REPORT_QUERY, GRC_FIXES_QUERY

logger = logging.getLogger(__name__)

def get_assembly_report_contigs(assembly_report: str, query: str) -> list[str]:
    df = pd.read_csv(assembly_report, 
                      sep='\t',
                      comment='#',
                      low_memory=False,
                      header=None,
                      names=['Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn','Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name'])
    df = df.query(query)
    if df.empty:
        logger.warning("No contigs found in assembly report with the specified query.")
        return []
    
    contigs = df['UCSC-style-name'].unique().tolist()

    return sorted(set(contigs))

def get_grc_fixes_contigs(grc_fixes_assessment: str, query: str) -> list[str]:
    grc = pd.read_csv(grc_fixes_assessment, sep='\t', low_memory=False)
    grc_filtered = grc.query(query)

    if grc_filtered.empty:
        logger.warning("No clinically relevant genes with discrepancies found in GRC fixes assessment.")
        return []

    # Both fix and primary contigs need to be added. There are a few instances where the fix contig
    # is for a gene annotated on an alt contig e.g. chr19_KI270866v1_alt -> chr19_MU273386v1_fix (HG-2469, GPI)
    contigs_fix = grc_filtered['alt_chr_ucsc'].unique().tolist()
    logger.debug(f"Number of fix contigs: {len(contigs_fix)}")

    contigs_primary = grc_filtered['chr_ucsc'].unique().tolist()
    logger.debug(f"Number of primary contigs: {len(contigs_primary)}")

    contigs_fix.extend(contigs_primary)
    
    return sorted(set(contigs_fix))

def get_target_contigs(assembly_report: str,
                       grc_fixes_assessment: str,
                       assembly_report_query: str,
                       grc_fixes_query: str) -> None:
    contigs_set_1 = get_assembly_report_contigs(assembly_report, assembly_report_query)
    logger.debug(f"No. of contigs contigs_set_1: {len(contigs_set_1)}")

    contigs_set_2 = get_grc_fixes_contigs(grc_fixes_assessment, grc_fixes_query)
    logger.debug(f"No. of contigs contigs_set_2: {len(contigs_set_2)}")

    contigs_set_1.extend(contigs_set_2)
    logger.debug(f"No. of contigs after merging: {len(contigs_set_1)}")

    unique_contigs = sorted(set(contigs_set_1))
    logger.debug(f"No. of unique contigs after merging: {len(unique_contigs)}")

    print(' '.join(unique_contigs), end='')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate contigs to be extracts from FASTA and GTF.")
    parser.add_argument("assembly_report", help="Path to the GRC assembly report.")
    parser.add_argument("grc_fixes_assessment", help="Path to the GRC fixes assessment TSV file.")
    parser.add_argument("--assembly_report_query", nargs='*', default=ASSEMBLY_REPORT_QUERY, help="Additional contigs to include.")
    parser.add_argument("--grc_fixes_query", default=GRC_FIXES_QUERY, help="Query to filter GRC fixes assessment.")

    args = parser.parse_args()

    get_target_contigs(args.assembly_report, args.grc_fixes_assessment, args.assembly_report_query, args.grc_fixes_query)