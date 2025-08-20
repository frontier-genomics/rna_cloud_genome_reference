import logging

import pandas as pd

logger = logging.getLogger(__name__)

class GenomeReport:
    def __init__(self, assembly_report: str):
        self.assembly_report = assembly_report
        self.regions = self._load_assembly_report()

    def _load_assembly_report(self) -> pd.DataFrame:
        # Load the assembly report into a DataFrame
        logger.info(f"Loading assembly report from {self.assembly_report}")
        regions = pd.read_csv(self.assembly_report, 
                    sep='\t',
                    comment='#',
                    low_memory=False,
                    header=None,
                    names=['Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn','Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name'])
        return regions

    def get_contig_ranges(self, ucsc_contig_name: str) -> tuple[int, int]:
        # Get the start and end positions of a contig based on its UCSC style name
        logger.debug(f"Getting contig ranges for UCSC contig name {ucsc_contig_name}")
        contig_info = self.regions[self.regions['UCSC-style-name'] == ucsc_contig_name]
        
        if contig_info.empty:
            logger.error(f"UCSC contig name {ucsc_contig_name} not found in assembly report")
            raise ValueError(f"UCSC contig name {ucsc_contig_name} not found in assembly report")
        
        start = 1 # UCSC contigs are 1-based
        end = contig_info['Sequence-Length'].values[0]
        return start, end
