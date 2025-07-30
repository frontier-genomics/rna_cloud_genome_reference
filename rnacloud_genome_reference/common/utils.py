import logging
import pandas as pd

logger = logging.getLogger(__name__)

class ChromosomeConverter:
    def __init__(self, assembly_report: str):
        self.assembly_report = assembly_report
        self.refseq_to_ucsc_map = self._load_refseq_to_ucsc_map()
        self.ucsc_to_refseq_map = self._load_ucsc_to_refseq_map()

    def _load_refseq_to_ucsc_map(self) -> dict:
        # Load the RefSeq to UCSC map from the assembly report
        logger.info(f"Loading RefSeq to UCSC map from {self.assembly_report}")
        regions = pd.read_csv(self.assembly_report, 
                    sep='\t',
                    comment='#',
                    low_memory=False,
                    header=None,
                    names=['Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn','Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name'])
        
        chromosome_map = regions.set_index('RefSeq-Accn')['UCSC-style-name'].to_dict()

        logger.info(f"Loaded chromosome map: {chromosome_map}")
        return chromosome_map
    
    def _load_ucsc_to_refseq_map(self) -> dict:
        # Create a reverse map from UCSC to RefSeq
        logger.info("Creating UCSC to RefSeq map")
        ucsc_to_refseq_map = {v: k for k, v in self.refseq_to_ucsc_map.items()}
        logger.info(f"Created UCSC to RefSeq map: {ucsc_to_refseq_map}")
        return ucsc_to_refseq_map
    
    def refseq_to_ucsc(self, refseq_id: str) -> str:
        # Convert RefSeq ID to UCSC style name
        logger.debug(f"Converting RefSeq ID {refseq_id} to UCSC style name")
        ucsc_name = self.refseq_to_ucsc_map.get(refseq_id, None)
        if ucsc_name is None:
            logger.error(f"RefSeq ID {refseq_id} not found in RefSeq to UCSC map")
            raise ValueError(f"RefSeq ID {refseq_id} not found in RefSeq to UCSC map")
        return ucsc_name
    
    def ucsc_to_refseq(self, ucsc_name: str) -> str:        # Convert UCSC style name to RefSeq ID
        logger.debug(f"Converting UCSC style name {ucsc_name} to RefSeq ID")
        refseq_id = self.ucsc_to_refseq_map.get(ucsc_name, None)
        if refseq_id is None:
            logger.error(f"UCSC style name {ucsc_name} not found in UCSC to RefSeq map")
            raise ValueError(f"UCSC style name {ucsc_name} not found in UCSC to RefSeq map")
        return refseq_id
