import logging
import re
import gzip
from dataclasses import dataclass

import pysam

logger = logging.getLogger(__name__)

@dataclass
class Feature:
    chromosome: str
    start: int
    end: int
    strand: str
    sequence: str | None

    def __repr__(self) -> str:
        sequence_display = self.sequence[:10] + "..." if self.sequence and len(self.sequence) > 10 else self.sequence
        return f"{self.__class__.__name__} (chromosome={self.chromosome}, start={self.start}, end={self.end}, strand={self.strand}, sequence={sequence_display})"

@dataclass
class Exon(Feature):
    exon_no: int

@dataclass
class Intron(Feature):
    intron_no: int

def extract_protein_coding_genes(gtf_file_path: str, output_file_path: str) -> None:
    logger.info(f"Extracting protein coding genes from {gtf_file_path} to {output_file_path}")
    
    with gzip.open(gtf_file_path, 'rt') as f:
        with open(output_file_path, 'w') as out_file:
            # Write output file header
            out_file.write("chr\tstart\tend\tstrand\tgene_name\tentrez_gene_id\n")

            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                
                if fields[2] == 'gene' and re.match(r'.*gene_biotype \"protein_coding\".*', fields[8]):
                    out_file.write("{chr}\t{start}\t{end}\t{strand}\t{gene_name}\t{entrez_gene_id}\n".format(chr=fields[0],
                                                                                                      start=fields[3],
                                                                                                      end=fields[4],
                                                                                                      strand=fields[6],
                                                                                                      gene_name=re.search(r'gene \"(.+?)\";', fields[8]).group(1),
                                                                                                      entrez_gene_id=re.search(r'\"GeneID:(.+?)\";', fields[8]).group(1)))
    logger.info(f"Protein coding genes extracted to {output_file_path}")

class GTFHandler:
    def __init__(self, gtf_file_path: str):
        self.gtf_file_path = gtf_file_path
        self.tbx: pysam.TabixFile = None # type: ignore

        try:
            self.tbx = pysam.TabixFile(gtf_file_path)
        except Exception as e:
            logger.error(f"Error opening GTF file {gtf_file_path}: {e}")
            raise

    def get_transcript_for_gene(self, chromosome: str, start: int, end: int, entrez_gene_id: int, mane: bool = True) -> str | None:
        logger.info(f"Obtaining transcript for Entrez Gene ID: {entrez_gene_id} at location {chromosome}:{start}-{end} (MANE: {mane})")
        
        transcript_id = None
        transcripts = []

        for record in self.tbx.fetch(chromosome, start, end, parser=pysam.asGTF()):
            if mane:
                if record['feature'] == 'transcript' and re.match(f'.*tag \"MANE Select\";.*', record['attributes']) and re.match(f'.*\"GeneID:{entrez_gene_id}\";.*', record['attributes']):
                    match = re.search(r'.*\"GenBank:(.+?)\";.*', record['attributes'])
                    if match:
                        transcripts.append(match.group(1))

            else:
                if record['feature'] == 'transcript' and not re.match(f'.*tag \"MANE Select\";.*', record['attributes']) and re.match(f'.*\"GeneID:{entrez_gene_id}\";.*', record['attributes']):
                    match = re.search(r'.*\"GenBank:(.+?)\";.*', record['attributes'])
                    if match:
                        transcripts.append(match.group(1))
                
        if len(transcripts) == 1:
            transcript_id = transcripts[0]
        elif len(transcripts) > 1:
            logger.warning(f'Multiple transcripts found for Entrez Gene ID: {entrez_gene_id} at location {chromosome}:{start}-{end}. Returning first transcript.')
            transcript_id = transcripts[0]
        elif len(transcripts) == 0:
            logger.warning(f'No transcript found for Entrez Gene ID: {entrez_gene_id} at location {chromosome}:{start}-{end}')
            return None
        
        logger.debug(f"Obtained transcript: {transcript_id}")
        return transcript_id
                    
    def get_exons_by_transcript(self, chromosome: str, start: int, end: int, transcript_id: str) -> list[Exon]:
        logger.info(f"Obtaining exons for transcript: {transcript_id} at location {chromosome}:{start}-{end}")
        
        exons = []

        for record in self.tbx.fetch(chromosome, start, end, parser=pysam.asGTF()):
            if record['feature'] == 'exon' and re.match(f'.*\"GenBank:{transcript_id}\";.*', record['attributes']):
                match = re.search(r'exon_number \"(\d+)\"', record['attributes'])
                
                if match:
                    exon_no = int(match.group(1))
                else:
                    logger.error("No exon_no found for {record}. Please check")
                    raise

                exons.append(Exon(chromosome=record['contig'],
                                  # Adding 1 to start is necessary because pysam is "pythonic" i.e. treates everything as 0 based (https://pysam.readthedocs.io/en/latest/faq.html#pysam-coordinates-are-wrong)
                                  start=int(record['start'] + 1),
                                  end=int(record['end']),
                                  strand=record['strand'],
                                  exon_no=exon_no,
                                  sequence=None))
        
        logger.debug(f"Obtained Exons: {exons}")
        return exons
    
    def get_exons_by_gene(self, chromosome: str, start: int, end: int, entrez_gene_id: int) -> list[Exon]:
        logger.info(f"Obtaining exons for Entrez Gene ID: {entrez_gene_id} at location {chromosome}:{start}-{end}")
        
        exons = []

        for record in self.tbx.fetch(chromosome, start, end, parser=pysam.asGTF()):
            if record['feature'] == 'exon' and re.match(f'.*\"GeneID:{entrez_gene_id}\";.*', record['attributes']):
                match = re.search(r'exon_number \"(\d+)\"', record['attributes'])
                
                if match:
                    exon_no = int(match.group(1))
                else:
                    logger.error("No exon_no found for {record}. Please check")
                    raise

                exons.append(Exon(chromosome=record['contig'],
                                  # Adding 1 to start is necessary because pysam is "pythonic" i.e. treates everything as 0 based (https://pysam.readthedocs.io/en/latest/faq.html#pysam-coordinates-are-wrong)
                                  start=int(record['start'] + 1),
                                  end=int(record['end']),
                                  strand=record['strand'],
                                  exon_no=exon_no,
                                  sequence=None))
        
        logger.debug(f"Obtained Exons: {exons}")
        return exons
    
    @staticmethod
    def derive_introns_from_exons(exons: list[Exon]) -> list[Intron]:
        if len(exons) == 0:
            logger.warning("No exons provided to derive introns from.")
            return []
        
        sorted_exons = sorted(exons, key=lambda x: (x.chromosome, x.start, x.exon_no))

        chromosome = sorted_exons[0].chromosome
        strand = sorted_exons[0].strand

        introns = []

        if len(sorted_exons) < 2:
            return introns
        
        for i in range(len(sorted_exons) - 1):
            current_exon = sorted_exons[i]
            next_exon = sorted_exons[i + 1]

            intron_entry = Intron(
                chromosome=chromosome,
                start=current_exon.end + 1,
                end=next_exon.start - 1,
                strand=strand,
                intron_no=current_exon.exon_no - 1 if strand == '-' else current_exon.exon_no,
                sequence=None
            )
            introns.append(intron_entry)

        return introns
