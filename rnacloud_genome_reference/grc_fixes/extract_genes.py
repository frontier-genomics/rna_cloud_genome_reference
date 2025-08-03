import logging
import re
import gzip

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def extract_genes(gtf_file_path: str, output_file_path: str) -> None:
    logger.info(f"Extracting genes from {gtf_file_path} to {output_file_path}")

    with gzip.open(gtf_file_path, 'rt') as f:
        with open(output_file_path, 'w') as out_file:
            # Write output file header
            out_file.write("chr\tstart\tend\tstrand\tgene_name\tentrez_gene_id\tgene_biotype\n")

            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                
                if fields[2] == 'gene':
                    gene_biotype = re.search(r'gene_biotype \"(.+?)\";', fields[8]) # Extract gene biotype
                    gene_name = re.search(r'gene \"(.+?)\";', fields[8]) # Extract gene name
                    entrez_gene_id = re.search(r'\"GeneID:(.+?)\";', fields[8])

                    if not gene_biotype or not gene_name or not entrez_gene_id:
                        logger.error(f"Missing required fields in line: {line.strip()}")
                        raise ValueError("Missing required fields in GTF line")

                    out_file.write("{chr}\t{start}\t{end}\t{strand}\t{gene_name}\t{entrez_gene_id}\t{gene_biotype}\n".format(chr=fields[0],
                                                                                                      start=fields[3],
                                                                                                      end=fields[4],
                                                                                                      strand=fields[6],
                                                                                                      gene_name=gene_name.group(1),
                                                                                                      entrez_gene_id=entrez_gene_id.group(1),
                                                                                                      gene_biotype=gene_biotype.group(1)))
    logger.info(f"Genes extracted to {output_file_path}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Extract genes from a GTF file.')
    parser.add_argument('gtf_file', type=str, help='Path to the input GTF file (gzipped)')
    parser.add_argument('output_file', type=str, help='Path to the output file for genes')

    args = parser.parse_args()

    extract_genes(args.gtf_file, args.output_file)