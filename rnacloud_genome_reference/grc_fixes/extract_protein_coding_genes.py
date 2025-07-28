import logging
import re
import gzip

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

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

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Extract protein coding genes from a GTF file.')
    parser.add_argument('gtf_file', type=str, help='Path to the input GTF file (gzipped)')
    parser.add_argument('output_file', type=str, help='Path to the output file for protein coding genes')

    args = parser.parse_args()

    extract_protein_coding_genes(args.gtf_file, args.output_file)