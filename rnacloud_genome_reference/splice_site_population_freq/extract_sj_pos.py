
import logging

import pandas as pd

from rnacloud_genome_reference.common.gtf import GTFHandler

logger = logging.getLogger(__name__)

def extract_sj_positions_from_clinically_significant_genes(clinical_genes_path: str, gtf_file_path: str, output_path: str) -> None:
    logger.info("Extracting splice junction positions from clinically significant genes...")
    
    clinical_genes = pd.read_csv(clinical_genes_path, sep='\t', low_memory=False)
    
    gtf_file = GTFHandler(gtf_file_path=gtf_file_path)

    with open(output_path, 'w') as f:
        f.write("chrom\tchrom_refseq\tpos\tentrez_gene_id\tgene_name\ttranscript\ttranscript_is_mane_select\texon_no\tdist_from_annot\tcategory\n")

        for _, row in clinical_genes.iterrows():
            sj_positions = gtf_file.obtain_sj_positions(
                chrom=row['chrom_refseq'],
                start=row['start'],
                end=row['end'],
                entrez_gene_id=row['entrez_gene_id']
            )

            for sj_pos in sj_positions:
                f.write("{chrom}\t{chrom_refseq}\t{pos}\t{entrez_gene_id}\t{gene_name}\t{transcript}\t{transcript_is_mane_select}\t{exon_no}\t{dist_from_annot}\t{category}\n".format(
                    chrom=row['chrom'],
                    chrom_refseq=row['chrom_refseq'],
                    pos=sj_pos.pos,
                    entrez_gene_id=row['entrez_gene_id'],
                    gene_name=row['gene_name'],
                    transcript=sj_pos.transcript,
                    transcript_is_mane_select=sj_pos.transcript_is_mane_select,
                    dist_from_annot=sj_pos.dist_from_exon,
                    exon_no=sj_pos.exon_no,
                    category=sj_pos.category
                ))
    
    logger.info("Splice junction positions extraction completed. Output saved to %s", output_path)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract splice junction positions from clinically significant genes.")
    parser.add_argument("--clinical_genes_path", required=True, help="Path to the clinically significant genes file.")
    parser.add_argument("--gtf_file_path", required=True, help="Path to the GTF file.")
    parser.add_argument("--output_path", required=True, help="Path to save the output file.")

    args = parser.parse_args()

    extract_sj_positions_from_clinically_significant_genes(args.clinical_genes_path, args.gtf_file_path, args.output_path)