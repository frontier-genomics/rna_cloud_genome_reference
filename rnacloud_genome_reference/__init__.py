import os
from pathlib import Path
from rnacloud_genome_reference.config import Config

DATA_DIR = Config.get_str('folders', 'data_dir')
TEMP_DIR = Config.get_str('folders', 'temp_dir')
OUTPUT_DIR = Config.get_str('folders', 'output_dir')

ANNOTATION_URL = Config.get_str('genome','annotation')
ANNOTATION_DESTINATION_FOLDER = os.path.join(DATA_DIR, os.path.basename(os.path.dirname(ANNOTATION_URL)))
ANNOTATION_DESTINATION_FILE = os.path.basename(ANNOTATION_URL)
ANNOTATION_DESTINATION_SORTED_FILE = (Path(ANNOTATION_DESTINATION_FILE).with_suffix('')).with_suffix('.sorted.gtf.gz').name

GENOME_URL = Config.get_str('genome', 'fasta')
GENOME_DESTINATION_FOLDER = os.path.join(DATA_DIR, os.path.basename(os.path.dirname(GENOME_URL)))
GENOME_DESTINATION_FILE = os.path.basename(GENOME_URL)

GENOME_REPORT_URL = Config.get_str('genome', 'assembly_report')
GENOME_REPORT_DESTINATION_FOLDER = os.path.join(DATA_DIR, os.path.basename(os.path.dirname(GENOME_REPORT_URL)))
GENOME_REPORT_DESTINATION_FILE = os.path.basename(GENOME_REPORT_URL)

GRC_FIXES_URL = Config.get_str('references', 'grc_fixes')
GRC_FIXES_DESTINATION_FOLDER = os.path.join(DATA_DIR, 'grc_fixes', os.path.basename(os.path.dirname(GRC_FIXES_URL)))
GRC_FIXES_DESTINATION_FILE = os.path.basename(GRC_FIXES_URL)

CLINICALLY_RELEVANT_GENES_URL = Config.get_str('references', 'clinically_relevant_genes')
CLINICALLY_RELEVANT_GENES_DESTINATION_FOLDER = os.path.join(DATA_DIR, 'clinically_relevant_genes', os.path.basename(os.path.dirname(CLINICALLY_RELEVANT_GENES_URL)))
CLINICALLY_RELEVANT_GENES_DESTINATION_FILE = os.path.basename(CLINICALLY_RELEVANT_GENES_URL)

PROTEIN_CODING_GENES = os.path.join(TEMP_DIR, 'protein_coding_genes.tsv')
SIMPLIFIED_GRC_FIXES = os.path.join(TEMP_DIR, 'simplified_grc_fixes.tsv')

GENE_ALT_CONTIGS_MAPPING = os.path.join(TEMP_DIR, 'gene_alt_contigs_mapping.tsv')
GENE_ALT_CONTIGS_COMPARISON = os.path.join(TEMP_DIR, 'gene_alt_contigs_mapping_comparison.tsv')
