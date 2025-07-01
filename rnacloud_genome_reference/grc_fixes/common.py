import logging
import os
import subprocess
import urllib.request
from pathlib import Path
import sqlite3

import pandas as pd

logger = logging.getLogger(__name__)

def sort_index_gtf_file(file_path: str,
                        sorted_file_path: str):
    
    gtf_file = Path(file_path)
    sorted_gtf_file = Path(sorted_file_path)

    if sorted_gtf_file.with_suffix('.gz.tbi').exists():
        logger.info(f"Sorted and indexed file already exists: {sorted_file_path}")
        return
    
    try:
        logger.info(f"Processing gtf file: {file_path}")
            
        # Step 1: Sort GTF file using bedtools sort
        logger.info("Sorting GTF file with bedtools sort...")
        sort_cmd = f"bedtools sort -i {gtf_file} | bgzip > {sorted_gtf_file}"
        result = subprocess.run(
            sort_cmd,
            shell=True,
            check=True, 
            capture_output=True, 
            text=True
        )

        logger.info(f"Indexing sorted GTF file with tabix...")
        result = subprocess.run(
            ['tabix', str(sorted_gtf_file)], 
            check=True, 
            capture_output=True, 
            text=True
        )

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e.cmd}")
        logger.error(f"Return code: {e.returncode}")
        logger.error(f"Error output: {e.stderr}")
        raise
        
    except FileNotFoundError as e:
        logger.error(f"Command not found. Please ensure the required tools are installed:")
        logger.error("- samtools")
        logger.error("- gunzip")
        logger.error("- bgzip")
        raise
        
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise
        

def index_fasta_file(file_path: str):
    """
    Index a FASTA file using samtools faidx.
    
    If the file is .gz compressed, it will:
    1. Decompress with gunzip
    2. Recompress with bgzip
    3. Index with samtools faidx
    
    Args:
        file_path (str): Path to the FASTA file to index
        
    Returns:
        bool: True if indexing successful, False otherwise
        
    Raises:
        FileNotFoundError: If the input file doesn't exist
        subprocess.CalledProcessError: If any command fails
    """
    
    # Convert to Path object for easier manipulation
    file_path = Path(file_path)

    index_file = file_path.with_suffix(file_path.suffix + '.fai')
    if index_file.exists():
        logger.info(f"Index file already exists: {index_file}")
        return
    
    # Check if file exists
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    try:
        # Check if file is gzipped
        if file_path.suffix == '.gz':
            logger.info(f"Processing gzipped file: {file_path}")
            
            # Step 1: Decompress with gunzip
            logger.info("Decompressing with gunzip...")
            result = subprocess.run(
                ['gunzip', str(file_path)], 
                check=True, 
                capture_output=True, 
                text=True
            )
            
            # Update file_path to point to decompressed file
            decompressed_file = file_path.with_suffix('')
            
            # Step 2: Recompress with bgzip
            logger.info("Recompressing with bgzip...")
            result = subprocess.run(
                ['bgzip', str(decompressed_file)], 
                check=True, 
                capture_output=True, 
                text=True
            )
        
        # Step 3: Index with samtools faidx
        logger.info(f"Indexing file with samtools faidx: {file_path}")
        result = subprocess.run(
            ['samtools', 'faidx', str(file_path)], 
            check=True, 
            capture_output=True, 
            text=True
        )
        
        # Check if index file was created
        index_file = file_path.with_suffix(file_path.suffix + '.fai')
        if index_file.exists():
            logger.info(f"Successfully created index: {index_file}")
            return
        else:
            logger.error("Index file was not created")
            raise RuntimeError("Index file was not created")
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e.cmd}")
        logger.error(f"Return code: {e.returncode}")
        logger.error(f"Error output: {e.stderr}")
        raise
        
    except FileNotFoundError as e:
        logger.error(f"Command not found. Please ensure the required tools are installed:")
        logger.error("- samtools")
        logger.error("- gunzip")
        logger.error("- bgzip")
        raise
        
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise


def download_file(url: str, destination_folder: str, destination_file: str) -> None:
    # Check if destination file already exists, if so - log a warning
    file_path = os.path.join(destination_folder, destination_file)
    if os.path.exists(file_path):
        logger.warning(f"File {file_path} already exists, skipping...")
        return
    
    # Check if folder exists - if not create it
    Path(destination_folder).mkdir(parents=True, exist_ok=True)
    
    # Download the file - Throw an error if there is any error downloading
    logger.info(f"Starting download of {url} to {file_path}")
    try:
        urllib.request.urlretrieve(url, file_path)
        logger.info(f"Download completed: {file_path}")
    except Exception as e:
        logger.error(f"Failed to download {url}: {e}")
        raise


def combine_grc_fixes_and_protein_coding_genes(grc_fixes_file: str, protein_coding_genes_file: str, output_file: str) -> None:
    # Read GRC fixes and protein coding genes files
    grc_fixes_df = pd.read_csv(grc_fixes_file, sep='\t')
    protein_coding_genes_df = pd.read_csv(protein_coding_genes_file, sep='\t')

    # Write both to sqlite database
    conn = sqlite3.connect(":memory:")
    grc_fixes_df.to_sql('grc_fixes', conn, if_exists='replace', index=False)
    protein_coding_genes_df.to_sql('protein_coding_genes', conn, if_exists='replace', index=False)

    sql = """
        SELECT gf.chr_refseq,
        gf.chr_ucsc,
        pcg.start,
        pcg.end,
        pcg.strand,
        pcg.gene_name,
        pcg.entrez_gene_id,
        gf.issue_id,
        gf.alt_chr_refseq,
        gf.alt_chr_ucsc,
        gf.alt_scaf_start,
        gf.alt_scaf_stop
    FROM grc_fixes gf
    JOIN protein_coding_genes pcg
    ON pcg.end >= gf.parent_start
    AND pcg.start <= gf.parent_stop
    AND pcg.chr = gf.chr_refseq;
    """

    combined_df = pd.read_sql_query(sql, conn)
    combined_df.to_csv(output_file, sep='\t', index=False, header=True)
    logger.info(f"Combined GRC fixes and protein coding genes saved to {output_file}")

    conn.close()