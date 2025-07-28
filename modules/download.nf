process DOWNLOAD_AND_INDEX_GENOME {
    tag "${genome_fasta_url.tokenize('/')[-1]}"
    storeDir "${params.data_dir}"

    input:
    val genome_fasta_url

    output:
    path "${genome_fasta_url.tokenize('/')[-1]}", emit: fasta
    path "${genome_fasta_url.tokenize('/')[-1]}.fai", emit: fasta_fai_index
    path "${genome_fasta_url.tokenize('/')[-1]}.gzi", emit: fasta_gzi_index

    script:
    """
    set -euo pipefail

    file_gz=\$(basename ${genome_fasta_url})
    file_fa=\$(basename \${file_gz} .gz)

    echo "Downloading genome from: ${genome_fasta_url}"
    wget -qc ${genome_fasta_url}

    echo "Unzipping \${file_gz}..."
    gunzip -f \${file_gz}

    echo "Compressing with bgzip..."
    bgzip \${file_fa}

    echo "Indexing with samtools..."
    samtools faidx \${file_gz}
    """
}

process DOWNLOAD_AND_INDEX_GTF {
    tag "${gtf_file_url.tokenize('/')[-1]}"
    storeDir "${params.data_dir}"

    input:
    val gtf_file_url

    output:
    path "${gtf_file_url.tokenize('/')[-1].replace('.gtf.gz', '.sorted.gtf.gz')}", emit: gtf
    path "${gtf_file_url.tokenize('/')[-1].replace('.gtf.gz', '.sorted.gtf.gz.tbi')}", emit: gtf_index

    script:
    """
    set -euo pipefail

    file_gz=\$(basename ${gtf_file_url})
    file_basename=\$(basename \${file_gz} .gtf.gz)

    echo "Downloading GTF from: ${gtf_file_url}"
    wget -qc ${gtf_file_url}

    echo "Sorting GTF file..."
    bedtools sort -i \${file_gz} | bgzip > \${file_basename}.sorted.gtf.gz

    echo "Indexing GTF file..."
    tabix -p gff \${file_basename}.sorted.gtf.gz
    """
}

process DOWNLOAD_FILE {
    tag "${file_url.tokenize('/')[-1]}"
    storeDir "${params.data_dir}"

    input:
    val file_url

    output:
    path "${file_url.tokenize('/')[-1]}", emit: downloaded_file

    script:
    """
    set -euo pipefail

    echo "Downloading file from: ${file_url}"
    wget -qc ${file_url}
    """
}