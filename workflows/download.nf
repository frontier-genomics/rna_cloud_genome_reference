// Download workflow processes

nextflow.enable.dsl = 2

process DOWNLOAD_ANNOTATION {
    tag "annotation"
    publishDir "${params.folders.data_dir}/genomes", mode: 'copy'
    
    input:
    val annotation_url
    
    output:
    path "*.gtf.gz", emit: gtf
    
    script:
    def filename = annotation_url.tokenize('/').last()
    """
    echo "Downloading annotation file: ${annotation_url}"
    wget -O ${filename} ${annotation_url}
    """
}

process DOWNLOAD_GENOME {
    tag "genome"
    publishDir "${params.folders.data_dir}/genomes", mode: 'copy'
    
    input:
    val genome_url
    
    output:
    path "*.fna.gz", emit: fasta
    
    script:
    def filename = genome_url.tokenize('/').last()
    """
    echo "Downloading genome file: ${genome_url}"
    wget -O ${filename} ${genome_url}
    """
}

process DOWNLOAD_ASSEMBLY_REPORT {
    tag "assembly_report"
    publishDir "${params.folders.data_dir}/genomes", mode: 'copy'
    
    input:
    val assembly_report_url
    
    output:
    path "*.txt", emit: report
    
    script:
    def filename = assembly_report_url.tokenize('/').last()
    """
    echo "Downloading assembly report: ${assembly_report_url}"
    wget -O ${filename} ${assembly_report_url}
    """
}

process DOWNLOAD_GRC_FIXES {
    tag "grc_fixes"
    publishDir "${params.folders.data_dir}/grc_fixes", mode: 'copy'
    
    input:
    val grc_fixes_url
    
    output:
    path "*.tsv", emit: grc_fixes
    
    script:
    def filename = grc_fixes_url.tokenize('/').last()
    """
    echo "Downloading GRC fixes: ${grc_fixes_url}"
    wget -O ${filename} ${grc_fixes_url}
    """
}

process DOWNLOAD_CLINICALLY_RELEVANT_GENES {
    tag "clinically_relevant_genes"
    publishDir "${params.folders.data_dir}/clinically_relevant_genes", mode: 'copy'
    
    input:
    val clinically_relevant_genes_url
    
    output:
    path "*.tsv", emit: genes
    
    script:
    def filename = clinically_relevant_genes_url.tokenize('/').last()
    """
    echo "Downloading clinically relevant genes: ${clinically_relevant_genes_url}"
    wget -O ${filename} ${clinically_relevant_genes_url}
    """
}