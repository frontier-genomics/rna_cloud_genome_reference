// GRC fixes assessment workflow

nextflow.enable.dsl = 2

process ASSESS_GRC_FIXES {
    tag "assess_grc_fixes"
    publishDir "${params.folders.output_dir}", mode: 'copy'
    
    input:
    path gtf_file
    path genome_file
    path assembly_report
    path grc_fixes_file
    path clinically_relevant_genes
    
    output:
    path "gene_alt_contigs_mapping_clinically_relevant.tsv", emit: assessment
    
    script:
    def gtf_name = gtf_file.getName()
    def genome_name = genome_file.getName()
    def report_name = assembly_report.getName()
    def grc_name = grc_fixes_file.getName()
    def genes_name = clinically_relevant_genes.getName()
    """
    # Create necessary directory structure
    mkdir -p ${params.folders.data_dir}/genomes
    mkdir -p ${params.folders.data_dir}/grc_fixes/\$(dirname "${grc_name}")
    mkdir -p ${params.folders.data_dir}/clinically_relevant_genes/\$(dirname "${genes_name}")
    mkdir -p ${params.folders.temp_dir}
    mkdir -p ${params.folders.output_dir}
    
    # Copy input files to expected locations to match the original structure
    cp ${gtf_file} ${params.folders.data_dir}/genomes/${gtf_name}
    cp ${genome_file} ${params.folders.data_dir}/genomes/${genome_name}
    cp ${assembly_report} ${params.folders.data_dir}/genomes/${report_name}
    cp ${grc_fixes_file} ${params.folders.data_dir}/grc_fixes/${grc_name}
    cp ${clinically_relevant_genes} ${params.folders.data_dir}/clinically_relevant_genes/${genes_name}
    
    # Run the assessment using the existing Python module
    python -m rnacloud_genome_reference.grc_fixes.assess_grc_fixes \\
        -o gene_alt_contigs_mapping_clinically_relevant.tsv
    """
}