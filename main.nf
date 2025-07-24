#!/usr/bin/env nextflow

/*
 * RNACloud Genome Reference Workflow
 * 
 * This workflow processes genome reference data for the RNACloud project.
 * It includes GRC fixes assessment and splice site population frequency analysis.
 */

nextflow.enable.dsl = 2

// Import workflows
include { DOWNLOAD_ANNOTATION } from './workflows/download.nf'
include { DOWNLOAD_GENOME } from './workflows/download.nf'
include { DOWNLOAD_ASSEMBLY_REPORT } from './workflows/download.nf'
include { DOWNLOAD_GRC_FIXES } from './workflows/download.nf'
include { DOWNLOAD_CLINICALLY_RELEVANT_GENES } from './workflows/download.nf'
include { ASSESS_GRC_FIXES } from './workflows/grc_fixes.nf'
include { SPLICE_SITE_POPULATION_FREQ } from './workflows/splice_site.nf'

// Main workflow
workflow {
    
    // Define input parameters from config
    annotation_url = params.genome.annotation
    genome_url = params.genome.fasta
    assembly_report_url = params.genome.assembly_report
    grc_fixes_url = params.references.grc_fixes
    clinically_relevant_genes_url = params.references.clinically_relevant_genes
    
    // Download all required files
    DOWNLOAD_ANNOTATION(annotation_url)
    DOWNLOAD_GENOME(genome_url)
    DOWNLOAD_ASSEMBLY_REPORT(assembly_report_url)
    DOWNLOAD_GRC_FIXES(grc_fixes_url)
    DOWNLOAD_CLINICALLY_RELEVANT_GENES(clinically_relevant_genes_url)
    
    // Run GRC fixes assessment workflow
    grc_fixes_result = ASSESS_GRC_FIXES(
        DOWNLOAD_ANNOTATION.out.gtf,
        DOWNLOAD_GENOME.out.fasta,
        DOWNLOAD_ASSEMBLY_REPORT.out.report,
        DOWNLOAD_GRC_FIXES.out.grc_fixes,
        DOWNLOAD_CLINICALLY_RELEVANT_GENES.out.genes
    )
    
    // Run splice site population frequency workflow
    splice_site_result = SPLICE_SITE_POPULATION_FREQ(
        DOWNLOAD_ANNOTATION.out.gtf,
        DOWNLOAD_CLINICALLY_RELEVANT_GENES.out.genes
    )
    
    // Emit final outputs
    emit:
    grc_fixes_output = grc_fixes_result
    splice_site_output = splice_site_result
}

// Workflow completion notification
workflow.onComplete {
    println "RNACloud Genome Reference workflow completed!"
    println "Results available in: ${params.folders.output_dir}"
}