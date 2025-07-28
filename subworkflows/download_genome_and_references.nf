nextflow.enable.dsl=2

//
// Pull in all the DOWNLOAD_* processes (path is relative to this file)
//
include { DOWNLOAD_AND_INDEX_GENOME }   from '../modules/download.nf'
include { DOWNLOAD_AND_INDEX_GTF }      from '../modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_ASSEMBLY_REPORT }            from '../modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_GRC_FIXES }                  from '../modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_CLINICALLY_RELEVANT_GENES }  from '../modules/download.nf'

//
// A DSL2 sub-workflow that reads params.* directly
//
workflow DOWNLOAD_GENOME_AND_REFERENCES {
    main:
    println "Downloading genome and references with the following parameters:"
    println "Genome FASTA URL              : ${params.genome_fasta_url}"
    println "Genome annotation URL         : ${params.genome_annotation_url}"
    println "Genome assembly report URL    : ${params.genome_assembly_report}"
    println "Reference GRC fixes URL       : ${params.reference_grc_fixes}"
    println "Clinically relevant genes URL : ${params.reference_clinically_relevant_genes}"

    // build channels from params inside the sub-workflow
    DOWNLOAD_AND_INDEX_GENOME(Channel.from(params.genome_fasta_url))
    DOWNLOAD_AND_INDEX_GTF(Channel.from(params.genome_annotation_url))
    DOWNLOAD_ASSEMBLY_REPORT(Channel.from(params.genome_assembly_report))
    DOWNLOAD_GRC_FIXES(Channel.from(params.reference_grc_fixes))
    DOWNLOAD_CLINICALLY_RELEVANT_GENES(Channel.from(params.reference_clinically_relevant_genes))

    emit:
    // expose the exact same seven channels as before
    fasta                 = DOWNLOAD_AND_INDEX_GENOME.out.fasta
    fasta_index           = DOWNLOAD_AND_INDEX_GENOME.out.fasta_fai_index
    gtf                   = DOWNLOAD_AND_INDEX_GTF.out.gtf
    gtf_index             = DOWNLOAD_AND_INDEX_GTF.out.gtf_index
    assembly_report       = DOWNLOAD_ASSEMBLY_REPORT.out.downloaded_file
    grc_fixes             = DOWNLOAD_GRC_FIXES.out.downloaded_file
    clinically_relevant   = DOWNLOAD_CLINICALLY_RELEVANT_GENES.out.downloaded_file
}