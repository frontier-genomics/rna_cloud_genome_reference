nextflow.enable.dsl=2

process CALCULATE_MD5_SUMMARY {
    tag "CALCULATE_MD5_SUMMARY"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    // 'files' is a list of file paths collected from final outputs of subworkflows
    path files

    output:
    // Output a single file with all md5 checksums
    path 'md5checksums.txt'

    script:
    """
    set -euo pipefail
    
    # Compute md5 checksums for all input files
    # 'files.collect { it.name }' extracts just the filenames from the full paths
    # 'join(" ")' turns the list of filenames into a single space-separated string
    # The resulting string is passed to 'md5sum', which computes checksums
    # The output is redirected into 'md5checksums.txt'
    md5sum ${files.collect { it.name }.join(' ')} > md5checksums.txt
    """
}