#!/bin/bash

set -euo pipefail

log() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $*"
}

error_exit() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] ERROR: $*" >&2
    exit 1
}

if [ $# -ne 2 ]; then
    error_exit "Usage: $0 <input_folder> <output_file.gz>"
fi

input_folder="$1"
output_file="$2"

if [ ! -d "$input_folder" ]; then
    error_exit "Input folder '$input_folder' does not exist or is not a directory."
fi

shopt -s nullglob
files=("$input_folder"/*.tsv.gz)
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
    error_exit "No .tsv.gz files found in '$input_folder'."
fi

log "Combining ${#files[@]} files from '$input_folder' into '$output_file'..."

# Write header from the first file
header=$(gunzip -c "${files[0]}" | head -n 1)
if [ -z "$header" ]; then
    error_exit "Failed to read header from the first file '${files[0]}'."
fi
log "Writing header to '$output_file'..."
echo "$header" | gzip > "$output_file"

# Append the contents of all files except their headers
for file in "${files[@]}"; do
    log "Processing $file..."
    # Skip the first line (header) of every file
    gunzip -c "$file" | tail -n +2 | gzip >> "$output_file"
done

log "Successfully combined files into '$output_file'."