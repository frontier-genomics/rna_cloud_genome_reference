#!/bin/bash

set -e # exit on first error

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

if [ -d "/app/temp" ]; then
    TEMP_FOLDER="/app/temp"
    log "Using existing temp folder: $TEMP_FOLDER"
elif [ -d "/tmp" ]; then
    TEMP_FOLDER="/tmp"
    log "Using existing temp folder: $TEMP_FOLDER"
fi

input_folder="$1"
output_file="$2"
temp_file="${output_file}.tmp"

if [ ! -d "$input_folder" ]; then
    error_exit "Input folder '$input_folder' does not exist or is not a directory."
fi

shopt -s nullglob
files=("$input_folder"/*.tsv.gz)
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
    error_exit "No .tsv.gz files found in '$input_folder'."
fi

log "Combining ${#files[@]} files from '$input_folder' into '$temp_file'..."

# Append the contents of all files except their headers
for file in "${files[@]}"; do
    log "Processing $file..."
    # Skip the first line (header) of every file
    gunzip -c "$file" | tail -n +2 | gzip >> "$temp_file"
done

# Write header from the first file
header=$(gunzip -c "${files[0]}" | head -n 1)
if [ -z "$header" ]; then
    error_exit "Failed to read header from the first file '${files[0]}'."
fi
log "Writing header to '$output_file'..."
echo "$header" | gzip > "$output_file"

# Sort the temp file and remove duplicates
log "Sorting and removing duplicates..."
gunzip -c "$temp_file" | sort -u -T "$TEMP_FOLDER" | gzip >> "$output_file"
rm "$temp_file"

log "Successfully combined files into '$output_file'."