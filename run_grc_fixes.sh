#!/bin/bash
set -e # exit on first error

echo "----------------------------------------------------------------------"
echo "Starting rnacloud_genome_reference GRC fixes assessment"
echo "----------------------------------------------------------------------"

echo "Current folder: $(pwd)"

docker run \
  --rm \
  --name rnacloud_runner \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/temp:/app/temp \
  -v $(pwd)/output:/app/output \
  rnacloud_runner -m rnacloud_genome_reference.grc_fixes.assess_grc_fixes -o /app/output/gene_alt_contigs_mapping_clinically_relevant.tsv