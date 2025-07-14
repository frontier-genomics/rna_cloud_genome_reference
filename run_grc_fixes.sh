#!/bin/bash
set -e # exit on first error

echo "---------------------------"
echo "Starting"
echo "---------------------------"

echo "Current folder: $(pwd)"

docker run \
  --rm \
  --name rnacloud_runner \
  -v $(pwd)/data:/data \
  -v $(pwd)/temp:/temp \
  -v $(pwd)/output:/output \
  rnacloud_runner -m rnacloud_genome_reference.grc_fixes.assess_grc_fixes -o /output/gene_alt_contigs_mapping_clinically_relevant.tsv