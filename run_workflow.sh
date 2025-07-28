#!/bin/bash
set -e # exit on first error

echo "----------------------------------------------------------------------"
echo "Starting rnacloud_genome_reference GRC fixes assessment"
echo "----------------------------------------------------------------------"

echo "Current folder: $(pwd)"

docker run \
  --rm \
  --name rnacloud_runner \
  --entrypoint nextflow \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/temp:/app/temp \
  -v $(pwd)/output:/app/output \
  -v $(pwd)/work:/app/work \
  rnacloud_runner main.nf -o /app/output/grc_fixes_assessment.tsv