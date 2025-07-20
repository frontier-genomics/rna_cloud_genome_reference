#!/bin/bash
set -e # exit on first error

echo "------------------------------------------------------------------"
echo "Starting splice site population frequency derivation"
echo "------------------------------------------------------------------"

echo "Current folder: $(pwd)"

echo "Clearing temp and output directories"
rm -rf temp/*.tsv
rm -rf output/*.tsv

docker run \
  --rm \
  --name rnacloud_runner \
  --entrypoint ./rnacloud_genome_reference/splice_site_population_freq/scripts/splice_site_gnomad_freq.sh \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/temp:/app/temp \
  -v $(pwd)/output:/app/output \
  rnacloud_runner