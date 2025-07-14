#!/bin/bash
set -e # exit on first error

echo "---------------------------"
echo "Starting"
echo "---------------------------"

echo "Current folder: $(pwd)"

docker run \
  --rm \
  --name rnacloud_runner \
  --entrypoint ./rnacloud_genome_reference/splice_site_population_freq/scripts/sj_gnomad_freq.sh \
  -v $(pwd)/data:/data \
  -v $(pwd)/temp:/temp \
  -v $(pwd)/output:/output \
  rnacloud_runner