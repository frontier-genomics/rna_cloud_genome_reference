#!/bin/bash
set -e # exit on first error

echo "---------------------------"
echo "Starting"
echo "---------------------------"

echo "Current folder: $(pwd)"

docker run \
  --rm \
  --name rnacloud_runner \
  --entrypoint ./run_splice_site_population_freq.sh \
  -v $(pwd)/data:/data \
  -v $(pwd)/temp:/temp \
  -v $(pwd)/output:/output \
  rnacloud_runner