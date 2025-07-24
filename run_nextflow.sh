#!/bin/bash
set -e # exit on first error

echo "----------------------------------------------------------------------"
echo "Starting RNACloud Genome Reference Nextflow workflow"
echo "----------------------------------------------------------------------"

echo "Current folder: $(pwd)"

# Check if Nextflow is available
if ! command -v nextflow >/dev/null 2>&1; then
    echo "Nextflow not found. Installing Nextflow..."
    curl -s https://get.nextflow.io | bash
    export PATH=$PATH:$(pwd)
fi

# Run the Nextflow workflow
echo "Running Nextflow workflow..."
nextflow run main.nf \
    --profile docker \
    -with-timeline timeline.html \
    -with-report report.html \
    -with-trace trace.txt \
    -with-dag dag.svg

echo "----------------------------------------------------------------------"
echo "Nextflow workflow completed!"
echo "----------------------------------------------------------------------"