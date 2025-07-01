#!/bin/bash
set -e # exit on first error

echo "---------------------------"
echo "Running tests"
echo "---------------------------"

docker run \
  --rm \
  --name rnacloud_runner \
  -v $(pwd)/data:/data \
  -v $(pwd)/temp:/temp \
  -v $(pwd)/output:/output \
  rnacloud_runner -m pytest -v