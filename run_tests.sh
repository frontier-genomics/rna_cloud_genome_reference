#!/bin/bash
set -e # exit on first error

echo "---------------------------"
echo "Running tests"
echo "---------------------------"

docker run \
  --rm \
  --name rnacloud_runner \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/temp:/app/temp \
  -v $(pwd)/output:/app/output \
  rnacloud_runner -m pytest -v