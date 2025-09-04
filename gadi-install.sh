#!/bin/bash

## GITHUB CONFIGURATION
GITHUB_OWNER="kidsneuro-lab"
GITHUB_REPO="rna_cloud_genome_reference"
VERSION="0.0.12"

## GADI CONFIGURATION
BASE_DIR=/g/data/qe93/rna_cloud_pipeline
CUSTOM_ASSEMBLY="${BASE_DIR}/GRCh38/assemblies/custom/rna_cloud/$VERSION"
CUSTOM_ANNOTATION="${BASE_DIR}/GRCh38/annotations/custom/rna_cloud/$VERSION"
rm -rf ${CUSTOM_ASSEMBLY} && mkdir -p ${CUSTOM_ASSEMBLY}
rm -rf ${CUSTOM_ANNOTATION} && mkdir -p ${CUSTOM_ANNOTATION}

echo "🔵 Downloading custom assembly to ${CUSTOM_ASSEMBLY}"
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".fasta.gz"      $CUSTOM_ASSEMBLY
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".fasta.gz.fai"  $CUSTOM_ASSEMBLY
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".fasta.gz.gzi"  $CUSTOM_ASSEMBLY

echo "🔵 Downloading custom annotation to ${CUSTOM_ANNOTATION}"
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".gtf.gz"      $CUSTOM_ANNOTATION
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".gtf.gz.tbi"  $CUSTOM_ANNOTATION