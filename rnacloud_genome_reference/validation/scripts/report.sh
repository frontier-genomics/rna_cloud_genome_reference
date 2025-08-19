#!/bin/bash

set -uo pipefail

FASTA="output/GCF_000001405.40_GRCh38.p14_rna_cloud.fasta.gz"
FASTA_INDEX="output/GCF_000001405.40_GRCh38.p14_rna_cloud.fasta.gz.fai"
GTF="output/GCF_000001405.40_GRCh38.p14_rna_cloud.gtf.gz"
GTF_INDEX="output/GCF_000001405.40_GRCh38.p14_rna_cloud.gtf.gz.tbi"
MASKED_REGIONS_BED="output/mask_regions.bed"


echo ➕ "Number of masked regions: $(wc -l < $MASKED_REGIONS_BED)"