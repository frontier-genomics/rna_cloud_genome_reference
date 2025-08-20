#!/bin/bash

set -uo pipefail

if [ $# -ne 6 ]; then
  echo "Usage: $0 <FASTA> <FASTA_INDEX> <GTF> <GTF_INDEX> <MASKED_REGIONS_BED> <UNMASKED_REGIONS_BED>"
  exit 1
fi

FASTA="$1"
FASTA_INDEX="$2"
GTF="$3"
GTF_INDEX="$4"
MASKED_REGIONS_BED="$5"
UNMASKED_REGIONS_BED="$6"

echo "FASTA                 : $FASTA"
echo "FASTA Index           : $FASTA_INDEX"
echo "GTF                   : $GTF"
echo "GTF Index             : $GTF_INDEX"
echo "Masked Regions BED    : $MASKED_REGIONS_BED"
echo "Unmasked Regions BED  : $UNMASKED_REGIONS_BED"

# Initialize validation status
VALIDATION_STATUS=0
echo 🎬 Starting Genome and Annotation build validation process

echo 🕵️‍♀️ Comparing chromosomes in GTF and FASTA files...

GTF_CONTIGS=$(mktemp --suffix=.txt)
FASTA_CONTIGS=$(mktemp --suffix=.txt)

tabix -l "$GTF" | sort > "$GTF_CONTIGS"
cut -f 1 "$FASTA_INDEX" | grep -v 'chrEBV' | sort > "$FASTA_CONTIGS"

diff "$GTF_CONTIGS" "$FASTA_CONTIGS"

if [ $? -eq 0 ]; then
  echo ✅ Contigs in GTF and FASTA match!
else
  echo ⛔️ Contigs in GTF and FASTA do not match!
  VALIDATION_STATUS=1
fi

rm -f "$GTF_CONTIGS" "$FASTA_CONTIGS"

echo 🕵️‍♀️ Check that EBV contig is present in the FASTA file...

EBV_CONTIG_COUNT=$(cut -f 1 "$FASTA_INDEX" | grep 'chrEBV' | wc -l)

if [ "$EBV_CONTIG_COUNT" -eq 1 ]; then
  echo ✅ EBV contig is present in the FASTA file!
else
  echo ⛔️ EBV contig is not present in the FASTA file or contains multiple entries \(No. of EBV contigs: ${EBV_CONTIG_COUNT}\)!
  VALIDATION_STATUS=1
fi

echo 🕵️‍♀️ Check that masked regions only contain Ns

bedtools getfasta \
  -fi "$FASTA" \
  -bed "$MASKED_REGIONS_BED" \
  -tab \
  -fo masked_regions.fasta 2>&1 | tee masked_regions.fasta.log

ATCG_COUNT=$(awk -F"\t" '$2~/[ATCG]/' masked_regions.fasta | wc -l)

if [ "$ATCG_COUNT" -eq 0 ]; then
  echo ✅ Masked regions only contain Ns!
else
  echo ⛔️ Masked regions contain non-N bases \(Contigs with ATCG bases: ${ATCG_COUNT}\)!
  VALIDATION_STATUS=1
fi

WARNINGS_COUNT=$(grep 'WARNING' masked_regions.fasta.log | wc -l)

if [ "$WARNINGS_COUNT" -eq 0 ]; then
  echo ✅ No warnings found in masked regions FASTA extraction!
else
  echo ⛔️ Warnings found in masked regions FASTA extraction \(No. of warnings: ${WARNINGS_COUNT}\)!
  VALIDATION_STATUS=1
fi

echo 🕵️‍♀️ Check that unmasked regions only do not contain Ns

bedtools getfasta \
  -fi "$FASTA" \
  -bed "$UNMASKED_REGIONS_BED" \
  -tab \
  -fo unmasked_regions.fasta 2>&1 | tee unmasked_regions.fasta.log

N_COUNT=$(awk -F"\t" '$2~/N/' unmasked_regions.fasta | wc -l)

if [ "$N_COUNT" -eq 0 ]; then
  echo ✅ Unmasked regions only contain ATCGs!
else
  echo ⛔️ Unmasked regions contain N bases \(Regions with N bases: ${N_COUNT}\)!
  VALIDATION_STATUS=1
fi

WARNINGS_COUNT=$(grep 'WARNING' output/unmasked_regions.fasta.log | wc -l)

if [ "$WARNINGS_COUNT" -eq 0 ]; then
  echo ✅ No warnings found in unmasked regions FASTA extraction!
else
  echo ⛔️ Warnings found in unmasked regions FASTA extraction \(No. of warnings: ${WARNINGS_COUNT}\)!
  VALIDATION_STATUS=1
fi

if [ "$VALIDATION_STATUS" -ne 0 ]; then
  echo ❗ Validation failed!
  exit "$VALIDATION_STATUS"
fi

echo 🎉 All validations passed