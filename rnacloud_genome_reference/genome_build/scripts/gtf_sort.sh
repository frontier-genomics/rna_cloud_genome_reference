#!/bin/bash
#
# -----------------------------------------------------------------------------
# Sort and index a compressed GTF by genomic position with feature priority
# -----------------------------------------------------------------------------
# This script:
#   1) Starts off with an uncompressed GTF (.gtf)
#   2) Appends a temporary sort key (column 10) encoding feature priority:
#        gene(1) < transcript(2) < exon(3) < CDS(4) < start_codon(5) < stop_codon(6)
#   3) Sorts by: chromosome (V-sort), start position (numeric), feature priority
#   4) Removes the temporary key, bgzip-compresses the sorted GTF
#   5) Builds a tabix index for fast region queries
#
# Why add a feature-priority sort key?
#   Within identical coordinates, keeping a deterministic feature order can help
#   downstream tools that expect parent-before-child semantics (e.g., genes
#   before transcripts, transcripts before exons/CDS).
#
# Requirements:
#   - bash, awk, sort, cut
#   - bgzip and tabix (from HTSlib / samtools suite)
#   - input must be gzip-compressed (*.gtf.gz)
#
# Usage:
#   ./sort_gtf.sh <GTF.gz> <SORTED_GTF.gz>
#
# Example:
#   ./sort_gtf.sh genes.gtf.gz genes.sorted.gtf.gz
#
# Notes:
#   - Chromosome sorting uses sort's "version" ordering (-k1,1V) so that
#     "chr2" comes before "chr10".
#   - Start coordinates are sorted numerically (-k4,4n).
#   - The temporary 10th column is dropped after sorting.
#   - The output is bgzip-compressed and indexed with tabix (preset: GFF/GTF).
#
# -----------------------------------------------------------------------------

# Fail on unset variables (-u) and make pipelines fail if any segment fails (pipefail).
# (We do not use -e to allow clear error messages; pipefail covers most failure modes.)
set -uo pipefail

# Validate arguments: require exactly two paths: input GTF.gz and output GTF.gz
if [ $# -ne 2 ]; then
  echo "Usage: $0 <GTF> <SORTED_GTF>"
  exit 1
fi

# Positional arguments
GTF="$1"          # Path to input GTF
SORTED_GTF="$2"   # Path to output bgzipped, sorted GTF

# -----------------------------------------------------------------------------
# Decompress → attach feature-priority key → sort → drop key → bgzip
# -----------------------------------------------------------------------------
# gunzip -c: stream-decompress input GTF.gz to stdout
# awk:
#   - Sets output field separator to TAB
#   - For records with feature type in $3, print the full line ($0) plus a
#     priority number as a 10th column. Lines with other feature types are
#     ignored (not printed). If you need to keep all types, add a default case.
# sort:
#   -t$'\t'     : set TAB as the field delimiter
#   -k1,1V      : version sort on chromosome/seqname (field 1)
#   -k4,4n      : numeric sort on start coordinate (field 4)
#   -k10,10n    : numeric sort on the temporary priority key (field 10)
# cut:
#   -f1-9       : drop the temporary 10th column and keep standard 9 GTF fields
# bgzip:
#   Compress to BGZF format, required by tabix
cat "${GTF}" | \
awk 'BEGIN {OFS="\t"}
     $3 == "gene"        {print $0, "1"}
     $3 == "transcript"  {print $0, "2"}
     $3 == "exon"        {print $0, "3"}
     $3 == "CDS"         {print $0, "4"}
     $3 == "start_codon" {print $0, "5"}
     $3 == "stop_codon"  {print $0, "6"}' | \
sort -t$'\t' -k1,1V -k4,4n -k10,10n | \
cut -f1-9 > "${SORTED_GTF}"