#!/bin/bash
set -e # exit on first error

COMBINE_GNOMAD_FREQ_SCRIPT="rnacloud_genome_reference/splice_site_population_freq/scripts/combine_gnomad_freq.sh"
GNOMAD_DATA_PATH=data
GNOMAD_VERSION=gnomad_r4
GNOMAD_REFERENCE_GENOME=GRCh38
GNOMAD_COMBINED_FILE="data/gnomad/${GNOMAD_REFERENCE_GENOME}/${GNOMAD_VERSION}_freq.tsv.gz"

echo "🎬 Starting gnomAD frequency data download and combination script"
echo "  COMBINE_GNOMAD_FREQ_SCRIPT : ${COMBINE_GNOMAD_FREQ_SCRIPT}"
echo "  GNOMAD_DATA_PATH           : ${GNOMAD_DATA_PATH}"
echo "  GNOMAD_VERSION             : ${GNOMAD_VERSION}"
echo "  GNOMAD_REFERENCE_GENOME    : ${GNOMAD_REFERENCE_GENOME}"
echo "  GNOMAD_COMBINED_FILE       : ${GNOMAD_COMBINED_FILE}"

if [ -f ${GNOMAD_COMBINED_FILE} ]; then
    echo "✅ Combined gnomAD frequency file already exists."
else
    echo 🏃‍♂️ Starting gnomad frequency data download
    python3 -m rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq --chunk-size 100000
    echo 🏁 Finished download

    echo "🏃‍♂️ Combining gnomAD frequency data into a single file."
    ${COMBINE_GNOMAD_FREQ_SCRIPT} ${GNOMAD_DATA_PATH} ${GNOMAD_COMBINED_FILE}
    echo "🏁 Finished combining gnomAD frequency data"
fi

echo "🏁 Finished gnomAD frequency data download and combination script"