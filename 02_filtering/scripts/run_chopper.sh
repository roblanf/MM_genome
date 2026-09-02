#!/usr/bin/env bash
# Script: 02_filtering/scripts/run_chopper.sh
# Purpose: Filter reads by length and quality score using Chopper
# Usage: bash 02_filtering/scripts/run_chopper.sh [input_fastq] [output_fastq] [min_length] [min_qual] [threads]

set -euo pipefail

INPUT_FASTQ="${1:-02_filtering/results/tiny_test_porechop.fastq.gz}"
OUTPUT_FASTQ="${2:-02_filtering/results/tiny_test_filtered.fastq.gz}"
MIN_LENGTH="${3:-1000}"   # Minimum read length (bp)
MIN_QUAL="${4:-10}"       # Minimum average Phred quality score (Q10)
THREADS="${5:-4}"

# Ensure output directory exists
mkdir -p "$(dirname "${OUTPUT_FASTQ}")"

echo "============================================================"
echo "[Chopper] Starting Quality & Length Filtering"
echo "Input:      ${INPUT_FASTQ}"
echo "Output:     ${OUTPUT_FASTQ}"
echo "Min Length: ${MIN_LENGTH} bp"
echo "Min Qual:   Q${MIN_QUAL}"
echo "Threads:    ${THREADS}"
echo "============================================================"

# Process gzipped or uncompressed FASTQ directly
if [[ "${INPUT_FASTQ}" == *.gz ]]; then
  zcat "${INPUT_FASTQ}" | chopper -q "${MIN_QUAL}" -l "${MIN_LENGTH}" -t "${THREADS}" | gzip > "${OUTPUT_FASTQ}"
else
  chopper -i "${INPUT_FASTQ}" -q "${MIN_QUAL}" -l "${MIN_LENGTH}" -t "${THREADS}" | gzip > "${OUTPUT_FASTQ}"
fi

echo "============================================================"
echo "[Chopper] Finished successfully."
echo "============================================================"
