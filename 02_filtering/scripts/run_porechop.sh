#!/usr/bin/env bash
# Script: 02_filtering/scripts/run_porechop.sh
# Purpose: Trim Nanopore adapters using Porechop - it discovers unkown adapters in ONT reads for downstream trimming
# Usage: bash 02_filtering/scripts/run_porechop.sh [input_fastq] [output_fastq] [threads]

set -euo pipefail

INPUT_FASTQ="${1:-raw_data/tiny_test.fastq}"
OUTPUT_FASTQ="${2:-02_filtering/results/tiny_test_porechop.fastq.gz}"
THREADS="${3:-4}"

# Ensure output directory exists
mkdir -p "$(dirname "${OUTPUT_FASTQ}")"

echo "============================================================"
echo "[Porechop] Starting Adapter Trimming"
echo "Input:   ${INPUT_FASTQ}"
echo "Output:  ${OUTPUT_FASTQ}"
echo "Threads: ${THREADS}"
echo "============================================================"

# Porechop outputs to stdout when redirected or directly to file
porechop \
  -i "${INPUT_FASTQ}" \
  -o "${OUTPUT_FASTQ%.gz}" \
  --threads "${THREADS}" \
  --verbosity 1

# Compress intermediate if output ends in .gz
if [[ "${OUTPUT_FASTQ}" == *.gz ]]; then
  gzip -f "${OUTPUT_FASTQ%.gz}"
fi

echo "============================================================"
echo "[Porechop] Finished successfully."
echo "============================================================"
