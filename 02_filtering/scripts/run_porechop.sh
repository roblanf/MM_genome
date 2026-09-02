#!/usr/bin/env bash
# Script: 02_filtering/scripts/run_porechop.sh
# Purpose: Trim Oxford Nanopore adapters (start/end) and split chimeric reads with Porechop.
# Usage: bash 02_filtering/scripts/run_porechop.sh [input_fastq] [output_dir] [threads]

set -euo pipefail

# Inputs and defaults
INPUT_FASTQ="${1:-raw_data/raw_reads.fastq.gz}"
OUTDIR="${2:-02_filtering/results/porechop}"
THREADS="${3:-16}"

mkdir -p "${OUTDIR}"

TRIMMED_FASTQ="${OUTDIR}/trimmed_reads.fastq.gz"
LOG_FILE="${OUTDIR}/porechop_run.log"

echo "============================================================"
echo "[Porechop] Starting ONT Adapter Trimming & Chimera Splitting"
echo "Input FASTQ:   ${INPUT_FASTQ}"
echo "Output FASTQ:  ${TRIMMED_FASTQ}"
echo "Threads:       ${THREADS}"
echo "Log File:      ${LOG_FILE}"
echo "============================================================"

# Run Porechop with output logging
porechop \
  --input "${INPUT_FASTQ}" \
  --output "${TRIMMED_FASTQ}" \
  --threads "${THREADS}" \
  --verbosity 1 2>&1 | tee "${LOG_FILE}"

echo "============================================================"
echo "[Porechop] Adapter trimming complete!"
echo "Trimmed FASTQ: ${TRIMMED_FASTQ}"
echo "Run Summary:   ${LOG_FILE}"
echo "============================================================"
