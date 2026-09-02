#!/usr/bin/env bash
# Script: 02_filtering/scripts/run_chopper.sh
# Purpose: Filter ONT reads by length and mean quality score using Chopper.
# Usage: bash 02_filtering/scripts/run_chopper.sh [input_fastq] [output_dir] [min_length] [min_qual] [threads]

set -euo pipefail

# Inputs and defaults (defaults to output of Porechop)
INPUT_FASTQ="${1:-02_filtering/results/porechop/trimmed_reads.fastq.gz}"
OUTDIR="${2:-02_filtering/results/chopper}"
MIN_LENGTH="${3:-1000}"   # Filter out reads < 1kb
MIN_QUAL="${4:-10}"       # Filter out reads with mean Q-score < 10 (90% accuracy)
THREADS="${5:-16}"

mkdir -p "${OUTDIR}"

FILTERED_FASTQ="${OUTDIR}/filtered_reads.fastq.gz"
LOG_FILE="${OUTDIR}/chopper_run.log"

echo "============================================================"
echo "[Chopper] Starting Read Filtering"
echo "Input FASTQ:   ${INPUT_FASTQ}"
echo "Output FASTQ:  ${FILTERED_FASTQ}"
echo "Min Length:    >= ${MIN_LENGTH} bp"
echo "Min Quality:   >= Q${MIN_QUAL}"
echo "Threads:       ${THREADS}"
echo "============================================================"

# Stream gzipped FASTQ into Chopper and compress output directly
gunzip -c "${INPUT_FASTQ}" | \
  chopper -l "${MIN_LENGTH}" -q "${MIN_QUAL}" --threads "${THREADS}" | \
  gzip -c > "${FILTERED_FASTQ}"

echo "============================================================"
echo "[Chopper] Filtering complete!"
echo "Filtered FASTQ: ${FILTERED_FASTQ}"
echo "============================================================"
