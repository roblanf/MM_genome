#!/usr/bin/env bash
# Script: 01_qc/scripts/run_gc_check.sh
# Purpose: Calculate GC content distribution and overall read statistics for raw ONT reads.
# Usage: bash 01_qc/scripts/run_gc_check.sh [input_fastq] [output_dir] [threads]

set -euo pipefail

INPUT_FASTQ="${1:-raw_data/raw_reads.fastq.gz}"
OUTDIR="${2:-01_qc/results/gc_content}"
THREADS="${3:-16}"

mkdir -p "${OUTDIR}"

STATS_FILE="${OUTDIR}/read_stats.txt"
GC_PER_READ_FILE="${OUTDIR}/gc_per_read.tsv"

echo "============================================================"
echo "[GC Check] Starting sequence analysis"
echo "Input FASTQ:   ${INPUT_FASTQ}"
echo "Output Dir:    ${OUTDIR}"
echo "Threads:       ${THREADS}"
echo "============================================================"

# Step 1: Calculate global summary stats (including average GC%)
echo "[1/2] Computing global yield and mean GC% with SeqKit..."
seqkit stats -j "${THREADS}" -a "${INPUT_FASTQ}" > "${STATS_FILE}"

# Step 2: Extract per-read length and GC content for distribution plotting
echo "[2/2] Extracting per-read length and GC content..."
seqkit fx2tab \
  --threads "${THREADS}" \
  --name \
  --length \
  --gc \
  "${INPUT_FASTQ}" > "${GC_PER_READ_FILE}"

echo "============================================================"
echo "[GC Check] Complete!"
echo "Summary stats:    ${STATS_FILE}"
echo "Per-read metrics: ${GC_PER_READ_FILE}"
echo "============================================================"
