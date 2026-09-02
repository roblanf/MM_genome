#!/usr/bin/env bash
# Script: 01_qc/scripts/run_nanoplot.sh
# Purpose: Runs NanoPlot on raw ONT FASTQ reads (pre-filtering baseline).
# Usage: bash 01_qc/scripts/run_nanoplot.sh [input_fastq] [output_dir] [threads]

set -euo pipefail

INPUT_FASTQ="${1:-raw_data/raw_reads.fastq.gz}"
OUTDIR="${2:-01_qc/results/raw_nanoplot}"
THREADS="${3:-16}"

mkdir -p "${OUTDIR}"

echo "============================================================"
echo "[NanoPlot] Running QC on raw reads"
echo "Input:   ${INPUT_FASTQ}"
echo "Output:  ${OUTDIR}"
echo "Threads: ${THREADS}"
echo "============================================================"

NanoPlot \
  --fastq "${INPUT_FASTQ}" \
  --outdir "${OUTDIR}" \
  --threads "${THREADS}" \
  --loglength \
  --N50 \
  --title "Raw ONT reads QC" \
  --tsv_stats

echo "============================================================"
echo "[NanoPlot] QC complete! HTML report saved in ${OUTDIR}"
echo "============================================================"
