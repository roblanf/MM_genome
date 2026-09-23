#!/usr/bin/env bash
# Script: 01_qc/scripts/run_nanoplot.sh
# Purpose: Runs NanoPlot on raw ONT FASTQ reads (pre-filtering baseline).
# Usage: bash 01_qc/scripts/run_nanoplot.sh [input_fastq] [output_dir] [threads]

#!/usr/bin/env bash
set -euo pipefail

INPUT_FASTQ="${1:-raw_data}"
OUTDIR="${2:-01_qc/results/nanoplot_test_output}"
THREADS="${3:-16}"

mkdir -p "${OUTDIR}"

NanoPlot \
  --threads "${THREADS}" \
  --fastq "${INPUT_FASTQ}"/*fastq* \
  --outdir "${OUTDIR}" \
  --loglength \
  --N50 \
  --title "E. phylacis Raw ONT QC" \
  --tsv_stats
