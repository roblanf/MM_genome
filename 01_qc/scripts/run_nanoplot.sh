#!/usr/bin/env bash
# Runs NanoPlot on raw ONT fastq_pass reads (pre-filtering baseline).
# Usage: run_nanoplot.sh <fastq_pass_dir> <output_dir> [threads]
set -euo pipefail
fastq_dir=$1
outdir=$2
threads=${3:-16}

mkdir -p "${outdir}"
NanoPlot \
  --fastq ${fastq_dir}/*.fastq* \
  --outdir "${outdir}" \
  --threads "${threads}" \
  --loglength \
  --N50 \
  --title "Raw ONT reads QC" \
  --tsv_stats
