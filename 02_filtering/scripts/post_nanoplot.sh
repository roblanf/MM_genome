#!/usr/bin/env bash
# Script: 02_filtering/scripts/post_nanoplot.sh
# Purpose: Runs NanoPlot on filtered reads for comparison after QC and filtering.

set -euo pipefail
source config.sh

# Make directory in results for NanoPlot output
filter_dir="02_filtering"
nanoplot_dir="${filter_dir}/results/nanoplot"

mkdir -p "${nanoplot_dir}"

# Run NanoPlot again, now on chopper-filtered output
NanoPlot \
  --fastq "${filter_dir}/filtered_reads.fastq.gz" \
  -t 64 \
  --downsample 100000 \
  -o "${nanoplot_dir}" \
  --title "E. phylacis filtered reads NanoPlot (>Q15, >20Kb length)"
