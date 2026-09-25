#!/usr/bin/env bash
# Script: 01_qc/scripts/run_nanoplot.sh
# Purpose: Runs NanoPlot on raw ONT FASTQ reads (pre-filtering baseline).

set -euo pipefail

# loads config file to use raw or test data depending on user choice
source config.sh

nanoplot_dir="01_qc/results/nanoplot"

mkdir -p "${nanoplot_dir}"

# Run NanoPlot
NanoPlot -t 128 \
         --fastq "${raw_data}"/*.fastq.gz \
         --downsample 100000 \
         -o "${nanoplot_dir}" \
         --title "E. phylacis ONT reads NanoPlot"
#Uses 128 threads, change at your discretion depending on server capacity
# downsample 100000 means randomly sample max of 100K reads to generate nanoplot and summary stats, avoid computation overload.
