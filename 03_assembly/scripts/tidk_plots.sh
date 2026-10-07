#!/usr/bin/env bash
# Script: 03_assembly/scripts/tidk_plots.sh
# Purpose: Get plots of repeat distributions for primary, h1, h2 assemblies on identified motifs

set -euo pipefail

source config.sh

mkdir -p 03_assembly/results/tidk/plots

# Finds tsv files made in last script for all assemblies, plots canonical motif distribution with tidk plot:

for asm in "primary.p_ctg" "hap1.p_ctg" "hap2.p_ctg"; do
  tidk plot \
    --tsv "03_assembly/results/tidk/MM_assembly.${asm}_AAACCCT_telomeric_repeat_windows.tsv" \
    --output "03_assembly/results/tidk/plots/${asm}_AAACCCT_plot"
done

