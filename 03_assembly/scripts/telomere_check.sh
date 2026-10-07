#!/usr/bin/env bash
# Script: 03_assembly/scripts/telomere_check.sh
# Purpose: Identify and profile telomere repeat motif(s) in the hifiasm assemblies using tidk.

set -euo pipefail

source config.sh

# Target motifs confirmed from tidk explore (canonical plant telomere AAACCCT and others, maybe centromeres)
MOTIFS=("AAACCCT" "TTTAGGG" "ACCCGTC" "AAAAAAT" "AAAAAAG")
ASSEMBLIES=(
  "MM_assembly.primary.p_ctg.fa"
  "MM_assembly.hap1.p_ctg.fa"
  "MM_assembly.hap2.p_ctg.fa"
)

tidk_dir="03_assembly/results/tidk"
mkdir -p "${tidk_dir}"

# For each assembly, tidk searches for the motifs, pipes it to the directory and output is made as a tsv file to be viewed.

for fa in "${ASSEMBLIES[@]}"; do
  for m in "${MOTIFS[@]}"; do
    tidk search --string "$m" \
      --dir "${tidk_dir}" \
      --output "${fa%.fa}_${m}" \
      --extension tsv \
      "03_assembly/results/hifiasm/${fa}"
  done
done
