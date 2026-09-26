#!/usr/bin/env bash
# Script: 03_assembly/scripts/run_hifiasm.sh
# Purpose: Primary genome assembly of filtered ONT reads using hifiasm.

set -euo pipefail

source config.sh

# Making output directory for results
assembly_dir="03_assembly/results/hifiasm"
mkdir -p "${assembly_dir}"

# Run hifiasm on filtered reads
cd "${ramdisk_dir}"

# in config.sh the data is named filtered_fastq whether you set to test or full filtered reads
# so this will make the script work for both!
input_file=$(basename "${filtered_fastq}")

hifiasm \
  -o MM_assembly \
  -t 100 \
  --ont \
  -l 3 \
  --telo-m AAACCCT \
  --dual-scaf \
  "${input_file}" \
  2>&1 | tee hifiasm.log
