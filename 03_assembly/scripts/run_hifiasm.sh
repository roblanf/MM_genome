#!/usr/bin/env bash
# Script: 03_assembly/scripts/run_hifiasm.sh
# Purpose: Primary genome assembly of filtered ONT reads using hifiasm.

set -euo pipefail

source config.sh

# Making output directory for results
assembly_dir="03_assembly/results/hifiasm"
mkdir -p "${assembly_dir}"

# Run hifiasm on filtered reads

hifiasm \
  -o "${assembly_dir}/MM_assembly" \
  -t 100 \
  --ont \
  -l 3 \
  --telo-m AAACCCT \
  --dual-scaf \
  "${filtered_fastq}" \
  2>&1 | tee "${assembly_dir}/hifiasm.log"

# Note for those testing prior to running on full data:
# Please run the following line to remove outputs from the test run, BEFORE you run this script again on the full code
# rm -f 03_assembly/results/hifiasm/MM_assembly.*
