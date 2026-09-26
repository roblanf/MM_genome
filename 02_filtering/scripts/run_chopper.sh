#!/usr/bin/env bash
# Script: 02_filtering/scripts/run_chopper.sh
# Purpose: Decompress files with pigz, filter reads using >=Q15 and length >=15kb using chopper.
# These cutoffs are set to provide good quality whilst maintaining >30x coverage per haplotype.

set -euo pipefail

source config.sh

filter_dir="02_filtering"
mkdir -p "${filter_dir}"

out_fastq="${filter_dir}/filtered_reads.fastq.gz"

# Use piping to decompress with pigz, then filter with chopper, and compress again
find ${raw_data} -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | \
  xargs pigz -dc -p 64 | \
  chopper -q 15 -l 15000 | \
  pigz -p 64 > "${out_fastq}"
