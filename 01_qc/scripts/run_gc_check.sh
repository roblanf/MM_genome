#!/usr/bin/env bash
# Quick GC distribution check using SeqKit
# Usage: bash 01_qc/scripts/run_gc_check.sh <path_to_fastq_file> [num_reads]

set -euo pipefail

# Default to 100,000 reads if 2nd argument isn't provided
FASTQ_INPUT="${1:-raw_data/tiny_test.fastq}"
NUM_READS="${2:-100000}"

echo "Analyzing GC content for top ${NUM_READS} reads in: ${FASTQ_INPUT}"
echo "------------------------------------------------------------"

seqkit head -n "${NUM_READS}" "${FASTQ_INPUT}" | \
  seqkit fx2tab -n -g | \
  cut -f 2 | \
  awk '{printf "%.0f\n", $1}' | \
  sort | \
  uniq -c | \
  awk '{printf "%s%%\t%s\n", $2, $1}' | \
  sort -n
