#!/usr/bin/env bash
# Script: 01_qc/scripts/run_gc_check.sh
# Purpose: Calculate global statistics and per-read GC content distribution for ONT reads

set -euo pipefail
source config.sh

gc_dir="01_qc/results/gc"
mkdir -p "${gc_dir}"

# Global yield, coverage, and mean GC% stats
seqkit stats \
  -j 64 \
  -a \
  ${raw_data}/*.fastq.gz \
  > "${gc_dir}/read_stats.txt"

# Extract per-read length and GC content - used to make histogram in R afterwards
seqkit fx2tab \
  --threads 64 \
  --name \
  --length \
  --gc \
  ${raw_data}/*.fastq.gz \
  > "${gc_dir}/gc_pcts.tsv"
