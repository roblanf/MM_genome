#!/usr/bin/env bash
# Script: 01_qc/scripts/run_kmc.sh
# Purpose: Run KMC k-mer counting and GenomeScope 2.0 on long-read FASTQ data.

set -euo pipefail
source config.sh

kmc_dir="01_qc/results/kmc"
mkdir -p "${kmc_dir}"
tmp_dir="kmc_dir/tmp"
mkdir -p "${tmp_dir}"

find ${raw_data} -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) > "${tmp_dir}"/files.txt

# First, count kmers with KMC.
# Again please be careful to set threads according to server status.
# k=21 is standard for genome size estimation.
# including k-mers that occur at least once (ci1)
#don't track k-mers past 10000 occurences, saves RAM and disk space (cs10000)

kmc \
  -k21 \
  -t64 \
  -m256 \
  -ci1 \
  -cs10000 \
  @"${tmp_dir}/files.txt" \
   "${tmp_dir}/kmc_db" \
    "${tmp_dir}"/

# Generating k-mer frequency histogram
#separate histogram lines up to 10000 occurrences.
kmc_tools transform "${tmp_dir}/kmc_db" histogram "${kmc_dir}/lr_histogram.txt" -ci1 -cs1000000

# Making GenomeScope 2.0 profile
genomescope2 \
  -i "${kmc_dir}/lr_histogram.txt" \
  -o "${kmc_dir}" \
  -k 21 \
  -p 2 \
  --verbose

rm -rf "${tmp_dir}"
