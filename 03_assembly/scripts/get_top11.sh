#!/usr/bin/env bash
# Script: 03_assembly/scripts/get_top11.sh
# Purpose: Calculate top 11 contig metrics - useful because there are 11 base chromosomes, 2n = 22.

set -euo pipefail

source config.sh

stats_dir="03_assembly/results/hifiasm/post_hifi_stats"

# Ensure output directory exists
mkdir -p "${stats_dir}"

# --- Primary Assembly ---
tot_primary=$(grep "Total scaffold length:" "${stats_dir}/stats_primary.txt" | awk '{print $NF}')
top11_primary=$(tail -n +2 "${stats_dir}/stats_primary_segments.txt" | sort -k4,4nr | head -n 11 | cut -f4 | paste -sd+ | bc)
pct_primary=$(bc <<< "scale=2; (${top11_primary} * 100) / ${tot_primary}")

# --- Haplotype 1 ---
tot_hap1=$(grep "Total scaffold length:" "${stats_dir}/stats_hap1.txt" | awk '{print $NF}')
top11_hap1=$(tail -n +2 "${stats_dir}/stats_hap1_segments.txt" | sort -k4,4nr | head -n 11 | cut -f4 | paste -sd+ | bc)
pct_hap1=$(bc <<< "scale=2; (${top11_hap1} * 100) / ${tot_hap1}")

# --- Haplotype 2 ---
tot_hap2=$(grep "Total scaffold length:" "${stats_dir}/stats_hap2.txt" | awk '{print $NF}')
top11_hap2=$(tail -n +2 "${stats_dir}/stats_hap2_segments.txt" | sort -k4,4nr | head -n 11 | cut -f4 | paste -sd+ | bc)
pct_hap2=$(bc <<< "scale=2; (${top11_hap2} * 100) / ${tot_hap2}")

# --- Save to Results File & Print to Screen ---
out_file="${stats_dir}/top11_metrics.txt"

{
  echo "PRIMARY: Total=${tot_primary} bp | Top11=${top11_primary} bp | Pct=${pct_primary}%"
  echo "HAP1:    Total=${tot_hap1} bp | Top11=${top11_hap1} bp | Pct=${pct_hap1}%"
  echo "HAP2:    Total=${tot_hap2} bp | Top11=${top11_hap2} bp | Pct=${pct_hap2}%"
} | tee "${out_file}"
