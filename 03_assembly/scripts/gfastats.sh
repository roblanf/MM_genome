#!/usr/bin/env bash
# Script: 03_assembly/scripts/gfastats.sh
# Purpose: Get some key graphs from the hifiasm assemble by converting to FASTA files
set -euo pipefail

source config.sh

assembly_dir="03_assembly/results/hifiasm"
stats_dir="${assembly_dir}/post_hifi_stats"

mkdir -p "${stats_dir}"

gfastats "${assembly_dir}/MM_assembly.bp.hap1.p_ctg.gfa" --discover-paths --segment-report > "${stats_dir}/stats_hap1_segments.txt"
gfastats "${assembly_dir}/MM_assembly.bp.hap1.p_ctg.gfa" --discover-paths > "${stats_dir}/stats_hap1.txt"

gfastats "${assembly_dir}/MM_assembly.bp.hap2.p_ctg.gfa" --discover-paths --segment-report > "${stats_dir}/stats_hap2_segments.txt"
gfastats "${assembly_dir}/MM_assembly.bp.hap2.p_ctg.gfa" --discover-paths > "${stats_dir}/stats_hap2.txt"

gfastats "${assembly_dir}/MM_assembly.bp.p_ctg.gfa" --discover-paths --segment-report > "${stats_dir}/stats_primary_segments.txt"
gfastats "${assembly_dir}/MM_assembly.bp.p_ctg.gfa" --discover-paths > "${stats_dir}/stats_primary.txt"
