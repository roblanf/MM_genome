#!/usr/bin/env bash
# Script: 03_assembly/scripts/gfastats.sh
# Purpose: Get some key graphs from the hifiasm assemble by converting to FASTA files
set -euo pipefail

source config.sh

assembly_dir="03_assembly/results/hifiasm"
stats_dir="${assembly_dir}/post_hifi_stats"

mkdir -p "${stats_dir}"

# Compute gfastats whilst also converting to .fa for use in compleasm_busco.sh

# Primary
gfastats "${assembly_dir}/MM_assembly.primary.bp.p_ctg.gfa" -t 16 --discover-paths --out-fasta "${assembly_dir}/MM_assembly.primary.p_ctg.fa" > "${stats_dir}/stats_primary.txt"

# Hap1
gfastats "${assembly_dir}/MM_assembly.bp.hap1.p_ctg.gfa" -t 16 --discover-paths --out-fasta "${assembly_dir}/MM_assembly.hap1.p_ctg.fa" > "${stats_dir}/stats_hap1.txt"

# Hap2
gfastats "${assembly_dir}/MM_assembly.bp.hap2.p_ctg.gfa" -t 16 --discover-paths --out-fasta "${assembly_dir}/MM_assembly.hap2.p_ctg.fa" > "${stats_dir}/stats_hap2.txt"
