#!/usr/bin/env bash
# Script: 03_assembly/scripts/gfastats.sh
# Purpose: Get some key graphs from the hifiasm assemble by converting to FASTA files
set -euo pipefail

source config.sh

assembly_dir="03_assembly/results/hifiasm"
stats_dir="${assembly_dir}/post_hifi_stats"

mkdir -p "${stats_dir}"

#!/usr/bin/env bash
# Script: 03_assembly/scripts/gfastats.sh
# Purpose: Generate assembly stats via gfastats and convert GFAs to FASTA for downstream assessment

set -euo pipefail

source config.sh

assembly_dir="03_assembly/results/hifiasm"
stats_dir="${assembly_dir}/post_hifi_stats"

mkdir -p "${stats_dir}"

# Primary
gfastats "${assembly_dir}/MM_assembly.primary.bp.p_ctg.gfa" -t 16 --discover-paths > "${stats_dir}/stats_primary.txt"
# Convert to FASTA without raising --out-fasta error!
awk '/^S/{print ">"$2"\n"$3}' "${assembly_dir}/MM_assembly.primary.bp.p_ctg.gfa" > "${assembly_dir}/MM_assembly.primary.p_ctg.fa"

# Haplotype 1
gfastats "${assembly_dir}/MM_assembly.bp.hap1.p_ctg.gfa" -t 16 --discover-paths > "${stats_dir}/stats_hap1.txt"
awk '/^S/{print ">"$2"\n"$3}' "${assembly_dir}/MM_assembly.bp.hap1.p_ctg.gfa" > "${assembly_dir}/MM_assembly.hap1.p_ctg.fa"

# 3. Haplotype 2 Assembly
gfastats "${assembly_dir}/MM_assembly.bp.hap2.p_ctg.gfa" -t 16 --discover-paths > "${stats_dir}/stats_hap2.txt"
awk '/^S/{print ">"$2"\n"$3}' "${assembly_dir}/MM_assembly.bp.hap2.p_ctg.gfa" > "${assembly_dir}/MM_assembly.hap2.p_ctg.fa"
