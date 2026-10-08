#!/usr/bin/env bash
# Script: 03_assembly/scripts/h1h2_alignment.sh
# Purpose: Align all contigs of Hap1 against Hap2 using minimap2.

set -euo pipefail

source config.sh

# 1. Output directory setup
out_dir="03_assembly/results/h1h2"
mkdir -p "${out_dir}"

hap1="03_assembly/results/hifiasm/MM_assembly.hap1.p_ctg.fa"
hap2="03_assembly/results/hifiasm/MM_assembly.hap2.p_ctg.fa"
paf_file="${out_dir}/hap1_vs_hap2.paf"

# 2. Align Hap1 (Query) to Hap2 (Target/Reference)
# asm20 preset used, suitable for up to several % divergence, good because parents are distant and we already saw ab 4.25% in GenomeScope
# Wide search space for secondary alignments but only best primary alignment outputted.
minimap2 -x asm20 -t 64 -N 1000 --secondary=no \
    "${hap2}" \
    "${hap1}" \
    > "${paf_file}"
