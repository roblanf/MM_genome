#!/usr/bin/env bash
# Script: 04_parental_assignment/scripts/align_parents.sh
# Purpose: Align haplotype assemblies (hap1 and hap2) against parental genomes using minimap2

set -euo pipefail

# Define directories and file paths
pa_dir="04_parental_assignment"
genomes_dir="${pa_dir}/parental_genomes"

decipiens_fa="${genomes_dir}/E_decipiens.fa"
virginea_fa="${genomes_dir}/E_virginea.fa"

hap1_fa="03_assembly/results/hifiasm/MM_assembly.hap1.p_ctg.fa"
hap2_fa="03_assembly/results/hifiasm/MM_assembly.hap2.p_ctg.fa"

# Ensure output directory exists (will if you already ran get_parents.sh)
mkdir -p "${pa_dir}"

# Aligning Haplotype Assemblies Against Parent Genomes
# Like for hap1 vs hap2 alignmnet, asm20 preset used due to suitability for up to several % divergence
# N=1000 means wide search space for secondary alignments but only best primary alignment outputted using --secondary=no.

# Hap1 aligned to decipiens
minimap2 -x asm20 -N 1000 --secondary=no -t 64 "${decipiens_fa}" "${hap1_fa}" > "${pa_dir}/results/hap1_vs_decipiens.paf"

# Hap1 aligned to virginea
minimap2 -x asm20 -N 1000 --secondary=no -t 64 "${virginea_fa}" "${hap1_fa}" > "${pa_dir}/results/hap1_vs_virginea.paf"

# Hap2 aligned to decipiens
minimap2 -x asm20 -N 1000 --secondary=no -t 64 "${decipiens_fa}" "${hap2_fa}" > "${pa_dir}/results/hap2_vs_decipiens.paf"

# Hap2 aligned to virginea
minimap2 -x asm20 -N 1000 --secondary=no -t 64 "${virginea_fa}" "${hap2_fa}" > "${pa_dir}/results/hap2_vs_virginea.paf"

