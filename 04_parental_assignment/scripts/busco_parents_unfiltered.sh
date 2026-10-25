#!/usr/bin/env bash
# Script: 04_parental_assignment/scripts/busco_parents_unfiltered.sh
# Purpose: Concatenate unfiltered binned parental contigs and evaluate genome completeness using compleasm.

set -euo pipefail

source config.sh

parent_results_dir="04_parental_assignment/results"
raw_hap1_dir="${parent_results_dir}/binned_hap1"
raw_hap2_dir="${parent_results_dir}/binned_hap2"

decipiens_unfiltered="${parent_results_dir}/decipiens_unfiltered.fa"
virginea_unfiltered="${parent_results_dir}/virginea_unfiltered.fa"

mkdir -p "${parent_results_dir}"

# Concatenate unfiltered bins from both haplotypes
cat "${raw_hap1_dir}/decipiens.fa" \
    "${raw_hap2_dir}/decipiens.fa" \
    > "${decipiens_unfiltered}"

cat "${raw_hap1_dir}/virginea.fa" \
    "${raw_hap2_dir}/virginea.fa" \
    > "${virginea_unfiltered}"

# Define assemblies for compleasm run
ASSEMBLIES=(
    "decipiens_unfiltered.fa"
    "virginea_unfiltered.fa"
)

out_dir="${parent_results_dir}/compleasm_parents_unfiltered"
mkdir -p "${out_dir}"

lineage="eudicotyledons_odb12"

# Run compleasm on both unfiltered assemblies
for FASTA in "${ASSEMBLIES[@]}"; do
  BASE=$(echo "$FASTA" | sed 's/_unfiltered.fa//')

  echo "-------------------------------------------------------"
  echo "Processing Unfiltered: $BASE"
  echo "-------------------------------------------------------"

  compleasm run \
    -a "${parent_results_dir}/${FASTA}" \
    -o "${out_dir}/${BASE}" \
    -l "${lineage}" \
    -t 64 \
    2>&1 | tee "${out_dir}/compleasm_${BASE}_unfiltered.log"

  echo "Finished processing unfiltered $BASE"
done
