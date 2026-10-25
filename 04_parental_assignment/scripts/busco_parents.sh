#!/usr/bin/env bash
# Script: 03_assembly/scripts/busco_parents.sh
# Purpose: Concatenate binned parental contigs and evaluate genome completeness using compleasm.

set -euo pipefail

source config.sh

# Prepare directory and unified parental FASTA files
parent_results_dir="04_parental_assignment/results"
mkdir -p "${parent_results_dir}"

decipiens_fa="${parent_results_dir}/decipiens_assembly.fa"
virginea_fa="${parent_results_dir}/virginea_assembly.fa"

cat "${parent_results_dir}/binned_hap1_filtered/decipiens.fa" \
    "${parent_results_dir}/binned_hap2_filtered/decipiens.fa" \
    > "${decipiens_fa}"

cat "${parent_results_dir}/binned_hap1_filtered/virginea.fa" \
    "${parent_results_dir}/binned_hap2_filtered/virginea.fa" \
    > "${virginea_fa}"

# Define assemblies for compleasm run
ASSEMBLIES=(
    "decipiens_all.fa"
    "virginea_all.fa"
)

out_dir="04_parental_assignment/results/compleasm_parents"
mkdir -p "${out_dir}"

lineage="eudicotyledons_odb12"

# Run compleasm on both parental haplotype assemblies
for FASTA in "${ASSEMBLIES[@]}"; do
  # Extract label ("decipiens" or "virginea") by removing other replacing other parts of file names with nothing.
  BASE=$(echo "$FASTA" | sed 's/_all.fa//; s/.fa//')

    echo "-------------------------------------------------------"
    echo "Processing: $BASE"
    echo "-------------------------------------------------------"

# Now iterate compleasm through the decipiens, virginea files
  compleasm run \
    -a "${parent_results_dir}/${FASTA}" \
    -o "${out_dir}/${BASE}" \
    -l "${lineage}" \
    -t 64 \
    2>&1 | tee "${out_dir}/compleasm_${BASE}.log"

  echo "Finished processing $BASE"
done

