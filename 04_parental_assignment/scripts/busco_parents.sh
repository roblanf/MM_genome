#!/usr/bin/env bash
# Script: 03_assembly/scripts/busco_parents.sh
# Purpose: Concatenate binned parental contigs and evaluate genome completeness using compleasm.

set -euo pipefail

source config.sh

# 1. Prepare directory and unified parental FASTA files
parent_results_dir="04_parental_assignment/results"
mkdir -p "${parent_results_dir}"

decipiens_fa="${parent_results_dir}/decipiens_assembly.fa"
virginea_fa="${parent_results_dir}/virginea_assembly.fa"

echo "-------------------------------------------------------"
echo "Concatenating binned haplotigs into unified parental FASTAs..."
echo "-------------------------------------------------------"

cat "${parent_results_dir}/binned_hap1_filtered/decipiens.fa" \
    "${parent_results_dir}/binned_hap2_filtered/decipiens.fa" \
    > "${decipiens_fa}"

cat "${parent_results_dir}/binned_hap1_filtered/virginea.fa" \
    "${parent_results_dir}/binned_hap2_filtered/virginea.fa" \
    > "${virginea_fa}"

# 2. Define assemblies for compleasm run
declare -A ASSEMBLIES=(
    ["decipiens"]="${decipiens_fa}"
    ["virginea"]="${virginea_fa}"
)

# 3. Output directory for completeness results
compleasm_dir="03_assembly/results/compleasm_parents"
mkdir -p "${compleasm_dir}"

lineage="eudicotyledons_odb12"

# 4. Run compleasm on both consolidated parental sets
for BASE in "${!ASSEMBLIES[@]}"; do
    FASTA_PATH="${ASSEMBLIES[$BASE]}"

    echo "-------------------------------------------------------"
    echo "Processing parental assembly: $BASE"
    echo "-------------------------------------------------------"

    compleasm run \
        -a "${FASTA_PATH}" \
        -o "${compleasm_dir}/${BASE}" \
        -l "${lineage}" \
        -t 128 \
        2>&1 | tee "${compleasm_dir}/compleasm_${BASE}.log"

    echo "Finished processing $BASE"
done
