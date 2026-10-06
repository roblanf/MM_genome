#!/usr/bin/env bash
# Script: 03_assembly/scripts/run_compleasm.sh
# Purpose: Evaluating genome completeness with BUSCO for hifiasm assemblies using compleasm and eudicotyledons_odb12.

set -euo pipefail

source config.sh

# Define assemblies to be used for looping with compleasm
ASSEMBLIES=(
  "MM_assembly.primary.p_ctg.fa"
  "MM_assembly.hap1.p_ctg.fa"
  "MM_assembly.hap2.p_ctg.fa"
)

# Making output directory for results
compleasm_dir="03_assembly/results/compleasm"
mkdir -p "${compleasm_dir}"

lineage="eudicotyledons_odb12"

# Loop through each assembly file defined in ASSEMBLIES array
for FASTA in "${ASSEMBLIES[@]}"; do
  # Extract label ("primary", "hap1", or "hap2") by removing other replacing other parts of file names with nothing.
  BASE=$(echo "$FASTA" | sed 's/.*MM_assembly.//; s/.p_ctg.fa//; s/.fa//')

    echo "-------------------------------------------------------"
    echo "Processing: $BASE"
    echo "-------------------------------------------------------"
# Now iterate compleasm through the primary, hap1 and hap2 fasta files
  compleasm run \
    -a "03_assembly/results/hifiasm/${FASTA}" \
    -o "${compleasm_dir}/${BASE}" \
    -l "${lineage}" \
    -t 64 \
    2>&1 | tee "${compleasm_dir}/compleasm_${BASE}.log"

  echo "Finished processing $BASE"
done
