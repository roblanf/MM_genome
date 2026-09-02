#!/usr/bin/env bash
# Script: 01_qc/scripts/run_kmc.sh
# Purpose: Run KMC k-mer counting and GenomeScope 2.0 on long-read FASTQ data.
# Usage: bash 01_qc/scripts/run_kmc.sh [input_fastq] [out_dir] [kmer_len] [threads] [max_mem_gb]

set -euo pipefail

# Parameters & Defaults
INPUT_FASTQ="${1:-raw_data/raw_reads.fastq.gz}"
OUTDIR="${2:-01_qc/results/kmc_genomescope}"
KMER_LEN="${3:-21}"          # k=21 is standard for genome size estimation
THREADS="${4:-16}"
MAX_MEM="${5:-64}"           # Max RAM memory limit in GB for KMC

# Directories & Temporary Files
mkdir -p "${OUTDIR}"
TMP_DIR="${OUTDIR}/tmp"
mkdir -p "${TMP_DIR}"

KMC_DB="${OUTDIR}/kmc_db"
HIST_FILE="${OUTDIR}/kmc_k${KMER_LEN}.hist"
GS_OUTDIR="${OUTDIR}/genomescope_k${KMER_LEN}"

echo "============================================================"
echo "[KMC & GenomeScope2] Starting k-mer analysis"
echo "Input FASTQ:   ${INPUT_FASTQ}"
echo "Output Dir:    ${OUTDIR}"
echo "K-mer Length:  k=${KMER_LEN}"
echo "Threads:       ${THREADS}"
echo "Max Memory:    ${MAX_MEM} GB"
echo "============================================================"

# Step 1: Count k-mers with KMC
echo "[1/3] Running KMC k-mer counting..."
# -fm: FASTQ mode; -ci1: count min threshold 1; -cs10000: count max threshold
kmc \
  -k"${KMER_LEN}" \
  -t"${THREADS}" \
  -m"${MAX_MEM}" \
  -ci1 \
  -cs100000 \
  -fq \
  "${INPUT_FASTQ}" \
  "${KMC_DB}" \
  "${TMP_DIR}"

# Step 2: Generate k-mer frequency histogram
echo "[2/3] Generating k-mer frequency histogram..."
kmc_tools transform "${KMC_DB}" histogram "${HIST_FILE}" -cx10000

# Step 3: Run GenomeScope 2.0
echo "[3/3] Running GenomeScope 2.0..."
genomescope2 \
  -i "${HIST_FILE}" \
  -o "${GS_OUTDIR}" \
  -k "${KMER_LEN}" \
  -p 2 \
  -l "${READ_TYPE:-long}" \
  --verbose

# Clean up temporary KMC working files
rm -rf "${TMP_DIR}"

echo "============================================================"
echo "[KMC & GenomeScope2] Pipeline complete!"
echo "Histogram saved to:         ${HIST_FILE}"
echo "GenomeScope results saved:  ${GS_OUTDIR}"
echo "============================================================"
