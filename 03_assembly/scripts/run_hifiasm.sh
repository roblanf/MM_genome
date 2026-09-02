#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Script: run_hifiasm.sh
# Purpose: Run Hifiasm assembly in RAM disk using filtered raw_reads.fastq.gz
# Location: MM_genome/03_assembly/scripts/run_hifiasm.sh
# ==============================================================================

# 1. Setup Directories & Relative Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

DATA_DIR="${PROJECT_ROOT}/02_filtering"
OUT_DIR="${PROJECT_ROOT}/03_assembly"
LOG_DIR="${OUT_DIR}/logs"
RAMDISK="/mnt/ramdisk"

# Input reads file name from 02_filtering
FILTERED_READS="${DATA_DIR}/raw_reads.fastq.gz"
SAMPLE_NAME="MM_hifiasm_asm"
THREADS=160

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

echo "=== [1/5] Setting up RAM Disk ==="
sudo mkdir -p "${RAMDISK}"
sudo mount -t tmpfs -o size=1500G tmpfs "${RAMDISK}" || true
sudo chown "$USER" "${RAMDISK}"

echo "=== [2/5] Copying Filtered Reads to RAM Disk ==="
if [ ! -f "${FILTERED_READS}" ]; then
    echo "Error: Filtered read file not found at ${FILTERED_READS}"
    exit 1
fi

cp "${FILTERED_READS}" "${RAMDISK}/"
READS_BASENAME="$(basename "${FILTERED_READS}")"

echo "=== [3/5] Running Hifiasm Assembly ==="
cd "${RAMDISK}"

# Run Hifiasm 
# (Note: Keep '--ont' if using Oxford Nanopore ultra-long reads, remove if using PacBio HiFi)
hifiasm \
    -o "${SAMPLE_NAME}" \
    -t "${THREADS}" \
    --ont \
    -l 3 \
    --telo-m AAACCCT \
    --dual-scaf \
    "${READS_BASENAME}" \
    2>&1 | tee "${LOG_DIR}/hifiasm_run.log"

echo "=== [4/5] Syncing Assembly Results back to ${OUT_DIR} ==="
# Sync results back while excluding the raw input copy from RAM disk
rsync -av --progress --exclude='*.fastq.gz' "${RAMDISK}/" "${OUT_DIR}/"

echo "=== [5/5] Cleaning up RAM Disk ==="
# Remove temporary files and unmount RAM disk to release system memory
rm -rf "${RAMDISK:?}"/*
cd "${PROJECT_ROOT}"
sudo umount "${RAMDISK}" || true

echo "=== Assembly Complete ==="
echo "Outputs successfully saved to: ${OUT_DIR}"
