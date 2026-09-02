#!/usr/bin/env bash
# Long-read k-mer counting with KMC and GenomeScope 2.0 histogram
# Usage: bash 01_qc/scripts/run_kmc.sh [raw_data_dir] [qc_out_dir]

set -euo pipefail

# Set defaults to match current directory structure
RAW_DATA_DIR="${1:-raw_data}"
QC_DIR="${2:-01_qc/results}"

# Output directory for k-mer results - goes into 01_qc/results
KMER_DIR="${QC_DIR}/kmer_distribution"
TMP_DIR="tmp_kmc_processing" #will be deleted at the end

echo "=== Starting KMC K-mer Analysis ==="
mkdir -p "${KMER_DIR}"
mkdir -p "${TMP_DIR}"

# Step 1: Find all FASTQ files (compressed or uncompressed)
find "${RAW_DATA_DIR}" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) > "${TMP_DIR}/files.txt"

# Step 2: Run KMC
# Adjust threads (-t) and memory (-m) safely for test runs
kmc -k21 -t4 -m16 -ci1 -cs10000 \
    @"${TMP_DIR}/files.txt" \
    "${TMP_DIR}/kmc_db" \
    "${TMP_DIR}/"

# Step 3: Create histogram for GenomeScope2
kmc_tools transform "${TMP_DIR}/kmc_db" histogram "${KMER_DIR}/long_read_histogram.txt" -ci1 -cs1000000

# Step 4: Run GenomeScope 2.0 (diploid)
genomescope2 \
    -i "${KMER_DIR}/long_read_histogram.txt" \
    -o "${KMER_DIR}/genomescope_results" \
    -k 21 \
    -p 2

# Clean up temporary database files
rm -rf "${TMP_DIR}"

echo "=== KMC & GenomeScope Analysis Complete ==="
