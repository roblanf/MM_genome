#!/bin/bash
set -euo pipefail
# Make sure test_data directory exists
mkdir -p test_data
# Extract the first 10,000 reads (40,000 lines) into test_data/
zcat /data/raw_data/MM/2025_long_reads/ONT_gDNA192_1071_RL/E_Phylacis/20251117_1227_2C_PBE85256_e040940c/fastq_fail/PBE85256_fail_e040940c_53a5bbe6_41.fastq.gz | head -n 40000 > test_data/tiny_test.fastq
echo "Test dataset successfully generated at test_data/tiny_test.fastq"

bash make_test_data.sh
