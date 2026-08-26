#!/bin/bash
# Ensure test_data directory exists
mkdir -p test_data

# Grab the first 10,000 reads (40,000 lines) from the raw FASTQ
zcat /data/raw_data/MM/2025_long_reads/ONT_gDNA192_1071_RL/E_Phylacis/20251117_1227_2C_PBE85256_e040940c/fastq_pass/PBE85256_pass_e040940c_53a5bbe6_19.fastq.gz | head -n 40000 > test_data/tiny_test.fastq
echo "Test dataset successfully generated at test_data/tiny_test.fastq"
