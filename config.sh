#!/bin/bash
set -euo pipefail
#!/bin/bash

# Configure between using test data or real data for running the scripts!
# To switch between test/real data, remove the # symbol from what you want to use and add # to the other line.
# Happy coding!

# FOR QC AND FILTERING:

# 1. TEST DATA
#raw_data="test_data"

# 2. REAL DATA
raw_data="/data/raw_data/MM/2025_long_reads/ONT_gDNA192_1071_RL/E_Phylacis/20251117_1227_2C_PBE85256_e040940c/fastq_pass"

# -------------------------------------------------------------------------------------------------------------------------------

# FOR HIFIASM:
# remember to take # off the one you want to use, and put # on the one you don't want to use.

# TEST DATA
#filtered_fastq="02_filtering/test_filtered_reads.fastq.gz"

# FULL DATA
filtered_fastq="02_filtering/filtered_reads.fastq.gz"
