#!/usr/bin/env bash
# Script: 03_assembly/scripts/test_data_asm.sh
# Purpose: Makes small subset (40K) of the filtered reads to be used for testing hifiasm script/assembly.
#          In config.sh, see the For Hifiasm section and toggle between this test data and real data.
#          You have to use this instead of test_data because the Hifiasm script needs filtered_fastq not raw_data to run.

zcat 02_filtering/filtered_reads.fastq.gz | head -n 40000 | gzip > 02_filtering/test_filtered.fastq.gz
