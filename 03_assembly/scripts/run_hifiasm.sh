#!/usr/bin/env bash
# Script: 03_assembly/scripts/run_hifiasm.sh
# Purpose: Primary genome assembly of filtered ONT reads using hifiasm.

set -euo pipefail

source config.sh

filter_dir="02_filtering"
assembly_dir="03_assembly/results/hifiasm"
filtered_fastq="${filter_dir}/filtered_reads.fastq.gz"
ramdisk_dir="/mnt/ramdisk"

mkdir -p "${assembly_dir}"

#Set up the RAM disk with sudo:
sudo mkdir -p "${ramdisk_dir}"
sudo mount -t tmpfs -o size=1500G tmpfs "${ramdisk_dir}"
sudo chown "${USER}" "${ramdisk_dir}"

# Copy filtered data into the RAM disk
cp "${filtered_fastq}" "${ramdisk_dir}/"

# Run hifiasm on filtered reads in RAM disk
cd "${ramdisk_dir}"

hifiasm \
  -o MM_assembly \
  -t 100 \
  --ont \
  -l 3 \
  --telo-m AAACCCT \
  --dual-scaf \
  filtered_reads.fastq.gz \
  2>&1 | tee hifiasm.log

# Move results back into project directory:
# returns you to previous directory before copying results.
cd - > /dev/null
rsync -av --exclude='*.fastq.gz' "${ramdisk_dir}/" "${assembly_dir}/"

# Unmount RAM disk
sudo umount "${ramdisk_dir}" || true
