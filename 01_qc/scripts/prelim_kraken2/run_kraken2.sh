#!/usr/bin/env bash
# Runs Kraken2 on raw ONT fastq reads (pre-filtering baseline).
# Usage: run_kraken2.sh <fastq_file> <output_dir> <kraken2_db> [threads]
set -euo pipefail
fastq_file=$1
outdir=$2
db=$3
threads=${4:-8}

# Extracts filename without folder path or extension (e.g., 'sample1' from 'data/sample1.fastq')
sample_id=$(basename "${fastq_file}" .fastq)
sample_id=$(basename "${sample_id}" .gz)

mkdir -p "${outdir}"
kraken2 \
  --db ${db} \
  --threads "${threads}" \
  --output "${outdir}/${sample_id}.kraken" \
  --report "${outdir}/${sample_id}_report.txt" \
  ${fastq_file}
