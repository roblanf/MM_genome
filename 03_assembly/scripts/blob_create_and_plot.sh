#!/usr/bin/env bash
# Script: 03_assembly/blob_create_and_plot.sh
# Purpose: Integrate the key components FASTA, BAM coverage, and DIAMOND taxonomy into BlobTools datasets and generate blobplots

set -euo pipefail

source config.sh

# Directory structure
asm_dir="03_assembly/results/hifiasm"
bam_dir="03_assembly/results/blobtools/mapping"
tax_dir="03_assembly/results/blobtools/taxonomy"
out_dir="03_assembly/results/blobtools/plots"

mkdir -p "${out_dir}"

# Create database for haplotype 1
blobtools add \
    --create \
    --threads 32 \
    --fasta "${asm_dir}/MM_assembly.hap1.p_ctg.fa" \
    --cov "${bam_dir}/hap1_ont.bam" \
    "${out_dir}/hap1_blobDB"

# Add DIAMOND taxonomy hits and dumps
blobtools add \
    --threads 32 \
    --hits "${tax_dir}/hap1_diamond.out" \
    --taxfmt diamond \
    --nodes 00_databases/uniprot/nodes.dmp \
    --names 00_databases/uniprot/names.dmp \
    "${out_dir}/hap1_blobDB"

# Generate BlobPlot (GC vs coverage, colour coded by phylum) and contig taxonomy summary table using blobtools view
blobtools view \
    --view blob \
    --format png \
    --out "${out_dir}/hap1" \
    "${out_dir}/hap1_blobDB"

# Repeating the steps for haplotype 2
blobtools add \
    --threads 32 \
    --create \
    --fasta "${asm_dir}/MM_assembly.hap2.p_ctg.fa" \
    --cov "${bam_dir}/hap2_ont.bam" \
    "${out_dir}/hap2_blobDB"

blobtools add \
    --threads 32 \
    --hits "${tax_dir}/hap2_diamond.out" \
    --taxfmt diamond \
    --nodes 00_databases/uniprot/nodes.dmp \
    --names 00_databases/uniprot/names.dmp \
    "${out_dir}/hap2_blobDB"

blobtools view \
    --view blob \
    --format png \
    --out "${out_dir}/hap2" \
    "${out_dir}/hap2_blobDB"
