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
blobtools create \
    --input "${asm_dir}/MM_assembly.hap1.p_ctg.fa" \
    --bam "${bam_dir}/hap1_ont.bam" \
    --hits "${tax_dir}/hap1_diamond.out" \
    --taxfmt diamond \
    --nodes 00_databases/uniprot/nodes.dmp \
    --names 00_databases/uniprot/names.dmp \
    --out "${out_dir}/hap1_blobDB"

# Generate BlobPlot (GC vs Coverage colored by Phylum)
blobtools plot \
    --input "${out_dir}/hap1_blobDB.json" \
    --taxlevel phylum \
    --out "${out_dir}/hap1"

# Generate contig taxonomy summary table
blobtools view \
    --input "${out_dir}/hap1_blobDB.json" \
    --out "${out_dir}/hap1"

# Repeating the steps for haplotype 2
blobtools create \
    --input "${asm_dir}/MM_assembly.hap2.p_ctg.fa" \
    --bam "${bam_dir}/hap2_ont.bam" \
    --hits "${tax_dir}/hap2_diamond.out" \
    --taxfmt diamond \
    --nodes 00_databases/uniprot/nodes.dmp \
    --names 00_databases/uniprot/names.dmp \
    --out "${out_dir}/hap2_blobDB"

blobtools plot \
    --input "${out_dir}/hap2_blobDB.json" \
    --taxlevel phylum \
    --out "${out_dir}/hap2"

blobtools view \
    --input "${out_dir}/hap2_blobDB.json" \
    --out "${out_dir}/hap2"
