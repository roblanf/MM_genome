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
taxdump_dir="00_databases/uniprot"

mkdir -p "${out_dir}"
mkdir -p "${out_dir}/hap1"
mkdir -p "${out_dir}/hap2"

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
    --taxdump "${taxdump_dir}" \
# Uncomment these to test my bitscore cutoff - see which hits are less confident
#    --bitscore 250 \
#    --replace
    "${out_dir}/hap1_blobDB"

# Generate BlobPlot (GC vs coverage, colour coded by phylum)
blobtools view \
    --plot \
    --view blob \
    --format png \
    --out "${out_dir}/hap1" \
    "${out_dir}/hap1_blobDB"

# Repeating the steps for haplotype 2
blobtools add \
    --create \
    --threads 32 \
    --fasta "${asm_dir}/MM_assembly.hap2.p_ctg.fa" \
    --cov "${bam_dir}/hap2_ont.bam" \
    "${out_dir}/hap2_blobDB"

blobtools add \
    --threads 32 \
    --hits "${tax_dir}/hap2_diamond.out" \
    --taxdump "${taxdump_dir}" \
#    --bitscore 250 \
#    --replace
    "${out_dir}/hap2_blobDB"

blobtools view \
    --plot \
    --view blob \
    --format png \
    --out "${out_dir}/hap2" \
    "${out_dir}/hap2_blobDB"

# Making tables to investigate contigs in each haplotype assembly:
# These will show contig IDs, lengths, GC%, coverage and phylum identified. 
blobtools filter \
    --table 03_assembly/results/blobtools/plots/hap1_stats.tsv \
    03_assembly/results/blobtools/plots/hap1_blobDB

blobtools filter \
    --table 03_assembly/results/blobtools/plots/hap2_stats.tsv \
    03_assembly/results/blobtools/plots/hap2_blobDB

# See all the stats with this line:
# less -S 03_assembly/results/blobtools/plots/hap1_stats.tsv

# Or sort by specific columns, e.g. sequence length high to low:
# sort -k4,4nr -t$'\t' 03_assembly/results/blobtools/plots/hap1_stats.tsv | less -S
