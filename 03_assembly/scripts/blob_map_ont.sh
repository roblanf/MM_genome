#!/usr/bin/env bash
# Script: 03_assembly/blob_map_ont.sh
# Purpose: Map filtered ONT reads to Hap1 and Hap2 assemblies, needed for BlobTools coverage

set -euo pipefail

source config.sh

hap1_fa="03_assembly/results/hifiasm/MM_assembly.hap1.p_ctg.fa"
hap2_fa="03_assembly/results/hifiasm/MM_assembly.hap2.p_ctg.fa"

blob_dir="03_assembly/results/blobtools/mapping"
mkdir -p "${blob_dir}"

# map post-filtering reads to haplotype 1 assembly with minimap, samtools to get right format for BlobTools
# Uses -c to make .csi files not .bai for the later blob_create_and_plot.sh script
minimap2 -ax map-ont -t 64 "${hap1_fa}" "${filtered_fastq}" | \
    samtools sort -@ 32 -m 2G -T "${blob_dir}/hap1_tmp" -o "${blob_dir}/hap1_ont.bam" -
samtools index -c -@ 32 "${blob_dir}/hap1_ont.bam"

# same idea for haplotype 2 assembly
minimap2 -ax map-ont -t 64 "${hap2_fa}" "${filtered_fastq}" | \
    samtools sort -@ 32 -m 2G -T "${blob_dir}/hap2_tmp" -o "${blob_dir}/hap2_ont.bam" -
samtools index -c -@ 32 "${blob_dir}/hap2_ont.bam"
