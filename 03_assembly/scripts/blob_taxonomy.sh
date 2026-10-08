#!/usr/bin/env bash
# Script: 03_assembly/blob_taxonomy.sh
# Purpose: Run DIAMOND taxonomic search on Hap1 and Hap2 for BlobTools taxid assignment

set -euo pipefail

source config.sh

hap1_fa="03_assembly/results/hifiasm/MM_assembly.hap1.p_ctg.fa"
hap2_fa="03_assembly/results/hifiasm/MM_assembly.hap2.p_ctg.fa"

tax_dir="03_assembly/results/blobtools/taxonomy"
mkdir -p "${tax_dir}"

# refer to 00_databases/setup_diamond_sprot.sh to get this path.
db_path="00_databases/uniprot/uniprot_sprot.dmnd"

# BlobTools standard BLAST output format
blast_fmt="6 qseqid staxids bitscore qstart qend sstart send pident evalue length"

# Run diamond blastx for haplotype 1
diamond blastx \
    -query "${hap1_fa}" \
    -db "${db_path}" \
    -outfmt "${diamond_fmt}" \
    -evalue 1e-25 \
    -max_hsps 1 \
    -max_target_seqs 10 \
    -num_threads 16 \
    -out "${tax_dir}/hap1_diamond.out"

# Run for haplotype 2
diamond blastx \
    -query "${hap2_fa}" \
    -db "${db_path}" \
    -outfmt "${diamond_fmt}" \
    -evalue 1e-25 \
    -max_hsps 1 \
    -max_target_seqs 10 \
    -num_threads 16 \
    -out "${tax_dir}/hap2_diamond.out"
