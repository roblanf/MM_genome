#!/usr/bin/env bash
# Script: 04_parental_assignment/scripts/filter_haplotigs.sh
# Purpose: Filter assembly contigs based on length, coverage, and taxonomic rules.

set -euo pipefail

assembly_fa="$1"      # e.g., 03_assembly/results/hifiasm/MM_assembly.hap1.p_ctg.fa
blob_tsv="$2"        # e.g., 03_assembly/results/blobtools/plots/hap1_stats.tsv
out_fa="$3"          # e.g., 04_parental_assignment/results/filtered_hap1.fa

awk '
BEGIN {
    # 1. Load and filter BlobTools stats
    while ((getline < b_file) > 0) {
        if ($1 == "index") continue;
        sub(/\r$/, "", $6); # Strip trailing carriage return if present

        c_id = $2;
        c_gc = $3;
        c_len = $4;
        c_cov = $5;
        c_phyl = $6;

        # Threshold rules:
        # - Length >= 1,000,000 bp (1 Mbp)
        # - Coverage between 20x and 90x (2/3x and 3x of 30x coverage target per haplotype respectively)
        # - Phylum is Streptophyta OR (Artverviricota and length >= 2,000,000 bp to include the mislabelled large contigs)

        is_length_ok = (c_len >= 1000000);
        is_cov_ok = (c_cov >= 20 && c_cov <= 90);

        is_tax_ok = 0;
        if (c_phyl == "Streptophyta") {
            is_tax_ok = 1;
        } else if (c_phyl == "Artverviricota" && c_len >= 2000000) {
            is_tax_ok = 1; # Keep the large contigs mislabelled due to summation of small viral hits
        }

        if (is_length_ok && is_cov_ok && is_tax_ok) {
            passed_contigs[c_id] = 1;
        }
    }
    close(b_file);
}
/^>/ {
    header = $0;
    sub(/^>/, "", header);
    split(header, a, /[ \t]/);
    contig_id = a[1];

    keep = (contig_id in passed_contigs);
}
{
    if (keep) {
        print $0 >> out_fasta;
    }
}
' b_file="${blob_tsv}" out_fasta="${out_fa}" "${assembly_fa}"

