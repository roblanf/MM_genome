#!/usr/bin/env bash
# Script: 04_parental_assignment/scripts/parent_scoring.sh
# Purpose: Generate a unified side-by-side parental assignment table with coverage, identity, and ratios.

set -euo pipefail

source config.sh

paf_dir="04_parental_assignment/results"
out_dir="${paf_dir}/parent_mapping"
mkdir -p "${out_dir}"

out_table="${out_dir}/hap1_parental_comparison.md"

# Run a combined awk script to process both PAFs and merge them into one table
awk '
# Function to process a PAF file into associative arrays
function load_paf(filename,   is_primary, as, i, a, key) {
    while ((getline < filename) > 0) {
        is_primary = 0;
        for(i=13; i<=NF; i++) {
            if($i == "tp:A:P") { is_primary = 1; }
        }
        if (is_primary) {
            as = $10;
            for(i=13; i<=NF; i++) {
                if($i ~ /^AS:i:/) { split($i, a, ":"); as = a[3]; }
            }
            # Keep track of best alignment per contig based on total matches
            if ($10 > best_matches[$1]) {
                best_matches[$1] = $10;
                chrom[$1] = $6;
                q_len[$1] = $2;
                aln_len[$1] = $11; # Column 11: alignment block length for identity
            }
        }
    }
    close(filename);
}

BEGIN {
    decip_paf = ARGV[1];
    virgin_paf = ARGV[2];

    # Remove them from ARGV so awk doesn"t try to parse them as standard inputs
    ARGV[1] = ""; ARGV[2] = "";

    load_paf(decip_paf);

    # Save decipiens data into temporary arrays
    for (id in best_matches) {
        d_chrom[id] = chrom[id];
        d_matches[id] = best_matches[id];
        d_qlen[id] = q_len[id];
        d_aln[id] = aln_len[id];
        d_cov[id] = (d_matches[id] / d_qlen[id]) * 100;
        d_id[id] = (d_aln[id] > 0) ? (d_matches[id] / d_aln[id]) * 100 : 0;
    }

    # Clear arrays for virginea already recorded by decipiens
    delete best_matches; delete chrom; delete q_len; delete aln_len;

    load_paf(virgin_paf);

    # Print Markdown Header
    print "| S1 ID | decip ID | virgin ID | decip Matches | virgin Matches | decip Coverage | virgin Coverage | decip Identity | virgin Identity |";
    print "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |";

    # Combine and output for every contig found in decipiens or virginea
    for (id in d_matches) {
        v_chrom = (id in chrom) ? chrom[id] : "N/A";
        v_matches = (id in best_matches) ? best_matches[id] : 0;
        v_qlen = (id in q_len) ? q_len[id] : d_qlen[id];
        v_aln = (id in aln_len) ? aln_len[id] : 1;
        v_cov = (v_qlen > 0) ? (v_matches / v_qlen) * 100 : 0;
        v_id = (v_aln > 0) ? (v_matches / v_aln) * 100 : 0;

        # Highlight the higher coverage between the two parents
        d_cov_str = sprintf("%.2f%%", d_cov[id]);
        v_cov_str = sprintf("%.2f%%", v_cov);
        if (d_cov[id] >= v_cov) {
            d_cov_str = "**" d_cov_str "**";
        } else {
            v_cov_str = "**" v_cov_str "**";
        }

        # Highlight the higher identity
        d_id_str = sprintf("%.2f%%", d_id[id]);
        v_id_str = sprintf("%.2f%%", v_id);
        if (d_id[id] >= v_id) {
            d_id_str = "**" d_id_str "**";
        } else {
            v_id_str = "**" v_id_str "**";
        }

        printf "| %s | %s | %s | %d | %d | %s | %s | %s | %s |\n", 
            id, d_chrom[id], v_chrom, d_matches[id], v_matches, 
            d_cov_str, v_cov_str, d_id_str, v_id_str;
    }
}' "${paf_dir}/hap1_vs_decipiens.paf" "${paf_dir}/hap1_vs_virginea.paf" | sort -k1,1 > "${out_table}"
