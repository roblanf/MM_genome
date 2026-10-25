#!/usr/bin/env bash
# Script: 04_parental_assignment/scripts/parent_scoring.sh
# Purpose: Generate a unified parental assignment table with both Target and Query coverages, identities, ratios, and deltas.

set -euo pipefail

source config.sh

paf_dir="04_parental_assignment/results"
out_dir="${paf_dir}/parent_mapping"
mkdir -p "${out_dir}"

hap1_hap2_paf="03_assembly/results/h1h2/hap1_vs_hap2.paf"
combined_table="${out_dir}/combined_parental_comparison.md"

# Establish combined_table structure:
{
    echo "| Contig | *decipiens* chromosome match | *virginea* chromosome match | *decipiens* Target Cov | *virginea* Target Cov | *decipiens* Query Cov | *virginea* Query Cov | *decipiens* Identity | *virginea* Identity | Cov Ratio | Δ Identity | Hap Match |"
    echo "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
} > "${combined_table}"

# Function adapted to add rows to combined_table
append_rows() {
    local decip_paf="$1"
    local virgin_paf="$2"
    local hap_label="$3"

    awk -v hap="$hap_label" '
    # Parse hap1_vs_hap2.paf using aggregate total length for robust matching like I did in h1/h2 alignment dotplot Rscript
    NR == FNR {
        q = $1;
        t = $6;
        alen = $11;

        # Track best Hap2 for each Hap1
        sum_h1_h2[q, t] += alen;
        if (sum_h1_h2[q, t] > max_h1_h2[q]) {
            max_h1_h2[q] = sum_h1_h2[q, t];
            best_h2_for_h1[q] = t;
        }

        # Track best Hap1 for each Hap2
        sum_h2_h1[t, q] += alen;
        if (sum_h2_h1[t, q] > max_h2_h1[t]) {
            max_h2_h1[t] = sum_h2_h1[t, q];
            best_h1_for_h2[t] = q;
        }
        next;
    }

    function load_paf(filename,   is_primary, as, i, a) {
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
                if ($10 > best_matches[$1]) {
                    best_matches[$1] = $10;
                    chrom[$1] = $6;
                    t_len[$1] = $7;    # Target Chromosome Length (Col 7)
                    q_len[$1] = $2;    # Query Contig Length (Col 2)
                    aln_len[$1] = $11;  # Alignment Block Length for Identity (Col 11)
                }
            }
        }
        close(filename);
    }

    BEGIN {
        decip_paf = ARGV[1];
        virgin_paf = ARGV[2];

        ARGV[1] = ""; ARGV[2] = "";

        load_paf(decip_paf);

        # Store decipiens values
        for (id in best_matches) {
            d_chrom[id] = chrom[id];
            d_matches[id] = best_matches[id];
            d_tlen[id] = t_len[id];
            d_qlen[id] = q_len[id];
            d_aln[id] = aln_len[id];

            # Coverages
            d_tcov[id] = (d_tlen[id] > 0) ? (d_matches[id] / d_tlen[id]) * 100 : 0;
            d_qcov[id] = (d_qlen[id] > 0) ? (d_matches[id] / d_qlen[id]) * 100 : 0;

            # Identity
            d_id[id] = (d_aln[id] > 0) ? (d_matches[id] / d_aln[id]) * 100 : 0;
        }

        delete best_matches; delete chrom; delete t_len; delete q_len; delete aln_len;

        load_paf(virgin_paf);

        for (id in d_matches) {
            v_chrom = (id in chrom) ? chrom[id] : "N/A";
            v_matches = (id in best_matches) ? best_matches[id] : 0;
            v_tlen = (id in t_len) ? t_len[id] : 1;
            v_qlen = (id in q_len) ? q_len[id] : d_qlen[id];
            v_aln = (id in aln_len) ? aln_len[id] : 1;

            v_tcov = (v_tlen > 0) ? (v_matches / v_tlen) * 100 : 0;
            v_qcov = (v_qlen > 0) ? (v_matches / v_qlen) * 100 : 0;
            v_id = (v_aln > 0) ? (v_matches / v_aln) * 100 : 0;

            # Bold highlight higher Target Coverage
            d_tcov_str = sprintf("%.2f%%", d_tcov[id]);
            v_tcov_str = sprintf("%.2f%%", v_tcov);
            if (d_tcov[id] >= v_tcov) {
                d_tcov_str = "**" d_tcov_str "**";
            } else {
                v_tcov_str = "**" v_tcov_str "**";
            }

            # Bold highlight higher Query Coverage
            d_qcov_str = sprintf("%.2f%%", d_qcov[id]);
            v_qcov_str = sprintf("%.2f%%", v_qcov);
            if (d_qcov[id] >= v_qcov) {
                d_qcov_str = "**" d_qcov_str "**";
            } else {
                v_qcov_str = "**" v_qcov_str "**";
            }

            # Bold highlight higher Identity
            d_id_str = sprintf("%.2f%%", d_id[id]);
            v_id_str = sprintf("%.2f%%", v_id);
            if (d_id[id] >= v_id) {
                d_id_str = "**" d_id_str "**";
            } else {
                v_id_str = "**" v_id_str "**";
            }

            cov_ratio = (v_tcov > 0) ? (d_tcov[id] / v_tcov) : 0;
            id_delta = d_id[id] - v_id;

            # Get haplotype matches, fallback "-" if no match
            if (hap == "Hap1") {
                hap_match_val = (id in best_h2_for_h1) ? best_h2_for_h1[id] : "-";
            } else {
                hap_match_val = (id in best_h1_for_h2) ? best_h1_for_h2[id] : "-";
            }

            printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %.2f | %+.2f%% | %s |\n",
                id, d_chrom[id], v_chrom,
                d_tcov_str, v_tcov_str, d_qcov_str, v_qcov_str, d_id_str, v_id_str, cov_ratio, id_delta, hap_match_val;
        }
    }' "${hap1_hap2_paf}" "${decip_paf}" "${virgin_paf}" | sort -k3,3 >> "${combined_table}"
}

# Add rows for haplotype 1 contigs
append_rows "${paf_dir}/hap1_vs_decipiens.paf" "${paf_dir}/hap1_vs_virginea.paf" "Hap1"

# Haplotype 2 Table
append_rows "${paf_dir}/hap2_vs_decipiens.paf" "${paf_dir}/hap2_vs_virginea.paf" "Hap2"
