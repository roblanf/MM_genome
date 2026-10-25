#!/usr/bin/env bash
# Script: 04_parental_assignment/scripts/dv_scoring.sh
# Purpose: Generate a unified parental assignment table with all metrics of interest.

set -euo pipefail

source config.sh

paf_dir="04_parental_assignment/results"
out_dir="${paf_dir}/parent_mapping"
mkdir -p "${out_dir}"

hap1_hap2_paf="03_assembly/results/h1h2/hap1_vs_hap2.paf"
h1_dec_paf="${paf_dir}/hap1_vs_decipiens.paf"
h1_vir_paf="${paf_dir}/hap1_vs_virginea.paf"
h2_dec_paf="${paf_dir}/hap2_vs_decipiens.paf"
h2_vir_paf="${paf_dir}/hap2_vs_virginea.paf"

combined_table="${out_dir}/combined_parental_comparison.md"

# Establish combined_table structure:
{
    echo "| Contig | Hap Match | *decipiens* chrom | *virginea* chrom | *decipiens* T-Cov | *virginea* T-Cov | *decipiens* ID | *virginea* ID | Cov Ratio | Δ ID | Partner *decipiens* T-Cov | Partner *virginea* T-Cov | Partner *decipiens* ID | Partner *virginea* ID | Partner Cov Ratio | Partner Δ ID |"
    echo "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
} > "${combined_table}"

# Pre-compute cross-haplotype best matches and parental metrics for all contigs - saves compute time
h1_h2_map=$(mktemp)
h2_h1_map=$(mktemp)
metrics_cache=$(mktemp)

# a) Compute Hap1 vs Hap2 best matches
awk '
{
    q = $1; t = $6; alen = $11;
    sum_h1_h2[q, t] += alen;
    sum_h2_h1[t, q] += alen;
}
END {
    for (pair in sum_h1_h2) {
        split(pair, p, SUBSEP);
        if (sum_h1_h2[pair] > max_h1[p[1]]) {
            max_h1[p[1]] = sum_h1_h2[pair];
            best_h2[p[1]] = p[2];
        }
    }
    for (q in best_h2) print q, best_h2[q] > "'"${h1_h2_map}"'";

    for (pair in sum_h2_h1) {
        split(pair, p, SUBSEP);
        if (sum_h2_h1[pair] > max_h2[p[1]]) {
            max_h2[p[1]] = sum_h2_h1[pair];
            best_h1[p[1]] = p[2];
        }
    }
    for (q in best_h1) print q, best_h1[q] > "'"${h2_h1_map}"'";
}' "${hap1_hap2_paf}"

# b) Helper function to calculate parental scores (coverage ratios, identity differences), add to metrics cache for use in table
calc_parental_scores() {
    local dec_paf="$1"
    local vir_paf="$2"

    awk '
    function load_paf(filename, is_primary, as, i, a) {
        while ((getline < filename) > 0) {
            is_primary = 0;
            for(i=13; i<=NF; i++) if($i == "tp:A:P") is_primary = 1;
            if (is_primary) {
                as = $10;
                for(i=13; i<=NF; i++) if($i ~ /^AS:i:/) { split($i, a, ":"); as = a[3]; }
                if ($10 > best_matches[$1]) {
                    best_matches[$1] = $10;
                    chrom[$1] = $6;
                    t_len[$1] = $7;
                    aln_len[$1] = $11;
                }
            }
        }
        close(filename);
    }
    BEGIN {
        decip_paf = ARGV[1];
        vir_paf = ARGV[2];
        ARGV[1] = ""; ARGV[2] = "";

        load_paf(decip_paf);
        for (id in best_matches) {
            d_chrom[id] = chrom[id];
            d_matches[id] = best_matches[id];
            d_tlen[id] = t_len[id];
            d_tcov[id] = (d_tlen[id] > 0) ? (d_matches[id] / d_tlen[id]) * 100 : 0;
            d_id[id] = (aln_len[id] > 0) ? (d_matches[id] / aln_len[id]) * 100 : 0;
        }

        delete best_matches; delete chrom; delete t_len; delete aln_len;

        load_paf(vir_paf);
        for (id in d_tcov) {
            v_chrom = (id in chrom) ? chrom[id] : "N/A";
            v_matches = (id in best_matches) ? best_matches[id] : 0;
            v_tlen = (id in t_len) ? t_len[id] : 1;
            v_tcov = (v_tlen > 0) ? (v_matches / v_tlen) * 100 : 0;
            v_id = (aln_len[id] > 0) ? (v_matches / aln_len[id]) * 100 : 0;

            # Output format: id \t d_chrom \t v_chrom \t d_tcov \t v_tcov \t d_id \t v_id
            print id, d_chrom[id], v_chrom, d_tcov[id], v_tcov, d_id[id], v_id >> "'"${metrics_cache}"'";
        }
    }' "$dec_paf" "$vir_paf"
}

# Run parental assignment metrics for Hap1 and Hap2
calc_parental_scores "${h1_dec_paf}" "${h1_vir_paf}"
calc_parental_scores "${h2_dec_paf}" "${h2_vir_paf}"

# Process and output rows using the pre-computed metrics and maps
process_hap() {
    local hap_label="$1"

    awk -v hap="$hap_label" -v h1_h2_file="${h1_h2_map}" -v h2_h1_file="${h2_h1_map}" -v cache_file="${metrics_cache}" '
    BEGIN {
        # Load cross-hap mapping
        if (hap == "Hap1") {
            while ((getline < h1_h2_file) > 0) best_map[$1] = $2;
        } else {
            while ((getline < h2_h1_file) > 0) best_map[$1] = $2;
        }
        close(h1_h2_file); close(h2_h1_file);

        # Load parental metrics cache for ALL contigs
  paf_dir="04_parental_assignment/results"
out_dir="${paf_dir}/parent_mapping"
mkdir -p "${out_dir}"

hap1_hap2_paf="03_assembly/results/h1h2/hap1_vs_hap2.paf"
h1_dec_paf="${paf_dir}/hap1_vs_decipiens.paf"
h1_vir_paf="${paf_dir}/hap1_vs_virginea.paf"
h2_dec_paf="${paf_dir}/hap2_vs_decipiens.paf"
h2_vir_paf="${paf_dir}/hap2_vs_virginea.paf"

combined_table="${out_dir}/combined_parental_comparison.md"

# 1. Establish combined_table structure (Query coverage removed, Partner metrics expanded):
{
    echo "| Contig | Hap Match | *decipiens* chrom | *virginea* chrom | *decipiens* T-Cov | *virginea* T-Cov | *decipiens* ID | *virginea* ID | Cov Ratio | Δ ID | Partner *decipiens* T-Cov | Partner *virginea* T-Cov | Partner *decipiens* ID | Partner *virginea* ID | Partner Cov Ratio | Partner Δ ID |"
    echo "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
} > "${combined_table}"

# 2. Pre-compute cross-haplotype best matches AND parental metrics for ALL contigs upfront
echo "Pre-computing cross-haplotype maps and parental metrics..."
h1_h2_map=$(mktemp)
h2_h1_map=$(mktemp)
metrics_cache=$(mktemp)

# A. Compute Hap1 vs Hap2 best matches
awk '
{
    q = $1; t = $6; alen = $11;
    sum_h1_h2[q, t] += alen;
    sum_h2_h1[t, q] += alen;
}
END {
    for (pair in sum_h1_h2) {
        split(pair, p, SUBSEP);
        if (sum_h1_h2[pair] > max_h1[p[1]]) {
            max_h1[p[1]] = sum_h1_h2[pair];
            best_h2[p[1]] = p[2];
        }
    }
    for (q in best_h2) print q, best_h2[q] > "'"${h1_h2_map}"'";

    for (pair in sum_h2_h1) {
        split(pair, p, SUBSEP);
        if (sum_h2_h1[pair] > max_h2[p[1]]) {
            max_h2[p[1]] = sum_h2_h1[pair];
            best_h1[p[1]] = p[2];
        }
    }
    for (q in best_h1) print q, best_h1[q] > "'"${h2_h1_map}"'";
}' "${hap1_hap2_paf}"

# B. Helper function to calculate parental scores for a given set of PAFs and write to metrics cache
calc_parental_scores() {
    local dec_paf="$1"
    local vir_paf="$2"
    
    awk '
    function load_paf(filename, is_primary, as, i, a) {
        while ((getline < filename) > 0) {
            is_primary = 0;
            for(i=13; i<=NF; i++) if($i == "tp:A:P") is_primary = 1;
            if (is_primary) {
                as = $10;
                for(i=13; i<=NF; i++) if($i ~ /^AS:i:/) { split($i, a, ":"); as = a[3]; }
                if ($10 > best_matches[$1]) {
                    best_matches[$1] = $10;
                    chrom[$1] = $6;
                    t_len[$1] = $7;
                    aln_len[$1] = $11;
                }
            }
        }
        close(filename);
    }
    BEGIN {
        decip_paf = ARGV[1];
        vir_paf = ARGV[2];
        ARGV[1] = ""; ARGV[2] = "";

        load_paf(decip_paf);
        for (id in best_matches) {
            d_chrom[id] = chrom[id];
            d_matches[id] = best_matches[id];
            d_tlen[id] = t_len[id];
            d_tcov[id] = (d_tlen[id] > 0) ? (d_matches[id] / d_tlen[id]) * 100 : 0;
            d_id[id] = (aln_len[id] > 0) ? (d_matches[id] / aln_len[id]) * 100 : 0;
        }

        delete best_matches; delete chrom; delete t_len; delete aln_len;

        load_paf(vir_paf);
        for (id in d_tcov) {
            v_chrom = (id in chrom) ? chrom[id] : "N/A";
            v_matches = (id in best_matches) ? best_matches[id] : 0;
            v_tlen = (id in t_len) ? t_len[id] : 1;
            v_tcov = (v_tlen > 0) ? (v_matches / v_tlen) * 100 : 0;
            v_id = (aln_len[id] > 0) ? (v_matches / aln_len[id]) * 100 : 0;

            # Output format: id \t d_chrom \t v_chrom \t d_tcov \t v_tcov \t d_id \t v_id
            print id, d_chrom[id], v_chrom, d_tcov[id], v_tcov, d_id[id], v_id >> "'"${metrics_cache}"'";
        }
    }' "$dec_paf" "$vir_paf"
}

# Run calculation for Hap1 and Hap2
calc_parental_scores "${h1_dec_paf}" "${h1_vir_paf}"
calc_parental_scores "${h2_dec_paf}" "${h2_vir_paf}"

# 3. Process and output rows using the pre-computed metrics and maps
process_hap() {
    local hap_label="$1"

    awk -v hap="$hap_label" -v h1_h2_file="${h1_h2_map}" -v h2_h1_file="${h2_h1_map}" -v cache_file="${metrics_cache}" '
    BEGIN {
        # Load cross-hap mapping
        if (hap == "Hap1") {
            while ((getline < h1_h2_file) > 0) best_map[$1] = $2;
        } else {
            while ((getline < h2_h1_file) > 0) best_map[$1] = $2;
        }
        close(h1_h2_file); close(h2_h1_file);

        # Load parental metrics cache for all contigs
        while ((getline < cache_file) > 0) {
            id = $1;
            dc[id] = $2;
            vc[id] = $3;
            dtc[id] = $4;
            vtc[id] = $5;
            did[id] = $6;
            vid[id] = $7;
        }
        close(cache_file);

        # Iterate through all loaded contigs that belong to this haplotype group
        for (id in dtc) {
            d_chrom = dc[id];
            v_chrom = vc[id];
            d_tcov = dtc[id];
            v_tcov = vtc[id];
            d_id = did[id];
            v_id = vid[id];

            cov_ratio = (v_tcov > 0) ? (d_tcov / v_tcov) : 0;
            id_delta = d_id - v_id;

            hap_match = (id in best_map) ? best_map[id] : "-";

            # Partner metrics lookup and derived metrics
            p_dtc = (hap_match != "-" && hap_match in dtc) ? dtc[hap_match] : 0;
            p_vtc = (hap_match != "-" && hap_match in vtc) ? vtc[hap_match] : 0;
            p_did = (hap_match != "-" && hap_match in did) ? did[hap_match] : 0;
            p_vid = (hap_match != "-" && hap_match in vid) ? vid[hap_match] : 0;

            p_cov_ratio = (p_vtc > 0) ? (p_dtc / p_vtc) : 0;
            p_id_delta = p_did - p_vid;

            printf "| %s | %s | %s | %s | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f | %+.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f | %+.2f%% |\n",
                id, hap_match, d_chrom, v_chrom,
                d_tcov, v_tcov, d_id, v_id, cov_ratio, id_delta,
                p_dtc, p_vtc, p_did, p_vid, p_cov_ratio, p_id_delta;
        }
    }' | sort -k1,1 >> "${combined_table}"
}

process_hap "Hap1"
process_hap "Hap2"

# Clean up temporary files
rm -f "${h1_h2_map}" "${h2_h1_map}" "${metrics_cache}"      while ((getline < cache_file) > 0) {
            id = $1;
            dc[id] = $2;
            vc[id] = $3;
            dtc[id] = $4;
            vtc[id] = $5;
            did[id] = $6;
            vid[id] = $7;
        }
        close(cache_file);

        # Iterate through all loaded contigs that belong to this haplotype group
        for (id in dtc) {
            d_chrom = dc[id];
            v_chrom = vc[id];
            d_tcov = dtc[id];
            v_tcov = vtc[id];
            d_id = did[id];
            v_id = vid[id];

            cov_ratio = (v_tcov > 0) ? (d_tcov / v_tcov) : 0;
            id_delta = d_id - v_id;

            hap_match = (id in best_map) ? best_map[id] : "-";

            # Partner metrics lookup and derived metrics
            p_dtc = (hap_match != "-" && hap_match in dtc) ? dtc[hap_match] : 0;
            p_vtc = (hap_match != "-" && hap_match in vtc) ? vtc[hap_match] : 0;
            p_did = (hap_match != "-" && hap_match in did) ? did[hap_match] : 0;
            p_vid = (hap_match != "-" && hap_match in vid) ? vid[hap_match] : 0;

            p_cov_ratio = (p_vtc > 0) ? (p_dtc / p_vtc) : 0;
            p_id_delta = p_did - p_vid;

            printf "| %s | %s | %s | %s | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f | %+.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f%% | %.2f | %+.2f%% |\n",
                id, hap_match, d_chrom, v_chrom,
                d_tcov, v_tcov, d_id, v_id, cov_ratio, id_delta,
                p_dtc, p_vtc, p_did, p_vid, p_cov_ratio, p_id_delta;
        }
    }' | sort -k1,1 >> "${combined_table}"
}

process_hap "Hap1"
process_hap "Hap2"

# Clean up temporary files
rm -f "${h1_h2_map}" "${h2_h1_map}" "${metrics_cache}"
