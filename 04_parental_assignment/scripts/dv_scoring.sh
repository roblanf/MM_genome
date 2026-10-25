#!/usr/bin/env bash
# Script: 04_parental_assignment/scripts/dv_scoring.sh
# Purpose: Generate a unified parental assignment table with robust pre-computation.

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

# 1. Establish combined_table structure:
{
    echo "| Contig | Hap Match | *decipiens* chrom | *virginea* chrom | *decipiens* T-Cov | *virginea* T-Cov | *decipiens* ID | *virginea* ID | Cov Ratio | Δ ID | Partner *decipiens* T-Cov | Partner *virginea* T-Cov | Partner *decipiens* ID | Partner *virginea* ID | Partner Cov Ratio | Partner Δ ID |"
    echo "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
} > "${combined_table}"

# 2. Pre-compute cross-haplotype maps using sort/awk pipeline (Zero loops, zero quote errors)
echo "Pre-computing cross-haplotype maps and parental metrics..."
h1_h2_map=$(mktemp)
h2_h1_map=$(mktemp)
metrics_cache=$(mktemp)

awk '{sum[$1 "\t" $6] += $11} END {for (k in sum) print k, sum[k]}' "${hap1_hap2_paf}" \
    | sort -k1,1 -k3,3nr \
    | awk '!seen[$1]++ {print $1, $2}' > "${h1_h2_map}"

awk '{sum[$1 "\t" $6] += $11} END {for (k in sum) print k, sum[k]}' "${hap1_hap2_paf}" \
    | sort -k2,2 -k3,3nr \
    | awk '!seen[$2]++ {print $2, $1}' > "${h2_h1_map}"

# 3. Helper function to calculate parental scores
calc_parental_scores() {
    local dec_paf="$1"
    local vir_paf="$2"
    export metrics_file="${metrics_cache}"

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
        m_file = ENVIRON["metrics_file"];

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

            print id, d_chrom[id], v_chrom, d_tcov[id], v_tcov, d_id[id], v_id >> m_file;
        }
    }' "$dec_paf" "$vir_paf"
}

calc_parental_scores "${h1_dec_paf}" "${h1_vir_paf}"
calc_parental_scores "${h2_dec_paf}" "${h2_vir_paf}"

# 4. Process and output rows
process_hap() {
    local hap_label="$1"
    export h1_h2_file="${h1_h2_map}"
    export h2_h1_file="${h2_h1_map}"
    export cache_file="${metrics_cache}"

    awk -v hap="$hap_label" '
    BEGIN {
        m1 = ENVIRON["h1_h2_file"];
        m2 = ENVIRON["h2_h1_file"];
        c_file = ENVIRON["cache_file"];

        if (hap == "Hap1") {
            while ((getline < m1) > 0) best_map[$1] = $2;
        } else {
            while ((getline < m2) > 0) best_map[$1] = $2;
        }
        close(m1); close(m2);

        while ((getline < c_file) > 0) {
            id = $1;
            dc[id] = $2;
            vc[id] = $3;
            dtc[id] = $4;
            vtc[id] = $5;
            did[id] = $6;
            vid[id] = $7;
        }
        close(c_file);

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

rm -f "${h1_h2_map}" "${h2_h1_map}" "${metrics_cache}"
echo "Unified comparison report generated successfully at ${combined_table}!"
