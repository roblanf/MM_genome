#!/usr/bin/env bash
# Script: 04_parental_assignment/scripts/parent_scoring.sh
# Purpose: Extract best parental chromosome matches, coverage, and alignment scores silently.

set -euo pipefail

source config.sh

# Directory structure
paf_dir="04_parental_assignment/results"
out_dir="${paf_dir}/parent_mapping"
mkdir -p "${out_dir}"

# --- Haplotype 1 vs E. decipiens ---
awk '
BEGIN {
    print "| Hap1 ID | Decipiens ID | Total Matches | Coverage % | Alignment Score |";
    print "| :--- | :--- | :--- | :--- | :--- |";
}
{   #Extract alignment score tag by searching through columns, defaulting to match score count $10 if not found
    as = $10;
    #
    for(i=13; i<=NF; i++) {
        if($i ~ /^AS:i:/) { split($i, a, ":"); as = a[3]; }
     }
   # $1 is query, $6 is target chromosome from parent reference genome. Matches line adds matching bases in split blocks if alignment is broken up
    key = $1"|"$6;
    matches[key] += $10;
    t_len[key] = $7;
    score[key] = as;
}
#Score saved!

# Calculate coverage percentage to drop insignificant alignment fragments
END {
    for (k in matches) {
        split(k, names, "|");
        cov = (matches[k] / t_len[k]) * 100;
        if (cov > 1) {
            printf "| %s | %s | %d | %.2f%% | %d |\n", names[1], names[2], matches[k], cov, score[k];
        }
    } # Finally, sort by contig ID and coverage, descending, avoid duplication by using first (best) matches
}' "${paf_dir}/hap1_vs_decipiens.paf" | sort -t'|' -k1,1 -k3,3nr | awk -F'|' '!seen[$1]++' > "${out_dir}/hap1_best_vs_decipiens.md"

# --- Haplotype 1 vs E. virginea ---
awk '
BEGIN {
    print "| Hap1 ID | Virginea ID | Total Matches | Coverage % | Alignment Score |";
    print "| :--- | :--- | :--- | :--- | :--- |";
}
{
    as = $10;
    for(i=13; i<=NF; i++) { if($i ~ /^AS:i:/) { split($i, a, ":"); as = a[3]; } }
    key = $1"|"$6;
    matches[key] += $10;
    t_len[key] = $7;
    score[key] = as;
}
END {
    for (k in matches) {
        split(k, names, "|");
        cov = (matches[k] / t_len[key]) * 100;
        if (cov > 1) {
            printf "| %s | %s | %d | %.2f%% | %d |\n", names[1], names[2], matches[key], cov, score[key];
        }
    }
}' "${paf_dir}/hap1_vs_virginea.paf" | sort -t'|' -k1,1 -k3,3nr | awk -F'|' '!seen[$1]++' > "${out_dir}/hap1_best_vs_virginea.md"
