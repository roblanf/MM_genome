#!/usr/bin/env bash
# Script: 04_parental_assignment/scripts/bin_parental_contigs.sh
# Purpose: Classify and bin contigs into parental fastas using awk.

set -euo pipefail

source config.sh

assembly_fa="$1"
decip_paf="$2"
virg_paf="$3"
out_dir="$4"

# To run this, specify the above for hap1 files or hap2 files like so:
#bash 04_parental_assignment/scripts/bin_parental_contigs.sh \
#    03_assembly/results/hifiasm/MM_assembly.hap1.p_ctg.fa \
#    04_parental_assignment/results/hap1_vs_decipiens.paf \
#    04_parental_assignment/results/hap1_vs_virginea.paf \
#    04_parental_assignment/results/binned_hap1

mkdir -p "${out_dir}"

decip_out="${out_dir}/decipiens.fa"
virg_out="${out_dir}/virginea.fa"
ambig_out="${out_dir}/ambiguous.fa"

# Process PAFs in Awk, classify every contig ID, and store classifications in a lookup file
awk '
function load_paf(filename, is_decip,   is_primary, i) {
    while ((getline < filename) > 0) {
        is_primary = 0;
        for(i=13; i<=NF; i++) {
            if($i == "tp:A:P") { is_primary = 1; }
        }
        if (is_primary) {
            if ($10 > best_matches[$1, is_decip]) {
                best_matches[$1, is_decip] = $10;
                t_len[$1, is_decip] = $7;
                aln_len[$1, is_decip] = $11;
            }
        }
    }
    close(filename);
}

BEGIN {
    decip_paf = ARGV[1];
    virg_paf = ARGV[2];
    ARGV[1] = ""; ARGV[2] = "";

    load_paf(decip_paf, 1);
    load_paf(virg_paf, 0);

    # Get unique contig list
    for (key in best_matches) {
        split(key, parts, SUBSEP);
        contig = parts[1];
        contigs[contig] = 1;
    }

    for (q in contigs) {
        d_matches = best_matches[q, 1];
        d_tlen = t_len[q, 1];
        d_aln = aln_len[q, 1];

        v_matches = best_matches[q, 0];
        v_tlen = t_len[q, 0];
        v_aln = aln_len[q, 0];

        d_cov = (d_tlen > 0) ? (d_matches / d_tlen) * 100 : 0;
        v_cov = (v_tlen > 0) ? (v_matches / v_tlen) * 100 : 0;

        d_id = (d_aln > 0) ? (d_matches / d_aln) * 100 : 0;
        v_id = (v_aln > 0) ? (v_matches / v_aln) * 100 : 0;

        cov_ratio = (v_cov > 0) ? (d_cov / v_cov) : (d_cov > 0 ? 99 : 1);
        id_delta = d_id - v_id;

        # Thresholds: Cov Ratio >= 1.20 or ID Delta >= +5.0 -> decipiens
        #             Cov Ratio <= 0.833333 or ID Delta <= -5.0 -> virginea
        #             Otherwise -> ambiguous
        if (cov_ratio >= 1.20 || id_delta >= 5.0) {
            print q "\tdecipiens";
        } else if (cov_ratio <= 0.833333 || id_delta <= -5.0) {
            print q "\tvirginea";
        } else {
            print q "\tambiguous";
        }
    }
}' "$decip_paf" "$virg_paf" > "${out_dir}/contig_classes.tsv"

# Next Step: Stream through the FASTA file once using Awk to split sequences into respective files

# Step 2: Stream through FASTA and split securely
awk -v cfile="${out_dir}/contig_classes.tsv" \
    -v decip_file="${decip_out}" \
    -v virg_file="${virg_out}" \
    -v ambig_file="${ambig_out}" '
BEGIN {
    while ((getline < cfile) > 0) {
        class[$1] = $2;
    }
    close(cfile);
    target = ambig_file;
}
/^>/ {
    header = $0;
    sub(/^>/, "", header);
    split(header, a, /[ \t]/);
    contig_id = a[1];

    bin = (contig_id in class) ? class[contig_id] : "ambiguous";
    target = (bin == "decipiens") ? decip_file : ((bin == "virginea") ? virg_file : ambig_file);
}
{
    print $0 >> target;
}' "$assembly_fa"
