#!/usr/bin/env bash
# Script: 04_parental_assignment/get_parents.sh
# Purpose: Download parental genomes (E. decipiens and E. virginea) from NCBI datasets

set -euo pipefail

# Create output folder for parental genomes
pa_dir="04_parental_assignment"
mkdir -p "${pa_dir}/parental_genomes"
cd "${pa_dir}/parental_genomes"

# Download E. decipiens (GCA_014182575.1)"
datasets download genome accession GCA_014182575.1 --filename e_decipiens.zip
unzip e_decipiens.zip -d e_decipiens_dir
mv e_decipiens_dir/ncbi_dataset/data/GCA_014182575.1/*.fna E_decipiens_raw.fa


# Download E. virginea (GCA_014182375.1)"
datasets download genome accession GCA_014182375.1 --filename e_virginea.zip
unzip e_virginea.zip -d e_virginea_dir
mv e_virginea_dir/ncbi_dataset/data/GCA_014182375.1/*.fna E_virginea_raw.fa

# Rename chromosome IDs to chromosome numbers for easy interpretation of matches:

# Decipiens
awk '
BEGIN {
    map["CM024609.1"] = "1";
    map["CM024610.1"] = "2";
    map["CM024611.1"] = "3";
    map["CM024612.1"] = "4";
    map["CM024613.1"] = "5";
    map["CM024614.1"] = "6";
    map["CM024615.1"] = "7";
    map["CM024616.1"] = "8";
    map["CM024617.1"] = "9";
    map["CM024618.1"] = "10";
    map["CM024619.1"] = "11";
    map["CM024632.1"] = "Pltd";
}
/^>/ {
    header = $1;
    sub(/^>/, "", header);
    if (header in map) {
        print ">" map[header];
    } else {
        print $0; # Leaves unplaced scaffolds/contigs unchanged
    }
    next;
}
{ print }
' E_decipiens_raw.fa > E_decipiens.fa

# Virginea
awk '
BEGIN {
    map["CM024514.1"] = "1";
    map["CM024515.1"] = "2";
    map["CM024516.1"] = "3";
    map["CM024517.1"] = "4";
    map["CM024518.1"] = "5";
    map["CM024519.1"] = "6";
    map["CM024520.1"] = "7";
    map["CM024521.1"] = "8";
    map["CM024522.1"] = "9";
    map["CM024523.1"] = "10";
    map["CM024524.1"] = "11";
    map["CM024525.1"] = "Pltd";
}
/^>/ {
    header = $1;
    sub(/^>/, "", header);
    if (header in map) {
        print ">" map[header];
    } else {
        print $0;
    }
    next;
}
{ print }
' E_virginea_raw.fa > E_virginea.fa

# Cleanup zip and intermediate directory artifacts
rm -rf *_dir *.zip
