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
mv e_decipiens_dir/ncbi_dataset/data/GCA_014182575.1/*.fna E_decipiens.fa


# Download E. virginea (GCA_014182375.1)"
datasets download genome accession GCA_014182375.1 --filename e_virginea.zip
unzip e_virginea.zip -d e_virginea_dir
mv e_virginea_dir/ncbi_dataset/data/GCA_014182375.1/*.fna E_virginea.fa

# Cleanup zip and intermediate directory artifacts
rm -rf *_dir *.zip
