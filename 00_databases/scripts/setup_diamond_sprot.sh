#!/usr/bin/env bash
# Script: 00_databases/setup_diamond_sprot.sh
# Purpose: Download UniProt Swiss-Prot and NCBI taxonomy mapping, then build a taxonomy-enabled DIAMOND DB

set -euo pipefail

mkdir -p 00_databases/uniprot
cd 00_databases/uniprot

# 1. Download Swiss-Prot FASTA
wget -N --no-check-certificate https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip -f uniprot_sprot.fasta.gz

# 2. Download NCBI taxonomy mapping files via HTTPS
wget -N --no-check-certificate https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
wget -N --no-check-certificate https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

# 3. Extract all lineage and node .dmp files so they're available for BlobToolsv2
tar -xzf new_taxdump.tar.gz

# 4. Build binary DIAMOND database with integrated taxonomy (no hyphen in 'makedb')
diamond makedb \
    --in uniprot_sprot.fasta \
    -d uniprot_sprot \
    --taxonmap prot.accession2taxid.gz \
    --taxonnodes nodes.dmp \
    --taxonnames names.dmp

cd ../..
