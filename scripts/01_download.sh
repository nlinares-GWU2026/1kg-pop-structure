#!/bin/bash
# Step 1: Download chr22 VCF and population panel from IGSR FTP

BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

mkdir -p data/raw

curl -L -O --output-dir data/raw/ \
  "${BASE}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

curl -L -O --output-dir data/raw/ \
  "${BASE}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"

curl -L -O --output-dir data/ \
  "${BASE}/integrated_call_samples_v3.20130502.ALL.panel"

echo "Download complete."
