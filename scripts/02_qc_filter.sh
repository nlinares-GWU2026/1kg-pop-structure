#!/bin/bash
# Step 2: QC filtering with PLINK

mkdir -p results/plink

# Part 1: Convert VCF to PLINK binary format
plink \
  --vcf data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
  --make-bed \
  --out results/plink/chr22_raw \
  --double-id \
  --allow-extra-chr \
  --biallelic-only

# Part 2: Apply QC filters
plink \
  --bfile results/plink/chr22_raw \
  --maf 0.05 \
  --geno 0.05 \
  --mind 0.05 \
  --hwe 1e-6 \
  --make-bed \
  --out results/plink/chr22_qc

echo "QC filtering complete. Check results/plink/chr22_qc.log for summary."
