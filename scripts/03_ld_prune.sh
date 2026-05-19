#!/bin/bash
# Step 3: LD pruning

# Assign unique IDs (chr:pos:ref:alt) to all variants
plink \
  --bfile results/plink/chr22_qc \
  --set-missing-var-ids @:#:\$1:\$2 \
  --make-bed \
  --out results/plink/chr22_qc_ids

# Identify independent SNPs
plink \
  --bfile results/plink/chr22_qc_ids \
  --indep-pairwise 50 10 0.2 \
  --out results/plink/pruning2

# Extract pruned SNP set
plink \
  --bfile results/plink/chr22_qc_ids \
  --extract results/plink/pruning2.prune.in \
  --make-bed \
  --out results/plink/chr22_pruned

echo "LD pruning complete. $(wc -l < results/plink/pruning2.prune.in) SNPs retained."
