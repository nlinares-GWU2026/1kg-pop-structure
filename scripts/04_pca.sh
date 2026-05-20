#!/bin/bash
# Step 4a: Principal Component Analysis

mkdir -p results/pca

plink \
  --bfile results/plink/chr22_pruned \
  --pca 10 \
  --out results/pca/chr22_pca

echo "PCA complete."
echo "Eigenvalues:"
cat results/pca/chr22_pca.eigenval
