#!/bin/bash
# Step 3: LD pruning

# --- For PCA (uses colon-separated IDs: chr:pos:ref:alt) ---
# Assign unique IDs with colons
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

# Extract pruned SNP set for PCA
plink \
  --bfile results/plink/chr22_qc_ids \
  --extract results/plink/pruning2.prune.in \
  --make-bed \
  --out results/plink/chr22_pruned

echo "PCA-ready dataset: $(wc -l < results/plink/pruning2.prune.in) SNPs retained."

# --- For ADMIXTURE/sNMF (requires underscore-separated IDs) ---
# NOTE: ADMIXTURE 1.3.0 crashes silently on colon-containing SNP IDs.
# A separate pruned dataset with underscore IDs is created for sNMF.

# Assign unique IDs with underscores: chr_pos_ref_alt
plink \
  --bfile results/plink/chr22_qc \
  --set-missing-var-ids @_#_\$1_\$2 \
  --make-bed \
  --out results/plink/chr22_qc_ids_admix

# Re-run pruning on underscore-ID dataset
plink \
  --bfile results/plink/chr22_qc_ids_admix \
  --indep-pairwise 50 10 0.2 \
  --out results/plink/pruning_admix

# Extract pruned SNP set for sNMF
plink \
  --bfile results/plink/chr22_qc_ids_admix \
  --extract results/plink/pruning_admix.prune.in \
  --make-bed \
  --out results/plink/chr22_pruned_admix

echo "sNMF-ready dataset: $(wc -l < results/plink/pruning_admix.prune.in) SNPs retained."
