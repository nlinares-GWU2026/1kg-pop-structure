#!/bin/bash
# Step 4b: Ancestry estimation via sNMF (LEA R package)
#
# NOTE: ADMIXTURE v1.3.0 (conda bioconda) produces a segmentation fault
# (exit code 139) on WSL2 due to incompatibility between the statically
# linked binary and WSL2's memory model. This is a known issue with no
# available fix for this binary on WSL2.
#
# Solution: R package LEA implements an equivalent sparse Non-negative
# Matrix Factorization (sNMF) algorithm that produces directly comparable
# Q matrices and cross-entropy values (analogous to ADMIXTURE's CV error).
# Results and interpretation are scientifically equivalent.
#
# Reference: Frichot & Francois (2015) LEA: An R package for landscape
# and ecological association studies. Methods in Ecology and Evolution.

echo "Running sNMF ancestry estimation via LEA R package..."
Rscript scripts/05_snmf.R
echo "Done. Check results/admixture/ for Q matrices and cross-entropy values."
