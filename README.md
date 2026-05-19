# Population Structure Analysis - 1000 Genomes Phase 3

A reproducible population genetics pipeline for analyzing human genetic diversity using the 1000 Genomes Project Phase 3 dataset. This project reconstructs global population structure from genome-wide SNP data using Principal Component Analysis (PCA) and ADMIXTURE, which are two primary methods in modern population genomics. 

--- 

## Background

The 1000 Genomes Project sequenced 2,504 individuals from 26 populations across 5 contintental superpopulations (Africa-AFR, European-EUR, East Asian-EAS, South Asian-SAS, AMR-Admixed American). By analyzing patterns of common genetic variation across the genome, we can detect population structure and how the signatures of human migration, isolation, and admixture are written into our DNA. 

This project demonstrates a complete population structure workflow applicable to any large-scale SNP dataset. 

---

## Project Structure

```
1kg-pop-structure/
|-- README.md                  <- you are here
|-- WRITEUP.md                 <- full methods, results, and biological interpretation
|-- environment.yaml           <- conda environment (reproducible setup)
|-- config/
│   |-- samples.txt            <- sample IDs and population labels
|-- data/
│   |-- raw/                   <- VCF files (not tracked by git — see .gitignore)
│   |-- integrated_call_samples_v3.20130502.ALL.panel  <- population metadata
|-- scripts/
│   |-- 01_download.sh         <- data acquisition from IGSR FTP
│   |-- 02_qc_filter.sh        <- PLINK quality control and filtering
│   |-- 03_ld_prune.sh         <- LD pruning for independent SNPs
│   |-- 04_pca.sh              <- PCA via PLINK
│   |-- 05_admixture.sh        <- ADMIXTURE ancestry estimation
│   |-- 06_visualize.py        <- Python visualization (PCA + ADMIXTURE plots)
|-- workflow/
│   |-- Snakefile              <- Snakemake pipeline (end-to-end automation)
|-- results/
│   |-- plink/                 <- PLINK binary files and QC outputs
│   |-- pca/                   <- PCA eigenvectors and eigenvalues
│   |-- admixture/             <- ADMIXTURE Q matrices and log files
│   |-- plots/                 <- final figures
|-- notebook/
    |-- population_structure_analysis.ipynb  ← interactive analysis notebook
```

---

## Pipeline Overview

```
RAW VCF (1000 Genomes)                 ->

Step 1: Data Acquisition
curl from IGSR FTP server              ->

Step 2: Quality Control (PLINK)
MAF > 5%, geno <5%, HWE p > 1e-6
1,097,204 -> 68,355 SNPs               ->

Step 3: LD Pruning (PLINK)
Window 50, step 10, r^2 < 0.2          ->

Step 4a: PCA - PLINK --pca
Step 4b: ADMIXTURE - K = 2-8, CV error 
(combined PLINK + ADMIXTURE)           ->

Step 5: Visualization
Python / ggplot2
```

---

## Quickstart

### 1. Clone the repository

```bash
git clone https://github.com/nlinares-GWU2026/1kg-pop-structure.git
cd 1kg-pop-structure
```

### 2. Create and activate the conda environment
 
```bash
conda env create -f environment.yaml
conda activate popgen
```

### 3. Download the data
 
```bash
bash scripts/01_download.sh
```

This downloads chromosome 22 VCF (~ 197 MB) and the population panel file from the IGSR FTP server. To download all the chromosomes (1-22, ~ 150 GB total), edit the script to loop over chromosomes 1-22.
