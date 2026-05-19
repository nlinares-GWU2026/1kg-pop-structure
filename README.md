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

### 4. Run QC filtering

```bash
bash scripts/02_qc_filter.sh
```

### 5. Run LD pruning

```bash
bash scripts/03_ld_prune.sh
```

### 6. Run PCA and ADMIXTURE

```bash
bash scripts/04_pca.sh
bash scripts/05_admixture.sh
```

### 7. Generate plots

```bash
python scripts/06_visualize.py
```

---

## Data 

| File | Source | Size |
|------|--------|------|
| `ALL.chr22...vcf.gz` | IGSR FTP | ~ 197 MB |
| `ALL.chr22...vcf.gz.tbi` | IGSR FTP | ~ 36 KB |
| `integrated_call_samples_v3...panel` | IGSR FTP | ~ 54 KB |

**Note:** Raw VCF files are not tracked by git (see `.gitignore`). Download them using `scripts/01_download.sh` or directly from:
`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/`

---

## Key Results

*To be updated as analysis completes.* 

- **QC filtering:** 1,097,204 -> 68,355 SNPs retained on chromosome 22 after MAF, HWE, and missingness filters
- **PCA:** PC1 and PC2 encapsulate global geography, separating AFR, EUR, EAS, SAS, and AMR superpopulations
- **ADMIXTURE:** Best-fit K = 5, corresponding to the 5 superpopulations; admixed American samples show expected multi-ancestry patterns

---
 
## Environment
 
```yaml
name: popgen
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - plink=1.90
  - admixture
  - bcftools
  - jupyterlab
  - pandas
  - matplotlib
  - seaborn
```
 
---

## References
 
1. 1000 Genomes Project Consortium (2015). A global reference for human genetic variation. *Nature*, 526, 68–74.
2. Patterson et al. (2006). Population structure and eigenanalysis. *PLOS Genetics*, 2(12), e190.
3. Alexander et al. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655–1664.
4. Chang et al. (2015). Second-generation PLINK. *GigaScience*, 4(1).

---

## Author 
Nicole Linares | George Washington University
MS in Health Data Science - Milken Institute School of Public Health
GitHub: [@nlinares-GWU2026](https://github.com/nlinares-GWU2026)
