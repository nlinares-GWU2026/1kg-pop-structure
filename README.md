# Population Structure Analysis - 1000 Genomes Phase 3

A reproducible population genetics pipeline for analyzing human genetic diversity using the 1000 Genomes Project Phase 3 dataset. This project reconstructs global population structure from genome-wide SNP data using Principal Component Analysis (PCA) and sparse Non-negative Matrix Factorization, which are two primary methods in modern population genomics. 

--- 

## Background

The 1000 Genomes Project sequenced 2,504 individuals from 26 populations across 5 continental superpopulations (AFR — African, EUR — European, EAS — East Asian, SAS — South Asian, AMR — Admixed American). By analyzing patterns of common genetic variation across the genome, we can detect population structure and how the signatures of human migration, isolation, and admixture are written into our DNA. 

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
│   |-- 05_admixture.sh        <- ADMIXTURE ancestry estimation (see note in Troubleshooting)
|   |-- 05_admixture.R         <- sNMF ancestry estimation via LEA R package
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
RAW VCF (1000 Genomes)                          ->

Step 1: Data Acquisition
curl from IGSR FTP server                       ->

Step 2: Quality Control (PLINK)
MAF > 5%, geno <5%, HWE p > 1e-6
1,097,204 -> 68,355 SNPs                        ->

Step 3: LD Pruning (PLINK)
Two pruned datasets created:           
chr22_pruned (colon IDs for PCA)
chr22_pruned_admix (underscore IDs) for sNMF    ->

Step 4a: PCA
    PLINK --pca 10
    chr22_pruned

        +

Step 4b: sNMF (LEA)
    K = 2-8, cross-entropy
    chr22_pruned_admix                          ->

Step 5: Visualization
Python (matplotlib, seaborn)
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

### 6. Run PCA and ancestry estimation

```bash
bash scripts/04_pca.sh
bash scripts/05_admixture.sh   # see Troubleshooting if on WSL2
Rscript scripts/05_admixture.R
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
- **LD pruning:** 68,355 → 7,751 independent SNPs retained (r² < 0.2)
- **PCA:** PC1 and PC2 encapsulate global geography, separating AFR, EUR, EAS, SAS, and AMR superpopulations. PC1 and PC2 explain 72.8% of genetic variance; expected to recapitulate global geography separating the 5 superpopulations.
- **sNMF:** Best-fit K = 5 (cross-entropy = 0.67164), corresponding to the 5 superpopulations (AFR, EUR, EAS, SAS, AMR)

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
  - bcftools
  - r-base=4.3
  - r-biocmanager
  # LEA installed inside R: BiocManager::install('LEA')
  # ADMIXTURE v1.3.0 excluded — segfaults on WSL2 (exit code 139)
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
5. Frichot, E., & François, O. (2015). LEA: An R package for landscape and ecological association studies. *Methods in Ecology and Evolution*, 6(8), 925–929.

---

## Troubleshooting

**ADMIXTURE v1.3.0 on WSL2:** Produces a segmentation fault (exit code 139) immediately after reading genotype data on WSL2. This occurs even after `ulimit -s unlimited` and is caused by incompatibility between the statically linked binary and WSL2's memory model. Ancestry estimation was performed using R's LEA package (`snmf()`) instead. It is a mathematically equivalent algorithm producing directly comparable Q matrices and cross-entropy values. 

**SNP IDs with colons:** ADMIXTURE 1.3.0 also crashes silently when SNP IDs contain colons (for example: `NN:XXXXXXXX:A:C`). Two separate LD-pruned datasets were created: `chr22_pruned` with colon IDs for PCA, and `chr22_pruned_admix` with underscore IDs for sNMF.

---

## Author 
Nicole Linares | George Washington University
MS in Health Data Science - Milken Institute School of Public Health
GitHub: [@nlinares-GWU2026](https://github.com/nlinares-GWU2026)
