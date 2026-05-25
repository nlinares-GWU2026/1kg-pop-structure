# Population Structure Analysis - Writeup

**Project:** Population Structure Analysis of the 1000 Genomes Phase 3 Dataset
**Author:** Nicole Linares
**Date:** May 2026
**Repository:** https://github.com/nlinares-GWU2026/1kg-pop-structure

---

## 1. Motivation

Understanding population structure - the non-random distribution of genetic variation across human groups - is crucial to developing and understanding modern population genetics. As humans migrated out of Africa roughly 60,000-70,000 years ago and spread across the globe, geographically isolated populations accumulated distinct patterns of genetic variation through genetic drift, natural selection, and unique mutation histories. These differences are subtle at the individual level but statistically detectable across the genome. 

Population structure and analysis have direct practical importance beyond historical inference. In genome-wide association studies (GWAS), failure to account for population stratification leads to misleading associations between genetic variants and traits, like false positives driven by correlation between ancestry and disease prevalence rather than biology. The methods developed in this project (PCA and sNMF-based ancestry estimation) are standard tools used to detect and correct for population stratification in large-scale genetic studies. 

This project uses genome-wide SNP data from the 1000 Genomes Project Phase 3 release to reconstruct population structure across 26 globally distributed human populations, with the goals of: 1) visualizing ancestry-driven genetic clustering via PCA, and 2) estimating ancestral population proportions via sparse Non-negative Matrix Factorization (sNMF) using the R package LEA.

---

## 2. Data

### 2.1 The 1000 Genomes Project

The 1000 Genomes Project is an international effort to characterize human genetic variation at the population level. Phase 3, completed in 2015 and described in the 1000 Genomes Project Consortium (2015), represents the final and most comprehensive release of the project. It includes whole-genome sequencing data from **2,504 individuals** across **26 populations** grouped into **5 superpopulations**:

| Superpopulation Code | Description |
|----------------------|-------------|
| AFR | African |
| EUR | European |
| EAS | East Asian |
| SAS | South Asian |
| AMR | Admixed American |

All genomes were sequenced and aligned to the **GRCh37 (hg19)** human reference genome. Variant calling utilized a multi-sample pipeline that combined low-coverage whole-genome sequencing, exome sequencing, and SNP array data, resulting in a highly accurate and comprehensive genotype set. 

### 2.2 Data Access

Data were downloaded from the IGSR (International Genome Sample Resource) FTP server hosted by the European Bioinformatics Institute: 

```
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
```

The following files were used:

| File | Size | Purpose |
|------|------|---------|
| `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | 197 MB | Genotype VCF for chromosome 22 |
| `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi` | 36 KB | Tabix index for VCF |
| `integrated_call_samples_v3.20130502.ALL.panel` | 54 KB | Sample population labels |

This analysis begins with chromosome 22 as a proof-of-concept. Chromosome 22 is the smallest human autosome and provides a computationally workable starting point before scaling to the full genome (chromosomes 1-22). 

### 2.3 File Format Notes

VCF (Variant Call Format) files store genotype data as a matrix of variants * individuals. Each row represents one genomic position and contains the reference allele, alternate allele(s), and genotype calls for all individuals coded as `0/0` (homozygous reference), `0/1` (heterozygous), or `1/1` (homozygous alternate). The `.tbi` tabix index allows software to jump to any genomic coordinate without reading the entire file. 

---

## 3. Methods

### 3.1 Software and Environment

All analyses were conducted in a conda virtual environment (`popgen`) running on Ubuntu 24 (via Windows Subsystem for Linux). Key software:

| Tool | Version | Purpose |
|------|---------|---------|
| PLINK | v1.9.0-b.8 (Oct 2024) | QC filtering, LD pruning, PCA |
| LEA (R Package) | 3.24.0 | Ancestry estimation via sNMF |
| Python | 3.11 | Visualization |
| bcftools | 1.23.1 | VCF inspection |
| R | 4.6.0 | sNMF analysis | 

### 3.2 Step 1 - Data Acquisition

VCF files and the population panel were downloaded using `curl` from the IGSR FTP server into the `data/raw/` and `data/` directories, respectively. File integrity was verified by inspecting the VCF header with `bcftools view -h`, which confirmed the correct reference genome (GRCh37), pipeline source (1000GenomesPhase3Pipeline), and file date. 

### 3.3 Step 2 - Quality Control and Filtering

Raw VCF data was first converted from VCF format to PLINK binary format (`.bed/.bim/.fam`) using: 

```bash
plink \
  --vcf ALL.chr22...vcf.gz \
  --make-bed \
  --out chr22_raw \
  --double-id \
  --allow-extra-chr \
  --biallelic-only
```

The `--biallelic-only` flag excluded multi-allelic variants, which are less easily managed for standard population genetics methods. The `--double-id` flag was required because the VCF contains only a single sample identifier per individual, while PLINK expects both a family ID and an individual ID. 

Quality control filters were then applied in a single PLINK pass:

```bash
plink \
  --bfile chr22_raw \
  --maf 0.05 \
  --geno 0.05 \
  --mind 0.05 \
  --hwe 1e-6 \
  --make-bed \
  --out chr22_qc
```

**Filters applied and rationale:**

- **`--maf 0.05` (Minor Allele Frequency > 5%):** Rare variants contribute little information to population structure analysis and increase noise. Only SNPs where the minor allele is present in at least 5% of chromosomes in the sample were retained.

- **`--geno 0.05` (Genotype missingness < 5% per SNP):** SNPs with high rates of missing data often reflect technical failures in sequencing or genotype calling rather than true biological variation. SNPs missing in more than 5% of individuals were excluded.

- **`--mind 0.05` (Individual missingness < 5% per person):** Individuals with sparse genotype data across the genome were excluded to prevent their poor-quality data from distorting population-level statistics.

- **`--hwe 1e-6` (Hardy-Weinberg Equilibrium, p > 1×10⁻⁶):** Severe deviation from Hardy-Weinberg proportions at a SNP typically indicates genotyping error rather than real biology. A stringent threshold of p < 1*10^-6 was used to avoid accidentally removing SNPs showing modest HWE deviation due to population stratification (Wahlund effect), which is expected in a globally diverse dataset.

**Filtering results summary:**
 
| Filter | Variants Removed | Cumulative Remaining |
|--------|-----------------|----------------------|
| Starting variants (chr22) | — | 1,103,547 |
| `--biallelic-only` | 6,343 | 1,097,204 |
| `--mind` | 0 | 1,097,204 |
| `--geno` | 0 | 1,097,204 |
| `--hwe 1e-6` | 47,339 | 1,049,865 |
| `--maf 0.05` | 981,510 | 68,355 |
| **Final QC dataset** | — | **68,355 SNPs, 2,504 individuals** |

**Notable observations:**

 The `--mind` and `--geno` filters removed zero variants or individuals, consistent with the exceptionally complete genotyping in the 1000 Genomes Phase 3 release (total genotyping rate = 1.0). The large number of HWE violations (47,339) likely reflects the **Wahlund effect** - when individuals from genetically differentiated populations are pooled, apparent HWE deviations arise at SNPs where allele frequencies differ between groups, even if each population individually is in HWE. The dominant filter was MAF, which removed ~89% of variants. This reflects the well-established genomic observation that the vast majority of human genetic variants are rare, with common variants (MAF > 5%) representing a minority of total SNP diversity.

 ### 3.4 Step 3 - Linkage Disequilibrium Pruning

Nearby SNPs on the same chromosome tend to be statistically correlated due to linkage disequilibrium (LD), which is the tendency for alleles at physically close loci to be inherited together rather than shuffling independently at each generation. This correlation is a direct consequence of the recombination pattern. Long haplotype blocks are passed from parent to child intact, meaning SNPs within a block carry a lot of redundant information.

If correlated SNPs are included in PCA or ADMIXTURE without pruning, genomic regions with strong LD effectively receive an inflated weight in the analysis, where a dense LD block of "N" correlated SNPs contributes "N" times as much signal as a single SNP elsewhere. This can potentially drive the principal components that reflect local LD architecture rather than genome-wide ancestry patterns.  

LD pruning was performed in 3 steps. First, because all 68,355 post-QC variants carried missing IDs (`.`) in the original VCF (this is a common feature of 1000 Genomes data where not all variants have assigned rsIDs), unique identifiers were assigned using PLINK's `--set-missing-var-ids` flag with the format `chromosome:position:ref_allele:alt_allele`:

```bash
plink \
  --bfile chr22_qc \
  --set-missing-var-ids @:#:$1:$2 \
  --make-bed \
  --out chr22_qc_ids
```

The `chromosome:position:ref:alt` format was required rather than `chromosome:position` because multiple variants existed at identical genomic positions (for example: a SNP and an indel are at the same position of `22:17908596`), which would have produced duplicate IDs under the simpler format.

LD pruning was then performed using PLINK's sliding window algorithm:

```bash
plink \
  --bfile chr22_qc_ids \
  --indep-pairwise 50 10 0.2 \
  --out pruning2
```

The 3 parameters define the pruning action. A window of **50 SNPs** is examined at a time. Within each window, all pairwise r^2 values are calculated, and any SNP forming a pair with r^2 > **0.2** with another SNP in the window is flagged for removal (the SNP having the lower minor allele frequency - MAF is removed preferentially). The window then advances by **10 SNPs** and the process repeats across the chromosome. An r^2 threshold of 0.2 represents the standard for population structure analysis. It is stringent enough to break up meaningful LD blocks while retaining sufficient SNP density for reliable inference. 

The pruned SNP list was then used to extract the independent subset:

```bash
plink \
  --bfile chr22_qc_ids \
  --extract pruning2.prune.in \
  --make-bed \
  --out chr22_pruned
```

**LD pruning results:**
 
| Stage | SNPs |
|-------|------|
| Post-QC input | 68,355 |
| Removed by LD pruning | 60,604 |
| **Retained for analysis** | **7,751** |

The relatively aggressive reduction (89% of SNPs removed) is consistent with the high LD structure expected on chromosome 22, which contains several large LD blocks. Full genome analysis across all 22 autosomes would yield a proportionally larger retained set (around 80,000 - 150,000 independent SNPs), providing greater precision for both PCA and sNMF ancestry estimation. 

**Note on duplicate pruned datasets:** Two separate LD-pruned datasets were created from the same QC-filtered input. The first (`chr22_pruned`) uses colon-separated SNP IDs (`chr:pos:ref:alt`) and was used for PCA, which has no restrictions on variant ID format. The second (`chr22_pruned_admix`) uses underscore-separated IDs (`chr_pos_ref_alt`) and was created specifically for sNMF ancestry estimation. ADMIXTURE v1.3.0 crashes silently when SNP IDs contain colons due to internal parsing behavior, and this issue was discovered during analysis (see Limitations). Both datasets contain identical SNP sets (7,751 variants, 2,504 individuals) - only the ID formatting differs.

### 3.5 Step 4a - Principal Component Analysis

PCA was performed on the LD-pruned dataset using PLINK:

```bash
plink \
  --bfile results/plink/chr22_pruned \
  --pca 10 \
  --out results/pca/chr22_pca
```

PLINK computes PCA by first constructing a genetic relatedness matrix (GRM) - a 2,504 * 2,504 matrix where each cell describes the genotypic covariance between two individuals across all 7,751 pruned SNPs. Eigendecomposition of this matrix yields principal components ranked by the amount of genetic variance they explain. The top 10 principal components were retained. Concurrent execution (27 threads) was used automatically by PLINK for the matrix operations. 

Output files:
- `results/pca/chr22_pca.eigenval` — 10 eigenvalues
- `results/pca/chr22_pca.eigenvec` — PC scores for all 2,504 individuals

### 3.6 Step 4b - Ancestry Estimation via sNMF

**Software note:** Ancestry estimation was initially planned using ADMIXTURE v1.3.0 (Alexander et al, 2009). However, this binary produces a segmentation fault (exit code 139) on WSL2 immediately after reading genotype data, regardless of stack size settings (`ulimit -s unlimited`). This is a known incompatibility between ADMIXTURE's statically linked binary and WSL2's memory model. Analysis was instead performed using the `snmf()` function from the R package LEA (Frichot & François, 2015), which implements sparse Non-negative Matrix Factorization. It is an algorithm that is mathematically equivalent to ADMIXTURE and produces comparable ancestry proportion matrices (Q matrices) and model fit statistics.

Ancestry estimation was run for K = 2 through K = 8 ancestral populations, with 3 independent repetitions per K value to ensure the best solution is found rather than a suboptimal one due to chance starting conditions. The best run per K was selected based on minimum cross-entropy, which is analogous to ADMIXTURE's cross-validation error. It measures the model's accuracy by checking how well the model predicts hidden DNA data. Lower scores mean a better fit.

```r
project <- snmf(
  "results/admixture/chr22.geno",
  K = 2:8,
  entropy = TRUE,
  repetitions = 3,
  project = "new",
  seed = 42
)
```

Prior to running sNMF, PLINK genotypes were exported in raw format (`--recodeA`) and converted to LEA's `.geno` format. It is a plain text matrix where rows represent SNPs, columns represent individuals, and values are 0, 1, 2, or 9 (representing missing data). All 7,751 SNPs had complete genotype data so no missing value encoding was required. 

---

## 4. Results

### 4.1 PCA

The top 10 principal components were computed from 7,751 LD-pruned SNPs across 2,504 individuals. Eigenvalues and percent variance explained are summarized below:

| PC | Eigenvalue | % Variance Explained |
|----|-----------|---------------------|
| PC1 | 143.411 | 48.7% |
| PC2 | 71.062 | 24.1% |
| PC3 | 27.363 | 9.3% |
| PC4 | 22.790 | 7.7% |
| PC5 | 5.511 | 1.9% |
| PC6 | 5.265 | 1.8% |
| PC7 | 4.956 | 1.7% |
| PC8 | 4.773 | 1.6% |
| PC9 | 4.702 | 1.6% |
| PC10 | 4.599 | 1.6% |

PC1 and PC2 together explain 72.8% of total genetic variance. This is a pronounced result demonstrating that the vast majority of population-level genetic signal is captured in just two dimensions. A sharp drop in variance is observed between PC4 (7.7%) and PC5 (1.9%), forming a clear scree plot elbow. This indicates that the four leading PCs capture the major axes of human population structure, while PC5 onward reflect finer-scale or residual variation. 

PC1 is expected to represent the African vs. Non-African divergence, the deepest split in genetic diversity, capturing the out-of-Africa bottleneck ~ 60,000 - 70,000 years ago. PC2 likely captures the East Asian vs. European divergence. PC3 and PC4 are expected to reflect finer substructure within subpopulations. Visual confirmation of these interpretations will be provided by the PCA scatterplot (Section 4.3).

Visual confirmation is provided in Figure 1 (Section 4.3). PC1 separates African (AFR) samples from all non-African populations, and PC2 separates East Asian (EAS) from European (EUR) samples, with South Asian (SAS) samples in an intermediate position. These interpretations are consistent with the eigenvalue distribution and known human demographic history. Admixed American (AMR) samples do not form a discrete cluster but instead scatter between the non-African groups, reflecting their mixed Indigenous American, European, and African ancestries in variable proportions across individuals. These interpretations are also consistent with the eigenvalue distribution and known human demographic history.

### 4.2 Ancestry Estimation (sNMF)

Cross-entropy values for K = 2 through K = 8 are summarized below:

| K | Cross-Entropy | Change from K-1 |
|---|--------------|-----------------|
| 2 | 0.69549 | — |
| 3 | 0.68053 | -0.01496 |
| 4 | 0.67573 | -0.00480 |
| 5 | 0.67164 | -0.00409 |
| 6 | 0.67177 | +0.00013 |
| 7 | 0.67151 | -0.00026 |
| 8 | 0.67138 | -0.00013 |

K = 2 is the minimum tested value, as K = 1 (a single undifferentiated ancestral population) is biologically uninformative and serves no meaningful baseline for comparison. 

**Best-fit K = 5** Cross-entropy decreases meaningfully from K = 2 to K = 5 (total reduction of 0.02385). At K = 6 the value increases slightly before plateauing through K = 7 and K = 8, where improvements become negligible ( < 0.0003 per step), indicating the model is overfitting beyond K = 5. K = 5 corresponds directly to the 5 superpopulations in 1000 Genomes dataset (AFR, EUR, EAS, SAS, AMR), consistent with the known demographic history of globally sampled human populations. Full visualization and biological interpretation of Q matrices are found in Step 4.3. 

### 4.3 Figures

All figures were generated using Python (matplotlib, seaborn) from the PCA and sNMF output files. Scripts are available in `scripts/06_visualize.py`. 

---

**Figure 1 - PCA Scatter Plot (PC1 vs PC2)**

*File: `results/plots/pca_pc1_pc2.png`*

PC1 (48.7% variance explained) cleanly separates African (AFR) samples from all non-African populations. This is the deepest split in human genetic diversity, revealing the genetic bottleneck caused by the early human migration out of Africa approximately 60,000 - 70,000 years ago. When a small founding population left Africa, they carried only a subset of African genetic diversity. It dominates nearly half of all variance in the dataset, making it a profound event. 

PC2 (24.1% variance explained) separates the East Asian (EAS) samples from European (EUR) samples within the non-African cluster. South Asian (SAS) samples occupy an intermediate position between EUR and EAS, consistent with their geographic location and the dual ancestry of South Asian populations from ancient West Eurasian and distinct South Asian lineages. 

African samples do not form a tight cluster. They trail broadly across PC1. This demonstrates greater within-Africa genetic diversity, a consequence of longer time depth and the absence of the major genetic bottleneck that shaped non-African groups. All non-African populations descend from a small ancestral group that left Africa, and their reduced diversity is visible as tighter clustering on the plot. 

Admixed American (AMR) samples scatter between the non-African clusters rather than forming a discrete group, exhibiting their mixed Indigenous American, European, and African ancestries in differing proportions across individuals. Some AMR individuals are located near the EUR cluster, some near AFR, and others in intermediate space. 

---

**Figure 2 - PCA Scree Plot**

*File: `results/plots/pca_scree.png`*

PC1 explains 48.7% of genetic variance, followed by PC2 at 24.1%. Together, they capture 72.8% of the total variance. This means that the majority of population-level genetic signal is essentially compressed into 2 dimensions. PC3 (9.3%) and PC4 (7.7%) likely capture finer sub-structures within superpopulations. A clear plateau begins at PC5 (1.9%), where variance explained drops sharply and remains flat through PC10 (~ 1.6%). This plateau defines the noise floor where PCs beyond this point show sampling variations and linkage disequilibrium (LD) residuals rather than meaningful population structure. The scree plot justifies focusing on PC1 and PC2 as the primary axes of analysis. 

---

**Figure 3 — sNMF Ancestry Bar Chart (K = 5)**

*File: `results/plots/admixture_K5.png`*

At K = 5, the 5 ancestral components map cleanly onto the 5 superpopulations. Ancestral component colors were assigned using the Hungarian algorithm for the best one-to-one mapping between components and superpopulations, ensuring each group receives a unique color regardless of component numbering by the algorithm. 

**AFR (red):** African samples show predominantly African ancestry with visible within-group variation in bar heights, showing greater genetic diversity within Africa compared to other continents. 

**EAS (green):** East Asian samples show the cleanest, most uniform ancestry signal of all 5 groups. The bars are nearly solid green with minimal mixing. This is consistent with the relative geographic isolation of East Asian populations and reduced historical mixing with distant populations. 

**EUR (blue):** European samples show a dominant European ancestry component with some variation. The relatively clean signal shows the post-genetic bottleneck homogeneity of European populations compared to Africans. 

**SAS (purple):** South Asian samples show mixed bars combining a dominant South Asian component (purple) with secondary contributions from other components. This visible mixing reflects the known dual ancestry of South Asian populations - a mixture of West Eurasian (related to EUR) and a distinct ancient South Indian lineage. 

**AMR (orange dominant but highly variable):** Admixed American samples are the most noticeable group. Individual bars  vary substantially, with some individuals showing predominantly orange (Indigenous American) ancestry, others showing substantial red (African), and others showing substantial blue (European). This variability accurately exhibits the colonial history of the Americas. Different individuals within the AMR superpopulation descend from populations with different proportions of Indigenous American, European, and African ancestry, depending on their country and community of origin. The sNMF result shows the genomic record of that history.

---

## 5. Limitations

- This analysis currently uses chromosome 22 only. While sufficient for demonstrating the pipeline and validating methodology, full-genome analysis across all 22 autosomes would provide more precise population structure estimates and greater statistical power for ADMIXTURE.
- The HWE filter was applied to the pooled multi-population sample. A more rigorous approach would apply HWE filters within each population separately, avoiding Wahlund effect false positives while catching true genotyping errors.
- Sex chromosomes (X, Y) were excluded. Population structure on the X chromosome can reveal additional signals of sex-biased migration and demographic history.
- PCA eigenvalues reflect chr22 only and are therefore not directly comparable to published whole-genome results. The high variance explained by PC1 (48.7%) and PC2 (24.1%) is inflated relative to full-genome estimates, where PC1  typically explains 10-20% of variance, because fewer total PCs are competing to explain the same population signal across a single chromosome.
- ADMIXTURE v1.3.0 could not be run on this system due to a segmentation fault (exit code 139) specific to WSL2. The LEA package's `snmf()` function was used as a direct substitute. While the algorithms are mathematically equivalent, results may differ slightly from published ADMIXTURE analyses due to differences in optimization implementation. Future work on a native Linux system or HPC cluster should validate results using ADMIXTURE directly. 

---

## 6. References
 
1. 1000 Genomes Project Consortium (2015). A global reference for human genetic variation. *Nature*, 526, 68–74.
2. Patterson, N., Price, A. L., & Reich, D. (2006). Population structure and eigenanalysis. *PLOS Genetics*, 2(12), e190.
3. Alexander, D. H., Novembre, J., & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655–1664.
4. Frichot, E., & François, O. (2015). LEA: An R package for landscape and ecological association studies. *Methods in Ecology and Evolution*, 6(8), 925–929.
5. Chang, C. C., et al. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, 4(1).
