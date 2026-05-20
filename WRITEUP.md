# Population Structure Analysis - Writeup

**Project:** Population Structure Analysis of the 1000 Genomes Phase 3 Dataset
**Author:** Nicole Linares
**Date:** May 2026
**Repository:** https://github.com/nlinares-GWU2026/1kg-pop-structure

---

## 1. Motivation

Understanding population structure - the non-random distribution of genetic variation across human groups - is crucial to developing and understanding modern population genetics. As humans migrated out of Africa roughly 60,000-70,000 years ago and spread across the globe, geographically isolated populations accumulated distinct patterns of genetic variation through genetic drift, natural selection, and unique mutation histories. These differences are subtle at the individual level but statistically detectable across the genome. 

Population structure and analysis have direct practical importance beyond historical inference. In genome-wide association studies (GWAS), failure to account for population stratification leads to misleading associations between genetic variants and traits, like false positives driven by correlation between ancestry and disease prevalence rather than biology. The methods developed in this project (PCA and ADMIXTURE) are the standard tools used to detect and correct for population stratification in large-scale genetic studies. 

This project uses genome-wide SNP data from the 1000 Genomes Project Phase 3 release to reconstruct population structure across 26 globally distributed human populations, with the goals of: 1) visualizing ancestry-driven genetic clustering via PCA, and 2) estimating ancestral population proportions via ADMIXTURE analysis.

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

This analysis begins with chromosome 22 as a proof-of-concept. Chromosome 22 is the smallest human autosome and provides a computationally traceable starting point before scaling to the full genome (chromosomes 1-22). 

### 2.3 File Format Notes

VCF (Variant Call Format) files store genotype data as a matrix of variants * individuals. Each row represents one genomic position and contains the reference allele, alternate allele(s), and genotype calls for all individuals coded as `0/0` (homozygous reference), `0/1` (heterozygous), or `1/1` (homozygous alternate). The `.tbi` tabix index allows software to jump to any genomic coordinate without reading the entire file. 

---

## 3. Methods

### 3.1 Software and Environment

All analyses were conducted in a conda virtual environment (`popgen`) running on Ubuntu 24 (via Windows Subsystem for Linux). Key software:

| Tool | Version | Purpose |
|------|---------|---------|
| PLINK | v1.9.0-b.8 (Oct 2024) | QC filtering, LD pruning, PCA |
| ADMIXTURE | — | Ancestry estimation |
| Python | 3.11 | Visualization |
| bcftools | — | VCF inspection |

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

- **`--hwe 1e-6` (Hardy-Weinberg Equilibrium, p > 1×10⁻⁶):** Severe deviation from Hardy-Weinberg proportions at SNP typically indicates genotyping error rather than real biology. A stringent threshold of p < 1*10^-6 was used to avoid accidentally removing SNPs showing modest HWE deviation due to population stratification (Wahlund effect), which is expected in a globally diverse dataset.

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

The `chromosome:position:ref:alt` format was required rather than `chromosome:position` because multiple variants existed at identical genomic positions (for example: a SNP and an indel are at the same position of `NN:XXXXXXXX`), which would have produced duplicate IDs under the simpler format.

LD pruning was then performed using PLINK's sliding window algorithm:

```bash
plink \
  --bfile chr22_qc_ids \
  --indep-pairwise 50 10 0.2 \
  --out pruning2
```

The 3 parameters define the pruning action. A window of **50 SNPs** is examined at a time. Within each window, all pairwise r^2 values are calculated, and any SNP forming a pair with r^2 > **0.2** with another SNP in the window is flagged for removal (the SNP having the lower minor allele frequency - MAF is removed preferentially). The window then advances by **10 SNPs** and the process repeats across the chromosome. An r^2 threshold of 0.2 represents the standard for population structure analysis. It is stringent enough to break up meaningful LD blocks while retaining sufficient SNP density for reliable inference. 

The pruned SNP list was then used to extract the independent subest:

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

The relatively aggressive reduction (89% of SNPs removed) is consistent with the high LD structure expected on chromosome 22, which contains several large LD blocks. Full genome analysis across all 22 autosomes would yield a proportionally larger retained set (around 80,000-150,000 independent SNPs), providing greater precision for both PCA and ADMIXTURE. 

---

## 4. Results

*To be completed after Steps 3–5 (LD pruning, PCA, ADMIXTURE, and visualization).*

---

## 5. Limitations

- This analysis currently uses chromosome 22 only. While sufficient for demonstrating the pipeline and validating methodology, full-genome analysis across all 22 autosomes would provide more precise population structure estimates and greater statistical power for ADMIXTURE.
- The HWE filter was applied to the pooled multi-population sample. A more rigorous approach would apply HWE filters within each population separately, avoiding Wahlund effect false positives while catching true genotyping errors.
- Sex chromosomes (X, Y) were excluded. Population structure on the X chromosome can reveal additional signals of sex-biased migration and demographic history.

---

## 6. References
 
1. 1000 Genomes Project Consortium (2015). A global reference for human genetic variation. *Nature*, 526, 68–74.
2. Patterson, N., Price, A. L., & Reich, D. (2006). Population structure and eigenanalysis. *PLOS Genetics*, 2(12), e190.
3. Alexander, D. H., Novembre, J., & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655–1664.
4. Chang, C. C., et al. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, 4(1).
