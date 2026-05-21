library(LEA)

# --- Step 1: Convert PLINK raw to LEA .geno format ---
cat("Converting PLINK raw to .geno format...\n")

raw <- read.table("results/admixture/chr22.raw",
                  header = TRUE, stringsAsFactors = FALSE)

# Remove the first 6 PLINK metadata columns
# (FID, IID, PAT, MAT, SEX, PHENOTYPE)
geno <- t(as.matrix(raw[, -(1:6)]))

# LEA uses 9 to encode missing data (we have none, but good practice)
geno[is.na(geno)] <- 9

# Write .geno file — no spaces between values, one SNP per row
geno_file <- "results/admixture/chr22.geno"
write(apply(geno, 1, paste, collapse = ""), geno_file)
cat("Written:", geno_file, "\n")
cat("Dimensions: ", nrow(geno), "SNPs x", ncol(geno), "individuals\n\n")

# --- Step 2: Run sNMF for K = 2 through 8 ---
cat("Running sNMF (equivalent to ADMIXTURE)...\n")

project <- snmf(
  geno_file,
  K = 2:8,
  entropy = TRUE,     # cross-entropy = equivalent to CV error
  repetitions = 3,    # run each K 3 times, keep best
  project = "new",
  seed = 42
)

# --- Step 3: Print cross-entropy for each K ---
cat("\n--- Cross-entropy by K (lower = better fit) ---\n")
for (k in 2:8) {
  ce <- cross.entropy(project, K = k)
  cat(sprintf("K=%d: %.5f\n", k, min(ce)))
}

cat("\nDone! Q matrices saved in results/admixture/\n")
