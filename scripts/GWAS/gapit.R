#!/usr/bin/env Rscript

# =============================
# Define paths
# =============================

R_LIB_PATH <- "/path/to/Rlibs"
WORK_DIR <- "./results"
GAIT_FUNCTIONS_FILE <- "/path/to/software/FarmCPU/gapit_functions.txt"
PHENO_FILE <- "/path/to/phenotype/pheno_combined.txt"
GENO_FILE <- "/path/to/genotype/spg.even.100k.chrpos.raw"
GM_FILE <- "/path/to/genotype/spg.even.100k.GM.txt"

if (!dir.exists(R_LIB_PATH)) {
  dir.create(R_LIB_PATH, recursive = TRUE)
}

.libPaths(c(R_LIB_PATH, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!dir.exists(WORK_DIR)) {
  dir.create(WORK_DIR, recursive = TRUE)
}

setwd(WORK_DIR)

# -----------------------------
# Libraries
# -----------------------------
library(data.table)
# Source local FarmCPU code
# -----------------------------
source(GAIT_FUNCTIONS_FILE)

# -----------------------------
# Load data
# -----------------------------

## ---- Phenotype ----
myY <- read.table(
  PHENO_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Remove FID column
myY <- myY[, -1]

# Rename IID to Taxa
colnames(myY)[1] <- "Taxa"

## ---- Genotype (PLINK recode A) ----
GD_raw <- fread(
  GENO_FILE,
  data.table = FALSE
)

## Drop PLINK metadata columns:
## FID PAT MAT SEX PHENOTYPE
## Build GAPIT-safe GD:
## Taxa (IID) column MUST be first
duroc_GD <- data.frame(
  Taxa = GD_raw$IID,
  GD_raw[, -(1:6)],
  check.names = FALSE
)

## ---- Genotype map rebuilt from .bim ----
duroc_GM <- read.table(
  GM_FILE,
  head = TRUE,
)

## ============================================================
## 3. HARD SAFETY CHECKS (STEP 4)
## ============================================================
cat("Checking GD / GM consistency...\n")

cat("GD markers:", ncol(duroc_GD) - 1, "\n")
cat("GM markers:", nrow(duroc_GM), "\n")

## ============================================================
## FIX GM ORDER TO MATCH GD
## ============================================================

cat("Reordering GM to match GD SNP order...\n")

idx <- match(colnames(duroc_GD)[-1], duroc_GM$SNP)

## Hard check: all SNPs must be found
stopifnot(!any(is.na(idx)))

duroc_GM <- duroc_GM[idx, ]

## Final safety check
stopifnot(
  identical(colnames(duroc_GD)[-1], duroc_GM$SNP)
)

cat("GD and GM are now perfectly aligned.\n")

## Exact match required
stopifnot(
  ncol(duroc_GD) - 1 == nrow(duroc_GM),
  identical(colnames(duroc_GD)[-1], duroc_GM$SNP)
)

cat("GD and GM are perfectly aligned.\n")

## ============================================================
## 4. Run GAPIT (STEP 5)
## ============================================================

cat("Starting GAPIT FarmCPU / BLINK...\n")

myGAPIT <- GAPIT(
  Y = myY,
  GD = duroc_GD,
  GM = duroc_GM,
  model = "BLINK",
  Multiple_analysis = TRUE
)

cat("GAPIT finished successfully.\n")
