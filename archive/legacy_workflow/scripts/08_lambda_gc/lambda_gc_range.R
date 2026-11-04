#!/usr/bin/env Rscript
# Purpose
#   Summarize genomic inflation in central-null regions.
# Manuscript section or figure
#   Figure 6B
# Required software
#   R with data.table
# Required inputs
#   GWAS results and NULL window definitions
# Outputs
#   lambda GC summary tables
# Example command
#   Rscript scripts/08_lambda_gc/lambda_gc_range.R --project-root . --output-dir derived_data/lambda_gc_summary
# Note
#   Proprietary genotype and phenotype inputs are not distributed with this repository.

#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
})

setwd(".")

# -----------------------------
# Helper functions
# -----------------------------
p_to_chisq <- function(p, eps = 1e-300) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- NA_real_
  p[p <= 0] <- eps
  p[p > 1] <- NA_real_
  qchisq(1 - p, df = 1)
}

lambda_gc <- function(chisq) {
  x <- chisq[is.finite(chisq)]
  median(x, na.rm = TRUE) / qchisq(0.5, df = 1)
}

# -----------------------------
# Load your saved summary table
# -----------------------------
dt <- fread("inflation_summary_effect_vs_null_center2Mb.csv")

# Keep central null regions only
null_dt <- dt[Region == "Null_Center_2Mb"]

# -----------------------------
# λGC by phenotype and method
# -----------------------------
lambda_by_pheno <- null_dt[, .(
  lambda_GC = lambda_gc(rep(median_chisq, 1))
), by = .(Pheno, Method2)]

# Since the file already contains lambda_gc, this is better:
lambda_by_pheno <- null_dt[, .(
  lambda_GC = lambda_gc
), by = .(Pheno, Method2)]

fwrite(lambda_by_pheno, "lambdaGC_by_pheno_method.csv")

# -----------------------------
# Range across phenotypes
# -----------------------------
lambda_range_by_method <- lambda_by_pheno[, .(
  mean_lambda = mean(lambda_GC, na.rm = TRUE),
  min_lambda  = min(lambda_GC, na.rm = TRUE),
  max_lambda  = max(lambda_GC, na.rm = TRUE),
  sd_lambda   = sd(lambda_GC, na.rm = TRUE),
  n_pheno     = sum(is.finite(lambda_GC))
), by = Method2][order(Method2)]

print(lambda_range_by_method)

fwrite(lambda_range_by_method, "lambdaGC_range_by_method.csv")
