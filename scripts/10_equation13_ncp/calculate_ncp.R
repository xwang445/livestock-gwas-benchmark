#!/usr/bin/env Rscript
# Purpose
#   Calculate theoretical and observed non-centrality parameter summaries for Equation 13.
# Manuscript section or figure
#   Equation 13 and Figure 8
# Required software
#   R with data.table and ggplot2
# Required inputs
#   GWAS results, causal LD scores, and effect-size summaries
# Outputs
#   Equation 13 summary tables and diagnostic figures
# Example command
#   Rscript scripts/10_equation13_ncp/calculate_ncp.R --project-root . --output-dir derived_data/equation13_summary
# Note
#   Proprietary genotype and phenotype inputs are not distributed with this repository.

rm(list=ls())

library(data.table)
library(ggplot2)

setwd(".")

plink_dir <- "./"
bolt_dir <- "."

p1_blink_kansas <- fread(
  file.path(plink_dir, "GAPIT.Association.GWAS_Results.BLINK.pheno1(Kansas).csv")
)

p1_blink_kansas$chisq <- qchisq(1 - p1_blink_kansas$P.value, df = 1)
p1_blink_kansas$NCP <- p1_blink_kansas$chisq
plot(p1_blink_kansas$NCP)

p1_bolt <- fread(
  file.path(bolt_dir, "bolt_pheno1.stats")
)
p1_bolt$chisq <- qchisq(1 - p1_bolt$P_BOLT_LMM_INF, df = 1)
p1_bolt$NCP <- p1_bolt$chisq
plot(p1_bolt$NCP)

ggplot(p1_bolt, aes(x = BETA)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(title = "Distribution of BETA (BOLT)",
       x = "Effect size (BETA)",
       y = "Count")

region_df <- fread("eligible_regions.csv")

# Ensure data.table
setDT(p1_bolt)
setDT(region_df)

# Standardize region_df column names to match join keys
setnames(region_df, old = c("chr", "start", "end"), new = c("CHR", "Start", "End"))

# Make sure types are consistent
p1_bolt[, CHR := as.integer(CHR)]
p1_bolt[, BP  := as.integer(BP)]

region_df[, CHR   := as.integer(CHR)]
region_df[, Start := as.integer(Start)]
region_df[, End   := as.integer(End)]

# Default
p1_bolt[, Region := "Outside"]

# Range overlap join: mark SNPs that fall into ANY region
p1_bolt[
  region_df,
  on = .(CHR, BP >= Start, BP <= End),
  Region := "Causal_Window"
]

list.dirs(".")
list.files(".",
           pattern = "spg.even.100k.pheno.*.allchr.gwas.txt")
