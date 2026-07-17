#!/usr/bin/env Rscript
# Purpose
#   Summarize significant SNP and 500 kb locus counts across methods.
# Manuscript section or figure
#   Significant SNP and locus summaries, Figure 3 / Figure 4
# Required software
#   R with data.table, ggplot2, dplyr, and stringr
# Required inputs
#   Method-specific GWAS outputs
# Outputs
#   Summary tables and optional comparison figures
# Example command
#   Rscript scripts/06_significant_loci/summarize_significant_snps_and_loci.R --project-root . --output-dir derived_data/method_summary
# Note
#   Proprietary genotype and phenotype inputs are not distributed with this repository.

setwd(".")

threshold <- 5e-8
log_threshold <- -log10(threshold)

methods <- c("SLEMM", "BOLT_lmmInfOnly", "REGENIE", "fastGWA", "FarmCPU", "BLINK")

results <- data.frame(Method=character(),
                      File=character(),
                      Sig_SNPs=integer(),
                      stringsAsFactors=FALSE)

for (m in methods) {
  
  files <- list.files(path=m, recursive=TRUE, full.names=TRUE)
  
  for (f in files) {
    
    sig_count <- NA
    
    # -------- SLEMM --------
    if (grepl("allchr\\.gwas\\.txt$", f)) {
      
      gwas <- read.table(f, header=TRUE, sep="\t",
                         check.names=FALSE, comment.char="")
      gwas$P <- as.numeric(gwas$P)
      sig_count <- sum(gwas$P < threshold, na.rm=TRUE)
    }
    
    # -------- BOLT --------
    if (grepl("\\.stats$", f)) {
      
      gwas <- read.table(f, header=TRUE)
      gwas$P <- as.numeric(gwas$P_BOLT_LMM_INF)
      sig_count <- sum(gwas$P < threshold, na.rm=TRUE)
    }
    
    # -------- REGENIE --------
    if (grepl("\\.regenie$", f)) {
      
      gwas <- read.table(f, header=TRUE)
      gwas$LOG10P <- as.numeric(gwas$LOG10P)
      sig_count <- sum(gwas$LOG10P > log_threshold, na.rm=TRUE)
    }
    
    # -------- fastGWA --------
    if (grepl("\\.fastGWA$", f)) {
      
      gwas <- read.table(f, header=TRUE)
      gwas$P <- as.numeric(gwas$P)
      sig_count <- sum(gwas$P < threshold, na.rm=TRUE)
    }
    
    # -------- FarmCPU / BLINK --------
    if (grepl("\\.csv$", f) &&
        (grepl("FarmCPU", f) || grepl("BLINK", f))) {
      
      gwas <- read.csv(f, header=TRUE)
      gwas$P.value <- as.numeric(gwas$P.value)
      sig_count <- sum(gwas$P.value < threshold, na.rm=TRUE)
    }
    
    if (!is.na(sig_count)) {
      results <- rbind(results,
                       data.frame(Method=m,
                                  File=basename(f),
                                  Sig_SNPs=sig_count))
    }
  }
}

write.csv(results, "GWAS_significant_SNP_summary.csv", row.names=FALSE)

print(results)

library(ggplot2)
library(dplyr)
library(stringr)

df <- read.csv("GWAS_significant_SNP_summary.csv")

# Extract phenotype
df$Phenotype <- str_extract(df$File, "pheno[0-9]+")

# Extract location (Kansas / NYC)
df$Location <- str_extract(df$File, "Kansas|NYC")

# Create new method label
df$Method2 <- df$Method
df$Method2[df$Method %in% c("FarmCPU", "BLINK")] <-
  paste(df$Method[df$Method %in% c("FarmCPU", "BLINK")],
        df$Location[df$Method %in% c("FarmCPU", "BLINK")])

# 🔹 Keep only pheno10
df_sub <- df %>%
  filter(Phenotype == "pheno1")

# Set common y-axis
ymax <- max(df_sub$Sig_SNPs) * 1.15

pdf("pheno1_comparison.pdf", width = 7, height = 5)

ggplot(df_sub,
       aes(x = Method2, y = Sig_SNPs, fill = Method2)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Sig_SNPs),
            vjust = -0.4,
            size = 4) +
  ylim(0, ymax) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(title = "Significant SNPs Comparison (pheno1)",
       y = "Number of Significant SNPs",
       x = "Method")

dev.off()

