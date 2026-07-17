#!/usr/bin/env Rscript
# Purpose
#   Estimate realized per-SNP PVE and genome-wide power from benchmark outputs.
# Manuscript section or figure
#   Realized per-SNP PVE and power
# Required software
#   R with data.table and ggplot2
# Required inputs
#   GWAS results, causal variant lists, phenotype simulation outputs
# Outputs
#   Power and realized PVE summary tables
# Example command
#   Rscript scripts/05_effect_size_power/power_analysis.R --project-root . --output-dir derived_data/method_summary
# Note
#   Proprietary genotype and phenotype inputs are not distributed with this repository.

#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# 0) Working directory + output
# ============================================================
setwd(".")
base_dir <- getwd()

out_dir <- file.path(base_dir, "benchmark_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

K_top <- 100
trim_q <- 0.999
save_plots <- TRUE

method_order <- c(
  "BOLT",
  "fastGWA",
  "REGENIE",
  "REGENIE_woLOCO",
  "SLEMM",
  "BLINK_Kansas",
  "BLINK_NYC",
  "FarmCPU_Kansas",
  "FarmCPU_NYC"
)

method_colors <- c(
  "BOLT"            = "#D55E00",
  "REGENIE"         = "#009E73",
  "fastGWA"         = "#0072B2",
  "BLINK_Kansas"    = "#E69F00",
  "BLINK_NYC"       = "#F0E442",
  "FarmCPU_Kansas"  = "#000000",
  "FarmCPU_NYC"     = "#999999",
  "REGENIE_woLOCO"  = "#56B4E9",
  "SLEMM"           = "red"
)

method_linetypes <- c(
  "BOLT"            = "solid",
  "REGENIE"         = "solid",
  "fastGWA"         = "solid",
  "BLINK_Kansas"    = "longdash",
  "BLINK_NYC"       = "solid",
  "FarmCPU_Kansas"  = "dashed",
  "FarmCPU_NYC"     = "solid",
  "REGENIE_woLOCO"  = "dotted",
  "SLEMM"           = "twodash"
)

# ============================================================
# 1) Helper functions
# ============================================================
extract_pheno <- function(x) as.integer(sub(".*pheno([0-9]+).*", "\\1", x))

extract_loc <- function(x) {
  if (grepl("\\(Kansas\\)", x)) return("Kansas")
  if (grepl("\\(NYC\\)", x))    return("NYC")
  NA_character_
}

p_to_chisq <- function(p, eps = 1e-300) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- NA_real_
  p[p <= 0] <- eps
  p[p > 1] <- NA_real_
  qchisq(1 - p, df = 1)
}

lambda_gc <- function(chisq) {
  x <- chisq[is.finite(chisq)]
  if (length(x) == 0) return(NA_real_)
  median(x, na.rm = TRUE) / qchisq(0.5, df = 1)
}

theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "right"
    )
}

save_plot <- function(p, filename, w = 12, h = 6) {
  if (isTRUE(save_plots)) {
    ggsave(file.path(out_dir, filename), plot = p, width = w, height = h, dpi = 300)
  }
}

# Build central 2 Mb null regions from gaps between adjacent effect windows
make_null_center_regions <- function(region_df, center_size = 2e6) {
  x <- copy(region_df)
  setorder(x, CHR, Start, End)
  
  x[, next_start := shift(Start, type = "lead"), by = CHR]
  
  null_df <- x[
    !is.na(next_start),
    .(
      CHR,
      gap_start = End + 1L,
      gap_end   = next_start - 1L
    ),
    by = CHR
  ]
  
  null_df <- null_df[gap_end >= gap_start]
  null_df[, gap_width := gap_end - gap_start + 1L]
  
  # keep only gaps that can contain a central 2 Mb segment
  null_df <- null_df[gap_width >= center_size]
  
  null_df[, mid := floor((gap_start + gap_end) / 2)]
  null_df[, Start := as.integer(mid - center_size / 2 + 1)]
  null_df[, End   := as.integer(mid + center_size / 2)]
  
  # clamp to gap boundaries
  null_df[Start < gap_start, Start := gap_start]
  null_df[End > gap_end, End := gap_end]
  
  null_df[, .(CHR, Start, End)]
}

# Safe region labeling
label_regions <- function(dt, effect_df, null_center_df) {
  dt <- copy(dt)
  setDT(dt)
  setDT(effect_df)
  setDT(null_center_df)
  
  dt[, `:=`(
    CHR = as.integer(CHR),
    BP  = as.integer(BP),
    Region = "Other_Null"
  )]
  
  effect_df[, `:=`(
    CHR   = as.integer(CHR),
    Start = as.integer(Start),
    End   = as.integer(End)
  )]
  
  null_center_df[, `:=`(
    CHR   = as.integer(CHR),
    Start = as.integer(Start),
    End   = as.integer(End)
  )]
  
  dt[
    effect_df,
    on = .(CHR, BP >= Start, BP <= End),
    Region := "Effect_Window"
  ]
  
  idx_null <- dt[Region == "Other_Null", which = TRUE]
  dt_sub <- dt[idx_null]
  
  dt_sub[
    null_center_df,
    on = .(CHR, BP >= Start, BP <= End),
    Region := "Null_Center_2Mb"
  ]
  
  dt[idx_null, Region := dt_sub$Region]
  dt
}

# ============================================================
# 2) Regions
# ============================================================
region_df <- fread(file.path(base_dir, "eligible_regions.csv"))
setDT(region_df)
setnames(region_df, c("chr", "start", "end"), c("CHR", "Start", "End"))
region_df[, `:=`(
  CHR   = as.integer(CHR),
  Start = as.integer(Start),
  End   = as.integer(End)
)]
setorder(region_df, CHR, Start, End)

null_center_df <- make_null_center_regions(region_df, center_size = 2e6)

cat("\nFirst effect windows:\n")
print(region_df[1:min(10, nrow(region_df))])

cat("\nFirst central null regions:\n")
print(null_center_df[1:min(10, nrow(null_center_df))])

# ============================================================
# 3) Readers
# Standardize to:
# SNP, CHR, BP, P, BETA, SE, Method, Location, Method2, Pheno, chisq, NCP
# ============================================================

read_gapit <- function(file, method_name) {
  dt <- fread(file)
  setDT(dt)
  
  effect_col <- if ("effect" %in% names(dt)) {
    "effect"
  } else if ("Effect" %in% names(dt)) {
    "Effect"
  } else {
    NA_character_
  }
  if (is.na(effect_col)) stop("GAPIT missing effect/Effect column: ", basename(file))
  
  loc <- extract_loc(basename(file))
  pheno <- extract_pheno(basename(file))
  if (is.na(pheno)) return(NULL)
  
  dt2 <- dt[, .(
    SNP  = as.character(SNP),
    CHR  = as.integer(Chr),
    BP   = as.integer(Pos),
    P    = as.numeric(P.value),
    BETA = as.numeric(get(effect_col)),
    SE   = NA_real_
  )]
  
  dt2[, `:=`(
    Method   = method_name,
    Location = loc,
    Method2  = paste0(method_name, "_", loc),
    Pheno    = pheno,
    chisq    = p_to_chisq(P),
    NCP      = p_to_chisq(P)
  )]
  
  dt2
}

read_bolt <- function(file) {
  dt <- fread(file)
  setDT(dt)
  
  pheno <- extract_pheno(basename(file))
  if (is.na(pheno)) return(NULL)
  
  dt2 <- dt[, .(
    SNP  = as.character(SNP),
    CHR  = as.integer(CHR),
    BP   = as.integer(BP),
    P    = as.numeric(P_BOLT_LMM_INF),
    BETA = as.numeric(BETA),
    SE   = as.numeric(SE)
  )]
  
  dt2[, `:=`(
    Method   = "BOLT",
    Location = NA_character_,
    Method2  = "BOLT",
    Pheno    = pheno,
    chisq    = (BETA / SE)^2,
    NCP      = (BETA / SE)^2
  )]
  
  dt2
}

read_fastgwa <- function(file) {
  dt <- fread(file)
  setDT(dt)
  
  pheno <- extract_pheno(basename(file))
  if (is.na(pheno)) return(NULL)
  
  dt2 <- dt[, .(
    SNP  = as.character(SNP),
    CHR  = as.integer(CHR),
    BP   = as.integer(POS),
    P    = as.numeric(P),
    BETA = as.numeric(BETA),
    SE   = as.numeric(SE)
  )]
  
  dt2[, `:=`(
    Method   = "fastGWA",
    Location = NA_character_,
    Method2  = "fastGWA",
    Pheno    = pheno,
    chisq    = (BETA / SE)^2,
    NCP      = (BETA / SE)^2
  )]
  
  dt2
}

read_regenie <- function(file, method_name) {
  dt <- fread(file)
  setDT(dt)
  
  pheno <- extract_pheno(basename(file))
  if (is.na(pheno)) return(NULL)
  
  dt2 <- dt[, .(
    SNP  = as.character(ID),
    CHR  = as.integer(CHROM),
    BP   = as.integer(GENPOS),
    P    = 10^(-as.numeric(LOG10P)),
    BETA = as.numeric(BETA),
    SE   = as.numeric(SE)
  )]
  
  dt2[, `:=`(
    Method   = method_name,
    Location = NA_character_,
    Method2  = method_name,
    Pheno    = pheno,
    chisq    = (BETA / SE)^2,
    NCP      = (BETA / SE)^2
  )]
  
  dt2
}

read_slemm <- function(file) {
  dt <- fread(file)
  setDT(dt)
  
  pheno <- extract_pheno(basename(file))
  if (is.na(pheno)) return(NULL)
  
  dt2 <- dt[, .(
    SNP  = as.character(ID),
    CHR  = as.integer(get("#CHROM")),
    BP   = as.integer(POS),
    P    = as.numeric(P),
    BETA = as.numeric(BETA),
    SE   = as.numeric(SE)
  )]
  
  dt2[, `:=`(
    Method   = "SLEMM",
    Location = NA_character_,
    Method2  = "SLEMM",
    Pheno    = pheno,
    chisq    = (BETA / SE)^2,
    NCP      = (BETA / SE)^2
  )]
  
  dt2
}

# ============================================================
# 4) Collect GWAS files
# ============================================================
all_list <- list()

blink_files <- list.files(
  file.path(base_dir, "BLINK"),
  pattern = "GAPIT\\.Association\\.GWAS_Results\\.BLINK\\.pheno.*\\.csv$",
  full.names = TRUE
)
for (f in blink_files) all_list[[length(all_list) + 1]] <- read_gapit(f, "BLINK")

farm_files <- list.files(
  file.path(base_dir, "FarmCPU"),
  pattern = "GAPIT\\.Association\\.GWAS_Results\\.FarmCPU\\.pheno.*\\.csv$",
  full.names = TRUE
)
for (f in farm_files) all_list[[length(all_list) + 1]] <- read_gapit(f, "FarmCPU")

bolt_files <- list.files(
  file.path(base_dir, "BOLTLMM"),
  pattern = "^bolt_pheno[0-9]+\\.stats$",
  full.names = TRUE
)
for (f in bolt_files) all_list[[length(all_list) + 1]] <- read_bolt(f)

fast_files <- list.files(
  file.path(base_dir, "fastGWA"),
  pattern = "pheno[0-9]+\\.fastgwa\\.fastGWA$",
  full.names = TRUE
)
for (f in fast_files) all_list[[length(all_list) + 1]] <- read_fastgwa(f)

reg_files <- list.files(
  file.path(base_dir, "REGENIE"),
  pattern = "^regenie_pheno[0-9]+_pheno[0-9]+\\.regenie$",
  full.names = TRUE
)
for (f in reg_files) all_list[[length(all_list) + 1]] <- read_regenie(f, "REGENIE")

regwl_files <- list.files(
  file.path(base_dir, "REGENIE_WO_LOCO"),
  pattern = "^regenie_pheno_wo_loco[0-9]+_pheno[0-9]+\\.regenie$",
  full.names = TRUE
)
for (f in regwl_files) all_list[[length(all_list) + 1]] <- read_regenie(f, "REGENIE_woLOCO")

slemm_files <- list.files(
  file.path(base_dir, "SLEMM"),
  pattern = "spg\\.even\\.100k\\.pheno[0-9]+\\.allchr\\.gwas\\.txt$",
  full.names = TRUE
)
for (f in slemm_files) all_list[[length(all_list) + 1]] <- read_slemm(f)

cat("\nFiles found:\n")
print(c(
  BLINK = length(blink_files),
  FarmCPU = length(farm_files),
  BOLTLMM = length(bolt_files),
  fastGWA = length(fast_files),
  REGENIE = length(reg_files),
  REGENIE_WO_LOCO = length(regwl_files),
  SLEMM = length(slemm_files)
))

all_list <- Filter(Negate(is.null), all_list)
gwas_all <- rbindlist(all_list, use.names = TRUE, fill = TRUE)

gwas_all <- gwas_all[!is.na(Pheno)]
gwas_all <- label_regions(gwas_all, region_df, null_center_df)
gwas_all[, Method2 := factor(Method2, levels = method_order)]

cat("\nLoaded records by Method2:\n")
print(gwas_all[, .N, by = .(Method2)][order(Method2)])

cat("\nCounts by Region:\n")
print(gwas_all[, .N, by = Region])

# ============================================================
# 4B) Manhattan plots by method and phenotype
#     all SNPs, stacked panels, PNG output
# ============================================================
method_order <- c(
    "BOLT",
    "fastGWA",
    "REGENIE",
    "REGENIE_woLOCO",
    "SLEMM",
    "BLINK_Kansas",
    "BLINK_NYC",
    "FarmCPU_Kansas",
    "FarmCPU_NYC"
)

method_labels <- c(
  "FarmCPU_Kansas"   = "FarmCPU (Kansas)",
  "FarmCPU_NYC"      = "FarmCPU (NYC)",
  "BLINK_Kansas"     = "BLINK (Kansas)",
  "BLINK_NYC"        = "BLINK (NYC)",
  "BOLT"             = "BOLT",
  "REGENIE"          = "REGENIE",
  "REGENIE_woLOCO"   = "REGENIE w/o LOCO",
  "fastGWA"          = "fastGWA",
  "SLEMM"            = "SLEMM"
)

gwas_all[, Method2 := factor(as.character(Method2), levels = method_order)]

sig_threshold <- 5e-8
target_phenos <- c(1, 2, 10)

plot_dt <- copy(gwas_all[Pheno %in% target_phenos])

# clean columns once
plot_dt[, `:=`(
  CHR    = as.integer(CHR),
  BP     = as.numeric(BP),
  P_plot = as.numeric(P)
)]

# guard against invalid P values
plot_dt[!is.finite(P_plot) | is.na(P_plot), P_plot := NA_real_]
plot_dt[P_plot <= 0, P_plot := 1e-300]
plot_dt[P_plot > 1,  P_plot := NA_real_]

# keep usable rows only
plot_dt <- plot_dt[!is.na(P_plot) & !is.na(CHR) & !is.na(BP)]

# re-apply factor after subsetting
plot_dt[, Method2 := factor(as.character(Method2), levels = method_order)]

if (nrow(plot_dt) == 0) {
  cat("\nNo SNPs available for Manhattan plotting.\n")
} else {
  
  # fixed chromosome list for pig autosomes
  all_chr <- data.table(CHR = 1:18)
  
  # chromosome lengths from GWAS data
  chr_len_gwas <- plot_dt[
    ,
    .(chr_len = max(BP, na.rm = TRUE)),
    by = CHR
  ]
  
  # fallback chromosome lengths from region file
  chr_len_region <- region_df[
    ,
    .(chr_len_region = max(as.numeric(End), na.rm = TRUE)),
    by = CHR
  ]
  
  chr_df <- merge(all_chr, chr_len_gwas, by = "CHR", all.x = TRUE)
  chr_df <- merge(chr_df, chr_len_region, by = "CHR", all.x = TRUE)
  
  chr_df[is.na(chr_len), chr_len := chr_len_region]
  chr_df[is.na(chr_len), chr_len := 1]
  
  setorder(chr_df, CHR)
  chr_df[, chr_len := as.numeric(chr_len)]
  chr_df[, tot := cumsum(shift(chr_len, fill = 0))]
  chr_df[, chr_end := tot + chr_len]
  chr_df[, center := tot + chr_len / 2]
  
  genome_max <- max(chr_df$chr_end, na.rm = TRUE)
  
  # merge cumulative offsets after cleaning
  plot_dt <- merge(
    plot_dt,
    chr_df[, .(CHR, tot)],
    by = "CHR",
    all.x = TRUE,
    sort = FALSE
  )
  
  plot_dt <- plot_dt[!is.na(tot)]
  plot_dt[, BPcum := BP + tot]
  plot_dt[, logP := -log10(P_plot)]
  plot_dt[, sig_group := fifelse(P_plot < sig_threshold, "Significant", "Non-significant")]
  
  # re-apply factor after merge
  plot_dt[, Method2 := factor(as.character(Method2), levels = method_order)]
  
  # scaffold to force all chromosomes to appear in each facet
  scaffold <- CJ(
    Pheno = target_phenos,
    Method2 = method_order,
    CHR = 1:18,
    unique = TRUE
  )
  scaffold[, Method2 := factor(Method2, levels = method_order)]
  scaffold <- merge(scaffold, chr_df[, .(CHR, center)], by = "CHR", all.x = TRUE, sort = FALSE)
  scaffold[, `:=`(BPcum = center, logP = 0)]
  
  cat("\nSNP counts by chromosome used for Manhattan plots:\n")
  print(plot_dt[, .N, by = .(Pheno, CHR)][order(Pheno, CHR)])
  
  for (ph in target_phenos) {
    dtm <- plot_dt[Pheno == ph & sig_group == "Significant"]
    scaffold_ph <- scaffold[Pheno == ph]
    
    if (nrow(dtm) == 0) {
      cat("No significant SNPs for pheno", ph, "\n")
      next
    }
    
    ymax <- ceiling(max(dtm$logP, na.rm = TRUE) * 1.05)
    
    p_manhattan <- ggplot() +
      geom_blank(
        data = scaffold_ph,
        aes(x = BPcum, y = logP)
      ) +
      geom_point(
        data = dtm,
        aes(x = BPcum, y = logP),
        color = "grey25",
        size = 0.65,
        alpha = 0.9
      ) +
      scale_x_continuous(
        limits = c(0, genome_max),
        breaks = chr_df$center,
        labels = as.character(chr_df$CHR),
        expand = expansion(mult = c(0.01, 0.04))
      ) +
      scale_y_continuous(
        limits = c(0, ymax),
        expand = expansion(mult = c(0, 0.02))
      ) +
      coord_cartesian(clip = "off") +
      facet_wrap(~ Method2, ncol = 1, scales = "fixed", labeller = as_labeller(method_labels)) +
      labs(
        x = "Chromosome",
        y = expression(-log[10](P)),
        title = paste0("Manhattan plot by GWAS method: pheno", ph)
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2),
        axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 11),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        plot.title = element_text(color = "black", face = "bold"),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(size = 13, face = "bold", color = "black"),
        legend.position = "none",
        plot.margin = margin(10, 30, 10, 10)
      )
    
    ggsave(
      filename = file.path(
        out_dir,
        paste0("Manhattan_significant_SNPs_by_method_pheno", ph, ".png")
      ),
      plot = p_manhattan,
      width = 20,
      height = 18,
      dpi = 300,
      bg = "white"
    )
  }
  
  cat("\nSaved Manhattan plots by method for pheno1, pheno2, and pheno10 to:\n")
  cat(out_dir, "\n")
}

# ============================================================
# 5) ALT SNPs for power analysis
# ============================================================
alt <- fread(file.path(base_dir, "snp.info.chr1_18.10k.csv"))
alt[, SNP_base := sub("_[^_]+$", "", SNP)]
alt[, c("CHR", "BP") := tstrsplit(SNP_base, "_", fixed = TRUE)]
alt[, `:=`(CHR = as.integer(CHR), BP = as.integer(BP))]

alt_key <- unique(alt[, .(CHR, BP)])
setkey(alt_key, CHR, BP)
setkey(gwas_all, CHR, BP)

gwas_all[, is_alt := FALSE]
gwas_all[alt_key, is_alt := TRUE]

# ============================================================
# 6) Power and FPR curves
# ============================================================
alpha_grid <- 10^seq(-10, -2, by = 0.25)

power_curve <- gwas_all[
  is_alt == TRUE,
  {
    pvals <- P
    data.table(
      alpha = alpha_grid,
      power = vapply(alpha_grid, function(a) mean(pvals <= a, na.rm = TRUE), numeric(1))
    )
  },
  by = .(Method2, Pheno)
]

fpr_curve <- gwas_all[
  Region == "Null_Center_2Mb",
  {
    pvals <- P
    data.table(
      alpha = alpha_grid,
      fpr = vapply(alpha_grid, function(a) mean(pvals <= a, na.rm = TRUE), numeric(1))
    )
  },
  by = .(Method2, Pheno)
]

power_curve[, Method2 := factor(Method2, levels = method_order)]
fpr_curve[, Method2 := factor(Method2, levels = method_order)]

p_power <- ggplot(power_curve, aes(x = alpha, y = power, color = Method2, linetype = Method2)) +
  geom_line(linewidth = 1.1) +
  scale_x_log10() +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_linetype_manual(values = method_linetypes, drop = FALSE) +
  facet_wrap(~ Pheno, ncol = 3) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  labs(
    title = "Empirical power curves (ALT SNPs)",
    x = "P-value threshold (alpha, log scale)",
    y = "Power",
    color = "Method",
    linetype = "Method"
  )
print(p_power)
save_plot(p_power, "power_curve_by_pheno.png", w = 14, h = 7)

p_fpr <- ggplot(fpr_curve, aes(x = alpha, y = fpr, color = Method2, linetype = Method2)) +
  geom_line(linewidth = 1.1) +
  scale_x_log10() +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_linetype_manual(values = method_linetypes, drop = FALSE) +
  facet_wrap(~ Pheno, ncol = 3) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  labs(
    title = "Empirical false positive curves (central 2 Mb of null windows)",
    x = "P-value threshold (alpha, log scale)",
    y = "Type I error",
    color = "Method",
    linetype = "Method"
  )
print(p_fpr)
save_plot(p_fpr, "fpr_curve_null_center2Mb_by_pheno.png", w = 14, h = 7)

# ============================================================
# 7) Significant SNP counts in effect windows vs deep-null regions
# ============================================================
alpha_sig <- 5e-8

# keep only the two target regions, then count
sig_summary_raw <- gwas_all[
  Region %in% c("Effect_Window", "Null_Center_2Mb"),
  .(n_sig = sum(P <= alpha_sig, na.rm = TRUE)),
  by = .(Pheno, Method2, Region)
]

# force exact region labels and factor levels
sig_summary_raw[, Region := trimws(as.character(Region))]
sig_summary_raw[, Method2 := as.character(Method2)]

# complete all combinations explicitly
sig_template <- CJ(
  Pheno = sort(unique(gwas_all$Pheno)),
  Method2 = method_order,
  Region = c("Effect_Window", "Null_Center_2Mb"),
  unique = TRUE
)

sig_summary <- merge(
  sig_template,
  sig_summary_raw,
  by = c("Pheno", "Method2", "Region"),
  all.x = TRUE,
  sort = FALSE
)

sig_summary[is.na(n_sig), n_sig := 0L]

# enforce factors after merge
sig_summary[, Method2 := factor(Method2, levels = method_order)]
sig_summary[, Region := factor(Region, levels = c("Effect_Window", "Null_Center_2Mb"))]

# sanity check
cat("\nSignificant SNP summary used for plotting:\n")
print(sig_summary[order(Pheno, Method2, Region)])

sig_summary_wide <- dcast(
  sig_summary,
  Pheno + Method2 ~ Region,
  value.var = "n_sig",
  fill = 0
)

setnames(
  sig_summary_wide,
  c("Effect_Window", "Null_Center_2Mb"),
  c("Sig_Effect", "Sig_DeepNull")
)

sig_summary_wide[, Spillover_Ratio := Sig_DeepNull / (Sig_Effect + 1)]

fwrite(
  sig_summary_wide,
  file.path(base_dir, "significant_snp_counts_by_region.csv")
)

p_sig <- ggplot(
  sig_summary,
  aes(x = Method2, y = n_sig, fill = Region)
) +
  geom_col(
    position = position_dodge2(width = 0.8, preserve = "single"),
    width = 0.7,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = n_sig),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "Effect_Window" = "#F8766D",
      "Null_Center_2Mb" = "#00BFC4"
    ),
    drop = FALSE
  ) +
  facet_wrap(~ Pheno, ncol = 3) +
  theme_pub() +
  labs(
    title = "Significant SNP counts by genomic region",
    x = NULL,
    y = paste0("Number of SNPs (P <= ", alpha_sig, ")"),
    fill = "Region"
  )

print(p_sig)
save_plot(p_sig, "06_significant_snps_effect_vs_deepnull.png", w = 14, h = 7)

# ============================================================
# 8) Inflation summary: effect vs central null
# ============================================================
inflation_summary <- gwas_all[
  Region %in% c("Effect_Window", "Null_Center_2Mb"),
  .(
    n_snp = .N,
    mean_chisq = mean(chisq, na.rm = TRUE),
    median_chisq = median(chisq, na.rm = TRUE),
    lambda_gc = lambda_gc(chisq),
    mean_log10p = mean(-log10(pmax(P, 1e-300)), na.rm = TRUE),
    median_log10p = median(-log10(pmax(P, 1e-300)), na.rm = TRUE)
  ),
  by = .(Pheno, Method2, Region)
][order(Pheno, Method2, Region)]

fwrite(inflation_summary, file.path(base_dir, "inflation_summary_effect_vs_null_center2Mb.csv"))

# ============================================================
# 9) Signal enrichment: effect vs central null
# ============================================================
enrich <- gwas_all[
  Region %in% c("Effect_Window", "Null_Center_2Mb"),
  .(
    mean_chisq = mean(chisq, na.rm = TRUE),
    median_chisq = median(chisq, na.rm = TRUE)
  ),
  by = .(Pheno, Method2, Region)
]

enrich_ratio <- enrich[
  ,
  {
    eff <- mean_chisq[Region == "Effect_Window"]
    nul <- mean_chisq[Region == "Null_Center_2Mb"]
    data.table(EnrichmentRatio = eff / nul)
  },
  by = .(Pheno, Method2)
]

p1 <- ggplot(enrich_ratio, aes(x = Method2, y = EnrichmentRatio)) +
  geom_col(width = 0.75) +
  facet_wrap(~ Pheno, ncol = 3) +
  theme_pub() +
  labs(
    title = "Signal enrichment: mean chi-square in effect windows / central null regions",
    x = NULL,
    y = "Enrichment ratio"
  )
print(p1)
save_plot(p1, "01_enrichment_ratio_by_pheno.png", w = 14, h = 7)

# ============================================================
# 10) Inflation in central null regions
# ============================================================
infl_null <- gwas_all[
  Region == "Null_Center_2Mb",
  .(lambda_GC = lambda_gc(chisq)),
  by = .(Pheno, Method2)
]

p2 <- ggplot(infl_null, aes(x = Method2, y = lambda_GC)) +
  geom_col(width = 0.75) +
  facet_wrap(~ Pheno, ncol = 3) +
  theme_pub() +
  labs(
    title = "Inflation in central 2 Mb null regions",
    x = NULL,
    y = expression(lambda[GC])
  )
print(p2)
save_plot(p2, "02_lambdaGC_null_center2Mb_by_pheno.png", w = 14, h = 7)

# ============================================================
# 11) Chi-square density by region and phenotype
# ============================================================
density_dt <- gwas_all[Region %in% c("Effect_Window", "Null_Center_2Mb")]
xmax <- as.numeric(quantile(density_dt$chisq, trim_q, na.rm = TRUE))

p3 <- ggplot(density_dt, aes(x = chisq, color = Method2)) +
  geom_density(linewidth = 1) +
  facet_grid(Region ~ Pheno) +
  theme_pub() +
  coord_cartesian(xlim = c(0, xmax)) +
  labs(
    title = expression(chi^2 ~ "distribution by region and phenotype"),
    x = expression(chi^2 ~ "(df = 1)"),
    y = "Density",
    color = "Method"
  )
print(p3)
save_plot(p3, "03_chisq_density_region_by_pheno.png", w = 27, h = 12)

# ============================================================
# 12) Mean chi-square in effect vs central null
#    Handle Inf/NA values before plotting
# ============================================================

# Diagnostic: count non-finite chisq values by method
cat("\nNon-finite chisq values by method:\n")
print(
  gwas_all[, .(
    nonfinite = sum(!is.finite(chisq)),
    total = .N
  ), by = Method2][order(-nonfinite)]
)

# Use only finite chisq values for summaries/plots
mean_chisq_plot <- gwas_all[
  Region %in% c("Effect_Window", "Null_Center_2Mb") & is.finite(chisq),
  .(mean_chisq = mean(chisq, na.rm = TRUE)),
  by = .(Pheno, Method2, Region)
]

# Save underlying table
fwrite(
  mean_chisq_plot,
  file.path(base_dir, "mean_chisq_effect_vs_null_by_pheno.csv")
)

# Print min/max by method to inspect truncated bars
cat("\nMean chi-square summary by method:\n")
print(
  mean_chisq_plot[, .(
    min_mean_chisq = min(mean_chisq, na.rm = TRUE),
    max_mean_chisq = max(mean_chisq, na.rm = TRUE),
    median_mean_chisq = median(mean_chisq, na.rm = TRUE)
  ), by = Method2][order(-max_mean_chisq)]
)

# Optional: print the most extreme rows
cat("\nLargest mean chi-square rows:\n")
print(
  mean_chisq_plot[order(-mean_chisq)][1:min(20, .N)]
)

# Safe y-axis upper limit
ylim_max <- max(mean_chisq_plot$mean_chisq, na.rm = TRUE) * 1.10

p4 <- ggplot(mean_chisq_plot, aes(x = Method2, y = mean_chisq, fill = Region)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  geom_text(
    aes(label = round(mean_chisq, 2)),
    position = position_dodge(width = 0.85),
    vjust = -0.25,
    size = 3
  ) +
  facet_wrap(~ Pheno, ncol = 3) +
  coord_cartesian(ylim = c(0, ylim_max)) +
  theme_pub() +
  labs(
    title = "Mean chi-square in effect windows vs central null regions",
    x = NULL,
    y = "Mean chi-square",
    fill = "Region"
  )

print(p4)
save_plot(p4, "04_mean_chisq_effect_vs_null_by_pheno.png", w = 14, h = 7)

# ============================================================
# 13) Top-hit localization
# ============================================================
top1 <- gwas_all[order(P), .SD[1], by = .(Pheno, Method2)]
top1[, TopInEffectWindow := (Region == "Effect_Window")]

p5a <- ggplot(top1, aes(x = Method2, y = factor(Pheno), fill = TopInEffectWindow)) +
  geom_tile() +
  theme_pub() +
  labs(
    title = "Top-1 localization: is the top hit in an effect window?",
    x = NULL,
    y = "Phenotype"
  )
print(p5a)
save_plot(p5a, "05A_top1_in_effect_window_tile.png", w = 14, h = 5)

topk <- gwas_all[order(P), .SD[1:K_top], by = .(Pheno, Method2)]
topk_sum <- topk[, .(
  PropTopKInEffectWindow = mean(Region == "Effect_Window")
), by = .(Pheno, Method2)]

p5b <- ggplot(topk_sum, aes(x = Method2, y = PropTopKInEffectWindow)) +
  geom_col(width = 0.75) +
  facet_wrap(~ Pheno, ncol = 3) +
  theme_pub() +
  labs(
    title = sprintf("Top-%d localization: fraction of top hits in effect windows", K_top),
    x = NULL,
    y = sprintf("Prop(top-%d in effect windows)", K_top)
  )
print(p5b)
save_plot(p5b, sprintf("05B_top%d_in_effect_window_by_pheno.png", K_top), w = 14, h = 7)

# ============================================================
# 14) Save key tables
# ============================================================
fwrite(enrich_ratio, file.path(base_dir, "window_enrichment_ratio_all_methods.csv"))
fwrite(infl_null, file.path(base_dir, "lambdaGC_null_center2Mb_all_methods.csv"))
fwrite(top1, file.path(base_dir, "top1_localization_all_methods.csv"))
fwrite(topk_sum, file.path(base_dir, sprintf("top%d_localization_all_methods.csv", K_top)))
fwrite(power_curve, file.path(base_dir, "power_curve_all_methods.csv"))
fwrite(fpr_curve, file.path(base_dir, "fpr_curve_null_center2Mb_all_methods.csv"))

cat("\nAll done. Plots saved to:\n", out_dir, "\n", sep = "")

alpha_sig <- 5e-8
sig_summary <- gwas_all[
  P <= alpha_sig & Region %in% c("Effect_Window", "Null_Center_2Mb"),
  .(
    n_sig = .N
  ),
  by = .(Pheno, Method2, Region)
]

sig_summary_wide <- dcast(
  sig_summary,
  Pheno + Method2 ~ Region,
  value.var = "n_sig",
  fill = 0
)

setnames(sig_summary_wide,
         c("Effect_Window", "Null_Center_2Mb"),
         c("Sig_Effect", "Sig_DeepNull"))

sig_summary_wide[, Spillover_Ratio := Sig_DeepNull / (Sig_Effect + 1)]

fwrite(
  sig_summary_wide,
  file.path(base_dir, "significant_snp_counts_by_region.csv")
)

p_sig <- ggplot(sig_summary, aes(x = Method2, y = n_sig, fill = Region)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Pheno, ncol = 3) +
  theme_pub() +
  labs(
    title = "Significant SNP counts by genomic region",
    x = NULL,
    y = paste0("Number of SNPs (P ≤ ", alpha_sig, ")")
  )

print(p_sig)

save_plot(p_sig, "significant_snps_effect_vs_deepnull.png", w = 14, h = 7)

