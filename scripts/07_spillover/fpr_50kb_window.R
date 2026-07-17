#!/usr/bin/env Rscript
# Purpose
#   Analyze 50 kb null-window false-positive rates for the spillover benchmark.
# Manuscript section or figure
#   Figures 6C and 6D
# Required software
#   R with data.table and ggplot2
# Required inputs
#   GWAS results and null-window definitions
# Outputs
#   Figure files and FPR summary tables
# Example command
#   Rscript scripts/07_spillover/fpr_50kb_window.R --project-root . --output-dir derived_data/spillover_summary
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

out_dir <- file.path(base_dir, "benchmark_plots_50kb_window_fpr")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plots <- TRUE
alpha_sig <- 5e-8
ncol_pheno <- 5

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

# ============================================================
# 1) Helper functions
# ============================================================
extract_pheno <- function(x) as.integer(sub(".*pheno([0-9]+).*", "\\1", x))

extract_loc <- function(x) {
  if (grepl("\\(Kansas\\)", x)) return("Kansas")
  if (grepl("\\(NYC\\)", x))    return("NYC")
  NA_character_
}

save_plot <- function(p, filename, w = 12, h = 6) {
  if (isTRUE(save_plots)) {
    ggsave(
      filename = file.path(out_dir, filename),
      plot = p,
      width = w,
      height = h,
      dpi = 300,
      bg = "white"
    )
  }
}

theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "right"
    )
}

# ============================================================
# 2) Null-window construction
#    eligible_regions.csv already contains the 2 Mb effect windows
#    Each 4 Mb null window is the gap between adjacent rows
# ============================================================
make_null_windows_from_adjacent_effects <- function(region_df) {
  x <- copy(region_df)
  setDT(x)
  setorder(x, CHR, Start, End)
  
  x[, next_start := shift(Start, type = "lead"), by = CHR]
  
  null_df <- x[
    !is.na(next_start),
    .(
      NULL_ID = .I,
      CHR     = as.integer(CHR),
      Start   = as.integer(End),
      End     = as.integer(next_start)
    )
  ]
  
  null_df[, width := End - Start]
  null_df <- null_df[width == 4000000L]
  null_df[, width := NULL]
  
  null_df[]
}

# Split each 4 Mb null window into ordered 50 kb bins
# left boundary -> center : order 1..n
# right boundary -> center: order 1..n
make_ordered_50kb_windows <- function(null_df, win_size = 50000L) {
  x <- copy(null_df)
  setDT(x)
  
  if (nrow(x) == 0L) return(data.table())
  
  out <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
    chr_i   <- x$CHR[i]
    null_id <- x$NULL_ID[i]
    s0      <- x$Start[i]
    e0      <- x$End[i]
    
    null_width <- e0 - s0
    half_width <- null_width / 2
    
    if (null_width != 4000000L) {
      stop("Null window is not exactly 4 Mb for NULL_ID = ", null_id)
    }
    if (half_width %% win_size != 0) {
      stop("Half null window is not divisible by 50 kb for NULL_ID = ", null_id)
    }
    
    n <- half_width / win_size
    
    left_start <- s0 + (0:(n - 1L)) * win_size
    left_end   <- left_start + win_size
    
    left_dt <- data.table(
      NULL_ID   = null_id,
      CHR       = chr_i,
      side      = "left",
      order     = 1:n,
      win_start = as.integer(left_start),
      win_end   = as.integer(left_end)
    )
    
    right_end   <- e0 - (0:(n - 1L)) * win_size
    right_start <- right_end - win_size
    
    right_dt <- data.table(
      NULL_ID   = null_id,
      CHR       = chr_i,
      side      = "right",
      order     = 1:n,
      win_start = as.integer(right_start),
      win_end   = as.integer(right_end)
    )
    
    rbind(left_dt, right_dt)
  }), use.names = TRUE)
  
  out[]
}

# ============================================================
# 3) Count SNP totals and significant SNPs in each 50 kb window
#    Then compute window-level false positive rate:
#    FPR = n_sig / n_total
# ============================================================
count_window_fpr <- function(gwas_df, ordered_null50_df, alpha = 5e-8) {
  gwas <- copy(gwas_df)
  win  <- copy(ordered_null50_df)
  
  setDT(gwas)
  setDT(win)
  
  if (nrow(win) == 0L) return(data.table())
  
  gwas[, `:=`(
    CHR = as.integer(CHR),
    BP  = as.integer(BP),
    P   = as.numeric(P),
    BP2 = as.integer(BP)
  )]
  
  win[, WIN_ID := .I]
  
  setkey(gwas, CHR, BP, BP2)
  setkey(win, CHR, win_start, win_end)
  
  hits <- foverlaps(
    gwas[, .(CHR, BP, BP2, P, Method2, Pheno)],
    win[, .(CHR, win_start, win_end, WIN_ID, NULL_ID, side, order)],
    by.x = c("CHR", "BP", "BP2"),
    by.y = c("CHR", "win_start", "win_end"),
    nomatch = 0L
  )
  
  counts <- hits[
    ,
    .(
      n_total = .N,
      n_sig   = sum(P <= alpha, na.rm = TRUE)
    ),
    by = .(Method2, Pheno, WIN_ID)
  ]
  
  full_template <- CJ(
    Method2 = unique(as.character(gwas_df$Method2)),
    Pheno   = unique(gwas_df$Pheno),
    WIN_ID  = win$WIN_ID,
    unique = TRUE
  )
  
  out <- merge(
    full_template,
    counts,
    by = c("Method2", "Pheno", "WIN_ID"),
    all.x = TRUE
  )
  
  out <- merge(
    out,
    win[, .(WIN_ID, NULL_ID, CHR, side, order, win_start, win_end)],
    by = "WIN_ID",
    all.x = TRUE
  )
  
  out[is.na(n_total), n_total := 0L]
  out[is.na(n_sig),   n_sig   := 0L]
  out[, fpr := fifelse(n_total > 0, n_sig / n_total, NA_real_)]
  
  out[]
}

# ============================================================
# 4) Average FPR across all windows with the same order
# ============================================================
summarize_same_order_fpr <- function(scored_df) {
  x <- copy(scored_df)
  setDT(x)
  
  x[
    ,
    .(
      mean_fpr = mean(fpr, na.rm = TRUE),
      se_fpr   = sd(fpr, na.rm = TRUE) / sqrt(sum(!is.na(fpr))),
      n_bins   = sum(!is.na(fpr))
    ),
    by = .(Method2, order)
  ][order(Method2, order)]
}

summarize_same_order_fpr_by_pheno <- function(scored_df) {
  x <- copy(scored_df)
  setDT(x)
  
  x[
    ,
    .(
      mean_fpr = mean(fpr, na.rm = TRUE),
      se_fpr   = sd(fpr, na.rm = TRUE) / sqrt(sum(!is.na(fpr))),
      n_bins   = sum(!is.na(fpr))
    ),
    by = .(Method2, Pheno, order)
  ][order(Method2, Pheno, order)]
}

# ============================================================
# 4b) Window-level lambda and observed/expected ratio
# ============================================================

compute_lambda <- function(pvals) {
  pvals <- as.numeric(pvals)
  pvals <- pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0 & pvals <= 1]
  
  if (length(pvals) == 0L) return(NA_real_)
  
  chisq <- qchisq(1 - pvals, df = 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
}

count_window_calibration <- function(gwas_df, ordered_null50_df, alpha = 5e-8) {
  gwas <- copy(gwas_df)
  win  <- copy(ordered_null50_df)
  
  setDT(gwas)
  setDT(win)
  
  if (nrow(win) == 0L) return(data.table())
  
  gwas[, `:=`(
    CHR = as.integer(CHR),
    BP  = as.integer(BP),
    P   = as.numeric(P),
    BP2 = as.integer(BP)
  )]
  
  win[, WIN_ID := .I]
  
  setkey(gwas, CHR, BP, BP2)
  setkey(win, CHR, win_start, win_end)
  
  hits <- foverlaps(
    gwas[, .(CHR, BP, BP2, P, Method2, Pheno)],
    win[, .(CHR, win_start, win_end, WIN_ID, NULL_ID, side, order)],
    by.x = c("CHR", "BP", "BP2"),
    by.y = c("CHR", "win_start", "win_end"),
    nomatch = 0L
  )
  
  counts <- hits[
    ,
    .(
      n_total   = .N,
      n_sig     = sum(P <= alpha, na.rm = TRUE),
      lambda    = compute_lambda(P)
    ),
    by = .(Method2, Pheno, WIN_ID)
  ]
  
  full_template <- CJ(
    Method2 = unique(as.character(gwas_df$Method2)),
    Pheno   = unique(gwas_df$Pheno),
    WIN_ID  = win$WIN_ID,
    unique = TRUE
  )
  
  out <- merge(
    full_template,
    counts,
    by = c("Method2", "Pheno", "WIN_ID"),
    all.x = TRUE
  )
  
  out <- merge(
    out,
    win[, .(WIN_ID, NULL_ID, CHR, side, order, win_start, win_end)],
    by = "WIN_ID",
    all.x = TRUE
  )
  
  out[is.na(n_total), n_total := 0L]
  out[is.na(n_sig),   n_sig   := 0L]
  
  out[, fpr := fifelse(n_total > 0, n_sig / n_total, NA_real_)]
  out[, exp_sig := n_total * alpha]
  out[, obs_exp_ratio := fifelse(exp_sig > 0, n_sig / exp_sig, NA_real_)]
  
  out[]
}

summarize_same_order_lambda <- function(scored_df) {
  x <- copy(scored_df)
  setDT(x)
  
  x[
    ,
    .(
      mean_lambda = mean(lambda, na.rm = TRUE),
      se_lambda   = sd(lambda, na.rm = TRUE) / sqrt(sum(!is.na(lambda))),
      n_bins      = sum(!is.na(lambda))
    ),
    by = .(Method2, order)
  ][order(Method2, order)]
}

summarize_same_order_lambda_by_pheno <- function(scored_df) {
  x <- copy(scored_df)
  setDT(x)
  
  x[
    ,
    .(
      mean_lambda = mean(lambda, na.rm = TRUE),
      se_lambda   = sd(lambda, na.rm = TRUE) / sqrt(sum(!is.na(lambda))),
      n_bins      = sum(!is.na(lambda))
    ),
    by = .(Method2, Pheno, order)
  ][order(Method2, Pheno, order)]
}

summarize_same_order_obs_exp <- function(scored_df) {
  x <- copy(scored_df)
  setDT(x)
  
  x[
    ,
    .(
      mean_ratio = mean(obs_exp_ratio, na.rm = TRUE),
      se_ratio   = sd(obs_exp_ratio, na.rm = TRUE) / sqrt(sum(!is.na(obs_exp_ratio))),
      n_bins     = sum(!is.na(obs_exp_ratio))
    ),
    by = .(Method2, order)
  ][order(Method2, order)]
}

summarize_same_order_obs_exp_by_pheno <- function(scored_df) {
  x <- copy(scored_df)
  setDT(x)
  
  x[
    ,
    .(
      mean_ratio = mean(obs_exp_ratio, na.rm = TRUE),
      se_ratio   = sd(obs_exp_ratio, na.rm = TRUE) / sqrt(sum(!is.na(obs_exp_ratio))),
      n_bins     = sum(!is.na(obs_exp_ratio))
    ),
    by = .(Method2, Pheno, order)
  ][order(Method2, Pheno, order)]
}

# ============================================================
# 5) Read effect windows
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

null_4mb_df <- make_null_windows_from_adjacent_effects(region_df)
null_50kb_ordered_df <- make_ordered_50kb_windows(null_4mb_df, win_size = 50000L)

cat("\nNumber of 4 Mb null windows:", nrow(null_4mb_df), "\n")
cat("Number of ordered 50 kb windows:", nrow(null_50kb_ordered_df), "\n")

if (nrow(null_4mb_df) == 0L) {
  stop("No 4 Mb null windows were found between adjacent effect windows.")
}

# ============================================================
# 6) Readers for GWAS results
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
    CHR = as.integer(Chr),
    BP  = as.integer(Pos),
    P   = as.numeric(P.value)
  )]
  
  dt2[, `:=`(
    Method2 = paste0(method_name, "_", loc),
    Pheno   = pheno
  )]
  
  dt2
}

read_bolt <- function(file) {
  dt <- fread(file)
  setDT(dt)
  
  pheno <- extract_pheno(basename(file))
  if (is.na(pheno)) return(NULL)
  
  dt2 <- dt[, .(
    CHR = as.integer(CHR),
    BP  = as.integer(BP),
    P   = as.numeric(P_BOLT_LMM_INF)
  )]
  
  dt2[, `:=`(
    Method2 = "BOLT",
    Pheno   = pheno
  )]
  
  dt2
}

read_fastgwa <- function(file) {
  dt <- fread(file)
  setDT(dt)
  
  pheno <- extract_pheno(basename(file))
  if (is.na(pheno)) return(NULL)
  
  dt2 <- dt[, .(
    CHR = as.integer(CHR),
    BP  = as.integer(POS),
    P   = as.numeric(P)
  )]
  
  dt2[, `:=`(
    Method2 = "fastGWA",
    Pheno   = pheno
  )]
  
  dt2
}

read_regenie <- function(file, method_name) {
  dt <- fread(file)
  setDT(dt)
  
  pheno <- extract_pheno(basename(file))
  if (is.na(pheno)) return(NULL)
  
  dt2 <- dt[, .(
    CHR = as.integer(CHROM),
    BP  = as.integer(GENPOS),
    P   = 10^(-as.numeric(LOG10P))
  )]
  
  dt2[, `:=`(
    Method2 = method_name,
    Pheno   = pheno
  )]
  
  dt2
}

read_slemm <- function(file) {
  dt <- fread(file)
  setDT(dt)
  
  pheno <- extract_pheno(basename(file))
  if (is.na(pheno)) return(NULL)
  
  dt2 <- dt[, .(
    CHR = as.integer(get("#CHROM")),
    BP  = as.integer(POS),
    P   = as.numeric(P)
  )]
  
  dt2[, `:=`(
    Method2 = "SLEMM",
    Pheno   = pheno
  )]
  
  dt2
}

# ============================================================
# 7) Collect GWAS files
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
gwas_all <- gwas_all[!is.na(Pheno) & !is.na(CHR) & !is.na(BP) & !is.na(P)]
gwas_all[, Method2 := factor(as.character(Method2), levels = method_order)]

cat("\nLoaded records by Method2:\n")
print(gwas_all[, .N, by = Method2][order(Method2)])

# ============================================================
# 8) Calculate window-level FPR
# ============================================================
window_fpr <- count_window_fpr(
  gwas_df = gwas_all,
  ordered_null50_df = null_50kb_ordered_df,
  alpha = alpha_sig
)

fwrite(window_fpr, file.path(out_dir, "null50kb_window_fpr_per_window.csv"))

# ============================================================
# 8) Calculate window-level calibration metrics
#    FPR, lambda, observed/expected ratio
# ============================================================
window_calibration <- count_window_calibration(
  gwas_df = gwas_all,
  ordered_null50_df = null_50kb_ordered_df,
  alpha = alpha_sig
)

fwrite(window_calibration, file.path(out_dir, "null50kb_window_calibration_per_window.csv"))

# ============================================================
# 9) Average FPR by order
# ============================================================
fpr_by_order <- summarize_same_order_fpr(window_fpr)
fpr_by_order[, Method2 := factor(as.character(Method2), levels = method_order)]

fpr_by_order_pheno <- summarize_same_order_fpr_by_pheno(window_fpr)
fpr_by_order_pheno[, Method2 := factor(as.character(Method2), levels = method_order)]

fwrite(fpr_by_order, file.path(out_dir, "null50kb_window_fpr_by_order_avg.csv"))
fwrite(fpr_by_order_pheno, file.path(out_dir, "null50kb_window_fpr_by_order_avg_by_pheno.csv"))

# ============================================================
# 9) Average calibration metrics by order
# ============================================================
fpr_by_order <- summarize_same_order_fpr(window_calibration)
fpr_by_order[, Method2 := factor(as.character(Method2), levels = method_order)]

fpr_by_order_pheno <- summarize_same_order_fpr_by_pheno(window_calibration)
fpr_by_order_pheno[, Method2 := factor(as.character(Method2), levels = method_order)]

lambda_by_order <- summarize_same_order_lambda(window_calibration)
lambda_by_order[, Method2 := factor(as.character(Method2), levels = method_order)]

lambda_by_order_pheno <- summarize_same_order_lambda_by_pheno(window_calibration)
lambda_by_order_pheno[, Method2 := factor(as.character(Method2), levels = method_order)]

obs_exp_by_order <- summarize_same_order_obs_exp(window_calibration)
obs_exp_by_order[, Method2 := factor(as.character(Method2), levels = method_order)]

obs_exp_by_order_pheno <- summarize_same_order_obs_exp_by_pheno(window_calibration)
obs_exp_by_order_pheno[, Method2 := factor(as.character(Method2), levels = method_order)]

fwrite(fpr_by_order, file.path(out_dir, "null50kb_window_fpr_by_order_avg.csv"))
fwrite(fpr_by_order_pheno, file.path(out_dir, "null50kb_window_fpr_by_order_avg_by_pheno.csv"))

fwrite(lambda_by_order, file.path(out_dir, "null50kb_window_lambda_by_order_avg.csv"))
fwrite(lambda_by_order_pheno, file.path(out_dir, "null50kb_window_lambda_by_order_avg_by_pheno.csv"))

fwrite(obs_exp_by_order, file.path(out_dir, "null50kb_window_obs_exp_ratio_by_order_avg.csv"))
fwrite(obs_exp_by_order_pheno, file.path(out_dir, "null50kb_window_obs_exp_ratio_by_order_avg_by_pheno.csv"))

# ============================================================
# 10) Plot overall average FPR
#     Plot 1: all methods, all 50 kb windows
# ============================================================
fpr_by_order_plot_all <- copy(fpr_by_order)

if (nrow(fpr_by_order_plot_all) == 0L) {
  stop("No data available for the all-methods FPR plot.")
}

y_max_all <- max(
  fpr_by_order_plot_all$mean_fpr + fpr_by_order_plot_all$se_fpr,
  na.rm = TRUE
)

p_order_avg_all <- ggplot(
  fpr_by_order_plot_all,
  aes(x = order, y = mean_fpr, color = Method2, group = Method2)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  #geom_errorbar(
  #  aes(ymin = pmax(mean_fpr - se_fpr, 0), ymax = mean_fpr + se_fpr),
  #  width = 0.3,
  #  linewidth = 0.3
  #) +
  scale_color_manual(
    values = method_colors,
    drop = FALSE,
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = seq(min(fpr_by_order_plot_all$order), max(fpr_by_order_plot_all$order), by = 2)
  ) +
  scale_y_continuous(limits = c(0, y_max_all * 1.05)) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  labs(
    title = paste0("Average false positive rate in ordered 50 kb null windows (P <= ", alpha_sig, ")"),
    subtitle = "All methods, all 50 kb windows",
    x = "50 kb window order",
    y = "False positive rate (FPR)",
    color = "Method"
  )

print(p_order_avg_all)
save_plot(
  p_order_avg_all,
  "null50kb_average_fpr_by_order_all_methods_all_windows_no_errbar.png",
  w = 10,
  h = 6
)

# ============================================================
# 11) Plot overall average FPR
#     Plot 2: exclude BOLT and REGENIE, all 50 kb windows
# ============================================================
fpr_by_order_plot_no_bolt_regenie <- fpr_by_order[
  !Method2 %in% c("BOLT", "REGENIE")
]

fpr_by_order_plot_no_bolt_regenie[, Method2 := droplevels(Method2)]

if (nrow(fpr_by_order_plot_no_bolt_regenie) == 0L) {
  stop("No data available for the FPR plot after removing BOLT and REGENIE.")
}

y_max_no_bolt_regenie <- max(
  fpr_by_order_plot_no_bolt_regenie$mean_fpr + fpr_by_order_plot_no_bolt_regenie$se_fpr,
  na.rm = TRUE
)

p_order_avg_no_bolt_regenie <- ggplot(
  fpr_by_order_plot_no_bolt_regenie,
  aes(x = order, y = mean_fpr, color = Method2, group = Method2)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  #geom_errorbar(
  #  aes(ymin = pmax(mean_fpr - se_fpr, 0), ymax = mean_fpr + se_fpr),
  #  width = 0.3,
  #  linewidth = 0.3
  #) +
  scale_color_manual(
    values = method_colors,
    drop = FALSE,
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = seq(min(fpr_by_order_plot_no_bolt_regenie$order), max(fpr_by_order_plot_no_bolt_regenie$order), by = 2)
  ) +
  scale_y_continuous(limits = c(0, y_max_no_bolt_regenie * 1.05)) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  labs(
    title = paste0("Average false positive rate in ordered 50 kb null windows (P <= ", alpha_sig, ")"),
    subtitle = "All 50 kb windows, excluding BOLT and REGENIE",
    x = "50 kb window order",
    y = "False positive rate (FPR)",
    color = "Method"
  )

print(p_order_avg_no_bolt_regenie)
save_plot(
  p_order_avg_no_bolt_regenie,
  "null50kb_average_fpr_by_order_no_BOLT_REGENIE_all_windows_no_errbar.png",
  w = 10,
  h = 6
)

# ============================================================
# 12) Plot overall average lambda
#     Plot 1: all methods, all 50 kb windows
# ============================================================
lambda_by_order_plot_all <- copy(lambda_by_order)

if (nrow(lambda_by_order_plot_all) == 0L) {
  stop("No data available for the all-methods lambda plot.")
}

y_max_lambda_all <- max(
  lambda_by_order_plot_all$mean_lambda + lambda_by_order_plot_all$se_lambda,
  na.rm = TRUE
)

p_lambda_avg_all <- ggplot(
  lambda_by_order_plot_all,
  aes(x = order, y = mean_lambda, color = Method2, group = Method2)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(ymin = pmax(mean_lambda - se_lambda, 0), ymax = mean_lambda + se_lambda),
    width = 0.3,
    linewidth = 0.3
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(
    values = method_colors,
    drop = FALSE,
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = seq(min(lambda_by_order_plot_all$order), max(lambda_by_order_plot_all$order), by = 2)
  ) +
  scale_y_continuous(limits = c(0, y_max_lambda_all * 1.05)) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  labs(
    title = "Average genomic inflation factor in ordered 50 kb null windows",
    subtitle = "All methods, all 50 kb windows",
    x = "50 kb window order",
    y = expression(lambda),
    color = "Method"
  )

print(p_lambda_avg_all)
save_plot(
  p_lambda_avg_all,
  "null50kb_average_lambda_by_order_all_methods_all_windows.png",
  w = 10,
  h = 6
)

# ============================================================
# 13) Plot overall average lambda
#     Plot 2: exclude BOLT and REGENIE, all 50 kb windows
# ============================================================
lambda_by_order_plot_no_bolt_regenie <- lambda_by_order[
  !Method2 %in% c("BOLT", "REGENIE")
]

lambda_by_order_plot_no_bolt_regenie[, Method2 := droplevels(Method2)]

if (nrow(lambda_by_order_plot_no_bolt_regenie) == 0L) {
  stop("No data available for the lambda plot after removing BOLT and REGENIE.")
}

y_max_lambda_no_bolt_regenie <- max(
  lambda_by_order_plot_no_bolt_regenie$mean_lambda + lambda_by_order_plot_no_bolt_regenie$se_lambda,
  na.rm = TRUE
)

p_lambda_avg_no_bolt_regenie <- ggplot(
  lambda_by_order_plot_no_bolt_regenie,
  aes(x = order, y = mean_lambda, color = Method2, group = Method2)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(ymin = pmax(mean_lambda - se_lambda, 0), ymax = mean_lambda + se_lambda),
    width = 0.3,
    linewidth = 0.3
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_color_manual(
    values = method_colors,
    drop = FALSE,
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = seq(min(lambda_by_order_plot_no_bolt_regenie$order), max(lambda_by_order_plot_no_bolt_regenie$order), by = 2)
  ) +
  scale_y_continuous(limits = c(0, y_max_lambda_no_bolt_regenie * 1.05)) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  labs(
    title = "Average genomic inflation factor in ordered 50 kb null windows",
    subtitle = "All 50 kb windows, excluding BOLT and REGENIE",
    x = "50 kb window order",
    y = expression(lambda),
    color = "Method"
  )

print(p_lambda_avg_no_bolt_regenie)
save_plot(
  p_lambda_avg_no_bolt_regenie,
  "null50kb_average_lambda_by_order_no_BOLT_REGENIE_all_windows.png",
  w = 10,
  h = 6
)

# ============================================================
# Global font settings (adjust here if needed)
# ============================================================
base_font_size <- 14
title_size <- 14
axis_title_size <- 14
axis_text_size <- 14
legend_title_size <- 14
legend_text_size <- 14

# ============================================================
# 14) Plot 1: all methods
# ============================================================
obs_exp_by_order_plot_all <- copy(obs_exp_by_order)

if (nrow(obs_exp_by_order_plot_all) == 0L) {
  stop("No data available for the all-methods observed/expected ratio plot.")
}

y_max_ratio_all <- max(
  obs_exp_by_order_plot_all$mean_ratio + obs_exp_by_order_plot_all$se_ratio,
  na.rm = TRUE
)

p_obs_exp_avg_all <- ggplot(
  obs_exp_by_order_plot_all,
  aes(x = order, y = mean_ratio, color = Method2, group = Method2)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(
    values = method_colors,
    drop = FALSE,
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = seq(min(obs_exp_by_order_plot_all$order),
                 max(obs_exp_by_order_plot_all$order), by = 2)
  ) +
  scale_y_continuous(limits = c(0, y_max_ratio_all * 1.05)) +
  theme_pub() +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(size = axis_text_size),
    plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    legend.text  = element_text(size = legend_text_size)
  ) +
  labs(
    title = paste0("Average observed/expected ratio in ordered 50 kb null windows (P <= ", alpha_sig, ")"),
    subtitle = "All methods, all 50 kb windows",
    x = "50 kb window order",
    y = "Observed / expected ratio",
    color = "Method"
  )

print(p_obs_exp_avg_all)

save_plot(
  p_obs_exp_avg_all,
  "null50kb_average_obs_exp_ratio_by_order_all_methods_all_windows_no_errorbar.png",
  w = 11,
  h = 7
)

# ============================================================
# 15) Plot 2: exclude BOLT and REGENIE
# ============================================================
obs_exp_by_order_plot_no_bolt_regenie <- copy(
  obs_exp_by_order[!Method2 %in% c("BOLT", "REGENIE")]
)

obs_exp_by_order_plot_no_bolt_regenie[, Method2 := droplevels(Method2)]

if (nrow(obs_exp_by_order_plot_no_bolt_regenie) == 0L) {
  stop("No data available after removing BOLT and REGENIE.")
}

y_max_ratio_no_bolt_regenie <- max(
  obs_exp_by_order_plot_no_bolt_regenie$mean_ratio +
    obs_exp_by_order_plot_no_bolt_regenie$se_ratio,
  na.rm = TRUE
)

p_obs_exp_avg_no_bolt_regenie <- ggplot(
  obs_exp_by_order_plot_no_bolt_regenie,
  aes(x = order, y = mean_ratio, color = Method2, group = Method2)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(
    values = method_colors,
    drop = FALSE,
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = seq(min(obs_exp_by_order_plot_no_bolt_regenie$order),
                 max(obs_exp_by_order_plot_no_bolt_regenie$order), by = 2)
  ) +
  scale_y_continuous(limits = c(0, y_max_ratio_no_bolt_regenie * 1.05)) +
  theme_pub() +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(size = axis_text_size),
    plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    legend.text  = element_text(size = legend_text_size)
  ) +
  labs(
    title = paste0("Average observed/expected ratio in ordered 50 kb null windows (P <= ", alpha_sig, ")"),
    subtitle = "All 50 kb windows, excluding BOLT and REGENIE",
    x = "50 kb window order",
    y = "Observed / expected ratio",
    color = "Method"
  )

print(p_obs_exp_avg_no_bolt_regenie)

save_plot(
  p_obs_exp_avg_no_bolt_regenie,
  "null50kb_average_obs_exp_ratio_by_order_no_BOLT_REGENIE_all_windows_no_errorbar.png",
  w = 11,
  h = 7
)

# ============================================================
# Combined plot with shared legend (large fonts)
# ============================================================
library(ggpubr)

# ============================================================
# Global font settings (adjust here if needed)
# ============================================================
base_font_size <- 18
title_size <- 20
axis_title_size <- 18
axis_text_size <- 18
legend_title_size <- 18
legend_text_size <- 18

# ============================================================
# 14) Plot 1: all methods
# ============================================================
obs_exp_by_order_plot_all <- copy(obs_exp_by_order)

if (nrow(obs_exp_by_order_plot_all) == 0L) {
  stop("No data available for the all-methods observed/expected ratio plot.")
}

y_max_ratio_all <- max(
  obs_exp_by_order_plot_all$mean_ratio + obs_exp_by_order_plot_all$se_ratio,
  na.rm = TRUE
)

p_obs_exp_avg_all <- ggplot(
  obs_exp_by_order_plot_all,
  aes(x = order, y = mean_ratio, color = Method2, group = Method2)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(
    values = method_colors,
    drop = FALSE,
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = seq(min(obs_exp_by_order_plot_all$order),
                 max(obs_exp_by_order_plot_all$order), by = 2)
  ) +
  scale_y_continuous(limits = c(0, y_max_ratio_all * 1.05)) +
  theme_pub() +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(size = axis_text_size),
    plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    legend.text  = element_text(size = legend_text_size)
  ) +
  labs(
    title = paste0("Average observed/expected ratio in ordered 50 kb null windows (P <= ", alpha_sig, ")"),
    subtitle = "All methods, all 50 kb windows",
    x = "50 kb window order",
    y = "Observed / expected ratio",
    color = "Method"
  )

print(p_obs_exp_avg_all)

save_plot(
  p_obs_exp_avg_all,
  "null50kb_average_obs_exp_ratio_by_order_all_methods_all_windows_no_errorbar.png",
  w = 11,
  h = 7
)

# ============================================================
# 15) Plot 2: exclude BOLT and REGENIE
# ============================================================
obs_exp_by_order_plot_no_bolt_regenie <- copy(
  obs_exp_by_order[!Method2 %in% c("BOLT", "REGENIE")]
)

obs_exp_by_order_plot_no_bolt_regenie[, Method2 := droplevels(Method2)]

if (nrow(obs_exp_by_order_plot_no_bolt_regenie) == 0L) {
  stop("No data available after removing BOLT and REGENIE.")
}

y_max_ratio_no_bolt_regenie <- max(
  obs_exp_by_order_plot_no_bolt_regenie$mean_ratio +
    obs_exp_by_order_plot_no_bolt_regenie$se_ratio,
  na.rm = TRUE
)

p_obs_exp_avg_no_bolt_regenie <- ggplot(
  obs_exp_by_order_plot_no_bolt_regenie,
  aes(x = order, y = mean_ratio, color = Method2, group = Method2)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(
    values = method_colors,
    drop = FALSE,
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = seq(min(obs_exp_by_order_plot_no_bolt_regenie$order),
                 max(obs_exp_by_order_plot_no_bolt_regenie$order), by = 2)
  ) +
  scale_y_continuous(limits = c(0, y_max_ratio_no_bolt_regenie * 1.05)) +
  theme_pub() +
  theme(
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(size = axis_text_size),
    plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = axis_text_size),
    legend.title = element_text(size = legend_title_size),
    legend.text  = element_text(size = legend_text_size)
  ) +
  labs(
    title = paste0("Average observed/expected ratio in ordered 50 kb null windows (P <= ", alpha_sig, ")"),
    subtitle = "All 50 kb windows, excluding BOLT and REGENIE",
    x = "50 kb window order",
    y = "Observed / expected ratio",
    color = "Method"
  )

print(p_obs_exp_avg_no_bolt_regenie)

save_plot(
  p_obs_exp_avg_no_bolt_regenie,
  "null50kb_average_obs_exp_ratio_by_order_no_BOLT_REGENIE_all_windows_no_errorbar.png",
  w = 11,
  h = 7
)

# ============================================================
# Combined plot with shared legend (large fonts)
# ============================================================

# remove legends from individual plots? not necessary for ggarrange
# but keep them normal

p_obs_exp_combined <- ggarrange(
  p_obs_exp_avg_all,
  p_obs_exp_avg_no_bolt_regenie,
  ncol = 1,
  nrow = 2,
  align = "v",
  common.legend = TRUE,
  legend = "bottom",
  labels = c("A", "B"),
  label.x = 0.01,   # move left/right (0 = far left)
  label.y = 1,   # move up/down (1 = top)
  font.label = list(size = 20, face = "bold")
)

print(p_obs_exp_combined)

ggsave(
  filename = "null50kb_average_obs_exp_ratio_combined_shared_legend.png",
  plot = p_obs_exp_combined,
  width = 13,
  height = 11,
  dpi = 300,
  bg = "white"
)
