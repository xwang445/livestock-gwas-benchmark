#!/usr/bin/env Rscript
# Purpose
#   Calculate spillover-locus rates in central-null 500 kb windows.
# Manuscript section or figure
#   Figure 6A
# Required software
#   R with data.table, ggplot2, and scales
# Required inputs
#   GWAS results and effect/null window definitions
# Outputs
#   Spillover summary tables and figures
# Example command
#   Rscript scripts/07_spillover/figure6A_spillover_locus_500kb.R --project-root . --output-dir derived_data/spillover_summary
# Note
#   Proprietary genotype and phenotype inputs are not distributed with this repository.

#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ============================================================
# COMPLETE STANDALONE SCRIPT FOR REVISED FIGURE 6A
#
# Metric:
#   Spillover-locus rate =
#   number of 500-kb windows within central-null regions that
#   contain at least one SNP with P <= 5e-8
#   divided by the total number of central-null 500-kb windows.
#
# The metric is calculated for each method and phenotype, then
# summarized as mean +/- SD across the 10 phenotypes.
#
# Outputs:
#   1. Figure6A_spillover_locus_rate_500kb.png
#   2. Figure6A_spillover_locus_rate_500kb.pdf
#   3. Figure6A_spillover_locus_rate_by_pheno.csv
#   4. Figure6A_spillover_locus_rate_summary.csv
#   5. Figure6A_spillover_500kb_window_status.csv
# ============================================================


# ============================================================
# 0. Working directory and output directory
# ============================================================

setwd(".")

base_dir <- getwd()

out_dir <- file.path(
  base_dir,
  "Figure6A_spillover_locus"
)

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 1. Analysis settings
# ============================================================

genome_wide_threshold <- 5e-8
central_null_size <- 2000000L
locus_window_size <- 500000L

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
  "FarmCPU_Kansas"  = "FarmCPU (Kansas)",
  "FarmCPU_NYC"     = "FarmCPU (NYC)",
  "BLINK_Kansas"    = "BLINK (Kansas)",
  "BLINK_NYC"       = "BLINK (NYC)",
  "BOLT"            = "BOLT-LMM",
  "REGENIE"         = "REGENIE (LOCO)",
  "REGENIE_woLOCO"  = "REGENIE w/o LOCO",
  "fastGWA"         = "fastGWA",
  "SLEMM"           = "SLEMM"
)


# ============================================================
# 2. Helper functions
# ============================================================

extract_pheno <- function(x) {
  ph <- sub(".*pheno([0-9]+).*", "\\1", x)
  
  if (identical(ph, x)) {
    return(NA_integer_)
  }
  
  as.integer(ph)
}

extract_loc <- function(x) {
  if (grepl("\\(Kansas\\)", x)) {
    return("Kansas")
  }
  
  if (grepl("\\(NYC\\)", x)) {
    return("NYC")
  }
  
  NA_character_
}


# Build the central 2-Mb null regions from gaps between
# consecutive effect windows.
make_null_center_regions <- function(
    region_df,
    center_size = 2000000L
) {
  x <- copy(region_df)
  setorder(x, CHR, Start, End)
  
  x[
    ,
    next_start := shift(Start, type = "lead"),
    by = CHR
  ]
  
  null_df <- x[
    !is.na(next_start),
    .(
      CHR,
      gap_start = End + 1L,
      gap_end = next_start - 1L
    )
  ]
  
  null_df <- null_df[
    gap_end >= gap_start
  ]
  
  null_df[
    ,
    gap_width := gap_end - gap_start + 1L
  ]
  
  null_df <- null_df[
    gap_width >= center_size
  ]
  
  null_df[
    ,
    midpoint := floor(
      (gap_start + gap_end) / 2
    )
  ]
  
  null_df[
    ,
    Start := as.integer(
      midpoint - center_size / 2 + 1
    )
  ]
  
  null_df[
    ,
    End := as.integer(
      midpoint + center_size / 2
    )
  ]
  
  null_df[
    Start < gap_start,
    Start := gap_start
  ]
  
  null_df[
    End > gap_end,
    End := gap_end
  ]
  
  null_df[
    ,
    .(
      CHR,
      Start,
      End
    )
  ]
}


# ============================================================
# 3. Read effect-window definition and construct central-null regions
# ============================================================

region_file <- file.path(
  base_dir,
  "eligible_regions.csv"
)

if (!file.exists(region_file)) {
  stop(
    "Cannot find eligible_regions.csv at: ",
    region_file
  )
}

region_df <- fread(region_file)

setnames(
  region_df,
  c("chr", "start", "end"),
  c("CHR", "Start", "End")
)

region_df[
  ,
  `:=`(
    CHR = as.integer(CHR),
    Start = as.integer(Start),
    End = as.integer(End)
  )
]

setorder(
  region_df,
  CHR,
  Start,
  End
)

null_center_df <- make_null_center_regions(
  region_df = region_df,
  center_size = central_null_size
)

null_center_df[
  ,
  NullRegionID := .I
]

cat(
  "\nNumber of central-null 2-Mb regions:",
  nrow(null_center_df),
  "\n"
)


# ============================================================
# 4. Reader functions for each GWAS method
# ============================================================

read_gapit <- function(
    file,
    method_name
) {
  dt <- fread(file)
  
  effect_col <- if ("effect" %in% names(dt)) {
    "effect"
  } else if ("Effect" %in% names(dt)) {
    "Effect"
  } else {
    NA_character_
  }
  
  if (is.na(effect_col)) {
    stop(
      "GAPIT file is missing effect/Effect column: ",
      basename(file)
    )
  }
  
  location <- extract_loc(
    basename(file)
  )
  
  phenotype <- extract_pheno(
    basename(file)
  )
  
  if (is.na(phenotype)) {
    return(NULL)
  }
  
  dt[
    ,
    .(
      SNP = as.character(SNP),
      CHR = as.integer(Chr),
      BP = as.integer(Pos),
      P = as.numeric(P.value),
      Method = method_name,
      Location = location,
      Method2 = paste0(
        method_name,
        "_",
        location
      ),
      Pheno = phenotype
    )
  ]
}


read_bolt <- function(file) {
  dt <- fread(file)
  
  phenotype <- extract_pheno(
    basename(file)
  )
  
  if (is.na(phenotype)) {
    return(NULL)
  }
  
  dt[
    ,
    .(
      SNP = as.character(SNP),
      CHR = as.integer(CHR),
      BP = as.integer(BP),
      P = as.numeric(P_BOLT_LMM_INF),
      Method = "BOLT",
      Location = NA_character_,
      Method2 = "BOLT",
      Pheno = phenotype
    )
  ]
}


read_fastgwa <- function(file) {
  dt <- fread(file)
  
  phenotype <- extract_pheno(
    basename(file)
  )
  
  if (is.na(phenotype)) {
    return(NULL)
  }
  
  dt[
    ,
    .(
      SNP = as.character(SNP),
      CHR = as.integer(CHR),
      BP = as.integer(POS),
      P = as.numeric(P),
      Method = "fastGWA",
      Location = NA_character_,
      Method2 = "fastGWA",
      Pheno = phenotype
    )
  ]
}


read_regenie <- function(
    file,
    method_name
) {
  dt <- fread(file)
  
  phenotype <- extract_pheno(
    basename(file)
  )
  
  if (is.na(phenotype)) {
    return(NULL)
  }
  
  dt[
    ,
    .(
      SNP = as.character(ID),
      CHR = as.integer(CHROM),
      BP = as.integer(GENPOS),
      P = 10^(-as.numeric(LOG10P)),
      Method = method_name,
      Location = NA_character_,
      Method2 = method_name,
      Pheno = phenotype
    )
  ]
}


read_slemm <- function(file) {
  dt <- fread(file)
  
  phenotype <- extract_pheno(
    basename(file)
  )
  
  if (is.na(phenotype)) {
    return(NULL)
  }
  
  dt[
    ,
    .(
      SNP = as.character(ID),
      CHR = as.integer(get("#CHROM")),
      BP = as.integer(POS),
      P = as.numeric(P),
      Method = "SLEMM",
      Location = NA_character_,
      Method2 = "SLEMM",
      Pheno = phenotype
    )
  ]
}


# ============================================================
# 5. Find and read all GWAS result files
# ============================================================

all_list <- list()


blink_files <- list.files(
  file.path(base_dir, "BLINK"),
  pattern = paste0(
    "GAPIT\\.Association\\.GWAS_Results\\.",
    "BLINK\\.pheno.*\\.csv$"
  ),
  full.names = TRUE
)

for (f in blink_files) {
  all_list[[length(all_list) + 1L]] <- read_gapit(
    f,
    "BLINK"
  )
}

farmcpu_files <- list.files(
  file.path(base_dir, "FarmCPU"),
  pattern = paste0(
    "GAPIT\\.Association\\.GWAS_Results\\.",
    "FarmCPU\\.pheno.*\\.csv$"
  ),
  full.names = TRUE
)

for (f in farmcpu_files) {
  all_list[[length(all_list) + 1L]] <- read_gapit(
    f,
    "FarmCPU"
  )
}


bolt_files <- list.files(
  file.path(base_dir, "BOLTLMM"),
  pattern = "^bolt_pheno[0-9]+\\.stats$",
  full.names = TRUE
)

for (f in bolt_files) {
  all_list[[length(all_list) + 1L]] <- read_bolt(f)
}


fastgwa_files <- list.files(
  file.path(base_dir, "fastGWA"),
  pattern = "pheno[0-9]+\\.fastgwa\\.fastGWA$",
  full.names = TRUE
)

for (f in fastgwa_files) {
  all_list[[length(all_list) + 1L]] <- read_fastgwa(f)
}


regenie_files <- list.files(
  file.path(base_dir, "REGENIE"),
  pattern = paste0(
    "^regenie_pheno[0-9]+_",
    "pheno[0-9]+\\.regenie$"
  ),
  full.names = TRUE
)

for (f in regenie_files) {
  all_list[[length(all_list) + 1L]] <- read_regenie(
    f,
    "REGENIE"
  )
}


regenie_no_loco_files <- list.files(
  file.path(base_dir, "REGENIE_WO_LOCO"),
  pattern = paste0(
    "^regenie_pheno_wo_loco[0-9]+_",
    "pheno[0-9]+\\.regenie$"
  ),
  full.names = TRUE
)

for (f in regenie_no_loco_files) {
  all_list[[length(all_list) + 1L]] <- read_regenie(
    f,
    "REGENIE_woLOCO"
  )
}


slemm_files <- list.files(
  file.path(base_dir, "SLEMM"),
  pattern = paste0(
    "spg\\.even\\.100k\\.pheno[0-9]+",
    "\\.allchr\\.gwas\\.txt$"
  ),
  full.names = TRUE
)

for (f in slemm_files) {
  all_list[[length(all_list) + 1L]] <- read_slemm(f)
}


cat("\nFiles found:\n")

print(
  c(
    BLINK = length(blink_files),
    FarmCPU = length(farmcpu_files),
    BOLTLMM = length(bolt_files),
    fastGWA = length(fastgwa_files),
    REGENIE = length(regenie_files),
    REGENIE_WO_LOCO = length(
      regenie_no_loco_files
    ),
    SLEMM = length(slemm_files)
  )
)


all_list <- Filter(
  Negate(is.null),
  all_list
)

if (length(all_list) == 0L) {
  stop(
    "No GWAS result files were found."
  )
}

gwas_all <- rbindlist(
  all_list,
  use.names = TRUE,
  fill = TRUE
)

gwas_all <- gwas_all[
  !is.na(Pheno) &
    !is.na(CHR) &
    !is.na(BP) &
    !is.na(P)
]

gwas_all[
  P <= 0,
  P := 1e-300
]

gwas_all[
  P > 1,
  P := NA_real_
]

gwas_all <- gwas_all[
  !is.na(P)
]

cat("\nLoaded records by method:\n")

print(
  gwas_all[
    ,
    .N,
    by = .(
      Method2,
      Pheno
    )
  ][
    order(
      Method2,
      Pheno
    )
  ]
)


# ============================================================
# 6. Divide every central-null 2-Mb region into 500-kb windows
# ============================================================

null_center_500kb <- null_center_df[
  ,
  {
    window_starts <- seq.int(
      from = Start,
      to = End,
      by = locus_window_size
    )
    
    data.table(
      WindowStart = window_starts,
      WindowEnd = pmin(
        window_starts +
          locus_window_size -
          1L,
        End
      )
    )
  },
  by = .(
    NullRegionID,
    CHR
  )
]

null_center_500kb[
  ,
  WindowWidth :=
    WindowEnd -
    WindowStart +
    1L
]

if (
  any(
    null_center_500kb$WindowWidth !=
    locus_window_size
  )
) {
  stop(
    "At least one central-null region did not divide ",
    "into exact 500-kb windows. Check eligible_regions.csv."
  )
}

null_center_500kb[
  ,
  WindowID := sprintf(
    "chr%s_null%s_%s_%s",
    CHR,
    NullRegionID,
    WindowStart,
    WindowEnd
  )
]

total_null_windows <- uniqueN(
  null_center_500kb$WindowID
)

cat(
  "\nTotal central-null 500-kb windows:",
  total_null_windows,
  "\n"
)

cat(
  "Windows per central-null region:",
  total_null_windows /
    uniqueN(
      null_center_500kb$NullRegionID
    ),
  "\n"
)


# ============================================================
# 7. Assign SNPs to central-null 500-kb windows
# ============================================================

gwas_intervals <- gwas_all[
  ,
  .(
    Method2 = as.character(Method2),
    Pheno = as.integer(Pheno),
    SNP,
    CHR = as.integer(CHR),
    BP = as.integer(BP),
    P = as.numeric(P),
    BP_start = as.integer(BP),
    BP_end = as.integer(BP)
  )
]

null_window_intervals <- null_center_500kb[
  ,
  .(
    WindowID,
    NullRegionID,
    CHR,
    WindowStart,
    WindowEnd
  )
]

setkey(
  null_window_intervals,
  CHR,
  WindowStart,
  WindowEnd
)

setkey(
  gwas_intervals,
  CHR,
  BP_start,
  BP_end
)

gwas_in_null_windows <- foverlaps(
  x = gwas_intervals,
  y = null_window_intervals,
  by.x = c(
    "CHR",
    "BP_start",
    "BP_end"
  ),
  by.y = c(
    "CHR",
    "WindowStart",
    "WindowEnd"
  ),
  type = "within",
  nomatch = 0L
)

cat(
  "\nGWAS records assigned to central-null windows:",
  nrow(gwas_in_null_windows),
  "\n"
)


# ============================================================
# 8. Identify genome-wide-significant spillover windows
# ============================================================

detected_window_results <- gwas_in_null_windows[
  ,
  .(
    minimum_P = min(
      P,
      na.rm = TRUE
    ),
    significant_SNP_count = sum(
      P <= genome_wide_threshold,
      na.rm = TRUE
    ),
    spillover_locus = as.integer(
      any(
        P <= genome_wide_threshold,
        na.rm = TRUE
      )
    )
  ),
  by = .(
    Method2,
    Pheno,
    WindowID
  )
]


# Build complete method x phenotype x window grid.
# Windows without significant SNPs remain in the denominator.
analysis_template <- CJ(
  Method2 = method_order,
  Pheno = sort(
    unique(
      gwas_all$Pheno
    )
  ),
  WindowID = unique(
    null_center_500kb$WindowID
  ),
  unique = TRUE
)

window_status <- merge(
  analysis_template,
  detected_window_results,
  by = c(
    "Method2",
    "Pheno",
    "WindowID"
  ),
  all.x = TRUE,
  sort = FALSE
)

window_status[
  is.na(spillover_locus),
  `:=`(
    spillover_locus = 0L,
    significant_SNP_count = 0L
  )
]


# ============================================================
# 9. Calculate spillover-locus rate for each phenotype
# ============================================================

spillover_by_pheno <- window_status[
  ,
  .(
    total_central_null_windows = .N,
    
    spillover_locus_count =
      sum(
        spillover_locus
      ),
    
    spillover_locus_rate =
      mean(
        spillover_locus
      )
  ),
  by = .(
    Method2,
    Pheno
  )
]

spillover_by_pheno[
  ,
  spillover_locus_percent :=
    100 *
    spillover_locus_rate
]

spillover_by_pheno[
  ,
  Method2 := factor(
    as.character(Method2),
    levels = method_order
  )
]

setorder(
  spillover_by_pheno,
  Method2,
  Pheno
)


# ============================================================
# 10. Summarize mean +/- SD across phenotypes
# ============================================================

spillover_summary <- spillover_by_pheno[
  ,
  .(
    n_phenotypes = .N,
    
    total_central_null_windows =
      unique(
        total_central_null_windows
      )[1],
    
    mean_spillover_locus_count =
      mean(
        spillover_locus_count
      ),
    
    sd_spillover_locus_count =
      sd(
        spillover_locus_count
      ),
    
    mean_spillover_locus_rate =
      mean(
        spillover_locus_rate
      ),
    
    sd_spillover_locus_rate =
      sd(
        spillover_locus_rate
      ),
    
    mean_spillover_locus_percent =
      mean(
        spillover_locus_percent
      ),
    
    sd_spillover_locus_percent =
      sd(
        spillover_locus_percent
      )
  ),
  by = Method2
]

spillover_summary[
  ,
  Method2 := factor(
    as.character(Method2),
    levels = method_order
  )
]

spillover_summary[
  ,
  plot_label := sprintf(
    "%.1f%%\n%.1f windows",
    mean_spillover_locus_percent,
    mean_spillover_locus_count
  )
]

setorder(
  spillover_summary,
  Method2
)


cat(
  "\nSpillover-locus results by phenotype:\n"
)

print(
  spillover_by_pheno
)


cat(
  "\nSpillover-locus summary across phenotypes:\n"
)

print(
  spillover_summary
)


# ============================================================
# 11. Save output tables
# ============================================================

fwrite(
  window_status,
  file.path(
    out_dir,
    "Figure6A_spillover_500kb_window_status.csv"
  )
)

fwrite(
  spillover_by_pheno,
  file.path(
    out_dir,
    "Figure6A_spillover_locus_rate_by_pheno.csv"
  )
)

fwrite(
  spillover_summary,
  file.path(
    out_dir,
    "Figure6A_spillover_locus_rate_summary.csv"
  )
)


# ============================================================
# 12. Make revised Figure 6A
#
# Bars: mean spillover-locus rate
# Error bars: SD across 10 phenotypes
# Points: individual phenotype values
# Labels: mean percentage and mean absolute count per genome
# ============================================================

plot_upper <- max(
  spillover_summary$mean_spillover_locus_rate +
    spillover_summary$sd_spillover_locus_rate,
  na.rm = TRUE
)

plot_upper <- max(
  plot_upper * 1.35,
  0.01
)


figure_6a <- ggplot(
  spillover_summary,
  aes(
    x = Method2,
    y = mean_spillover_locus_rate,
    fill = Method2
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(
        mean_spillover_locus_rate -
          sd_spillover_locus_rate,
        0
      ),
      ymax =
        mean_spillover_locus_rate +
        sd_spillover_locus_rate
    ),
    width = 0.18,
    linewidth = 0.5
  ) +
  geom_point(
    data = spillover_by_pheno,
    aes(
      x = Method2,
      y = spillover_locus_rate
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.10,
      height = 0
    ),
    size = 1.3,
    alpha = 0.70
  ) +
  geom_text(
    aes(
      label = plot_label
    ),
    vjust = -0.45,
    size = 3.3,
    lineheight = 0.9
  ) +
  scale_fill_manual(
    values = method_colors,
    drop = FALSE
  ) +
  scale_x_discrete(
    labels = method_labels,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(
      accuracy = 1
    ),
    limits = c(
      0,
      plot_upper
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  labs(
    title = paste0(
      "Genome-wide spillover loci "
    ),
    subtitle = paste0(
      "Central-null 500-kb windows (P ≤ ",
      format(genome_wide_threshold, scientific = TRUE),
      "); mean ± SD across 10 replicates"
    ),
    x = NULL,
    y = "Spillover-locus rate"
  ) +
  theme_bw(
    base_size = 14
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      linewidth = 0.2
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    axis.title = element_text(
      size = 16
    ),
    plot.title = element_text(
      size = 18,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    plot.margin = margin(
      12,
      12,
      8,
      8
    )
  )


print(
  figure_6a
)


ggsave(
  filename = file.path(
    out_dir,
    "Figure6A_spillover_locus_rate_500kb.png"
  ),
  plot = figure_6a,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    out_dir,
    "Figure6A_spillover_locus_rate_500kb.pdf"
  ),
  plot = figure_6a,
  width = 9,
  height = 6.5
)


cat(
  "\nFinished.\n",
  "Figure and result tables were saved to:\n",
  out_dir,
  "\n",
  sep = ""
)
