# L4: harmonize GWAS statistics with the LOCO-adjusted LD (L2) and run SuSiE-adj.
# Same harmonization, locus set, SNP filtering and fine-mapping parameters as scripts/04_susie_adj.R, except:
#   * LD / n_eff come from the method-specific LOCO adjustment (loco/ld/<METHOD>/<pheno>/<locus>/).
#   * PRIMARY fit: estimate_residual_variance = TRUE (SuSiE-adj authors' example and jwang259's production
#     script; susieR 0.14.2's own susie_rss default is FALSE). No fallback: failures are recorded as failures.
#   * Loci whose primary fit errors get a separate SENSITIVITY fit with estimate_residual_variance = FALSE,
#     written to sens_ERV_FALSE_* files and never mixed into the primary results.
# The causal SNP list is NOT read here.
# Usage: Rscript L4_susie_adj_loco.R [NCORES] [test]
suppressMessages({ library(data.table); library(susieR); library(parallel) })
stopifnot(as.character(packageVersion("susieR")) == "0.14.2")
setDTthreads(1)

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1]))
W <- dirname(dirname(dirname(script_path)))
L <- file.path(W, "loco")
BENCH <- dirname(W)
PAR <- list(L = 10L, coverage = 0.95, min_abs_corr = 0.5, estimate_residual_variance = TRUE, max_iter = 1000L)
AMBIG <- c("AT", "TA", "CG", "GC")
MDIR <- c("BOLT-LMM" = "BOLT_LMM", "REGENIE-LOCO" = "REGENIE")

bim <- fread(file.path(W, "genotype/spg.20k.bim"), header = FALSE, col.names = c("chr", "SNP", "cm", "position", "A1", "A2"))
bim[, cm := NULL]

read_gwas <- function(method, p) {
  if (method == "BOLT-LMM") {
    d <- fread(sprintf("%s/BOLTLMM/bolt_pheno%d.stats", BENCH, p))
    d <- d[, .(SNP, gw_chr = CHR, gw_pos = BP, eff = ALLELE1, oth = ALLELE0, beta = BETA, SE, P = P_BOLT_LMM_INF,
               log10P = -log10(P_BOLT_LMM_INF))]
  } else {
    d <- fread(sprintf("%s/REGENIE/regenie_pheno%d_pheno%d.regenie", BENCH, p, p))
    # REGENIE IDs lack the allele suffix: map to the panel ID by chr:pos, alleles checked below
    d <- d[, .(gw_chr = CHROM, gw_pos = GENPOS, gw_id = ID, eff = ALLELE1, oth = ALLELE0, beta = BETA, SE,
               log10P = LOG10P, P = 10^(-LOG10P))]
    d <- merge(d, bim[, .(gw_chr = chr, gw_pos = position, SNP)], by = c("gw_chr", "gw_pos"), all.x = TRUE)
  }
  d
}

classify <- function(msg) fcase(
  grepl("residual variance", msg, ignore.case = TRUE), "residual_variance_estimation",
  grepl("positive (semi-)?definite|eigen", msg, ignore.case = TRUE), "LD_matrix",
  grepl("NA|NaN|finite", msg), "non_finite",
  default = "other_error")

fit_locus <- function(loc, gw, ld_dir, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ld <- fread(file.path(ld_dir, "snps.tsv"))
  ldsum <- fread(file.path(ld_dir, "summary.tsv"))
  m <- nrow(ld); n_eff <- ldsum$n_eff
  R <- matrix(readBin(file.path(ld_dir, "R.bin"), "double", m * m), m, m)
  stopifnot(ldsum$n_snps == m, isSymmetric(R, tol = 1e-10), all(abs(diag(R) - 1) < 1e-12))
  req <- bim[chr == loc$chr & position >= loc$start & position <= loc$end]

  x <- merge(req, ld[, .(SNP, ld_index, counted)], by = "SNP", all.x = TRUE)
  x <- merge(x, gw[SNP %in% req$SNP], by = "SNP", all.x = TRUE)
  x[, dup := duplicated(SNP) | duplicated(SNP, fromLast = TRUE)]
  x[, reason := fcase(dup, "duplicated_id",
                      is.na(ld_index), "not_in_adjusted_LD(invariant)",
                      is.na(beta), "missing_in_GWAS",
                      gw_chr != chr | gw_pos != position, "position_mismatch",
                      !((eff == A1 & oth == A2) | (eff == A2 & oth == A1)), "allele_mismatch",
                      !is.finite(beta) | !is.finite(SE), "nonfinite_beta_or_SE",
                      SE <= 0, "invalid_SE",
                      default = "OK")]
  x[, ambiguous := paste0(A1, A2) %in% AMBIG]
  x[, flipped := reason == "OK" & eff != counted]
  k <- x[reason == "OK"][order(ld_index)]
  k[, `:=`(beta_h = ifelse(flipped, -beta, beta), effect_allele = counted,
           other_allele = ifelse(counted == A1, A2, A1))]
  k[, z := beta_h / SE]
  qc <- data.table(locus_id = loc$job_id, method = loc$method, phenotype = loc$phenotype, chr = loc$chr,
                   requested = nrow(req), matched = nrow(k), dropped = nrow(req) - nrow(k),
                   drop_reasons = if (any(x$reason != "OK")) x[reason != "OK", paste(names(table(reason)), table(reason), sep = "=", collapse = ";")] else "none",
                   allele_flips = sum(k$flipped), ambiguous_kept = sum(k$ambiguous), n_eff = n_eff, eta = ldsum$eta,
                   eta_source = ldsum$eta_source, max_abs_z = max(abs(k$z)))
  fwrite(x[, .(SNP, chr, position, A1, A2, ld_index, counted, eff, oth, reason, ambiguous, flipped)],
         file.path(out_dir, "harmonization_qc.tsv"), sep = "\t")

  Rk <- R[k$ld_index, k$ld_index, drop = FALSE]
  stopifnot(identical(ld$SNP[k$ld_index], k$SNP))                  # z and R in the same order
  run_susie <- function(erv) {
    warn <- character(0); t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(withCallingHandlers(suppressMessages(          # drops susieR's static "HINT: ... in-sample LD" message
      susie_rss(z = k$z, R = Rk, n = n_eff, L = PAR$L, coverage = PAR$coverage, min_abs_corr = PAR$min_abs_corr,
                estimate_residual_variance = erv, max_iter = PAR$max_iter)),
      warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    list(fit = fit, warn = paste(unique(warn), collapse = " | "), secs = proc.time()[["elapsed"]] - t0)
  }
  # Diagnostic only (does not alter the fit): z/LD consistency, 0 = consistent
  qc[, s_rss := tryCatch(suppressWarnings(estimate_s_rss(k$z, Rk, n = n_eff)), error = function(e) NA_real_)]

  write_fit <- function(fit, prefix) {
    cs <- fit$sets$cs; cs_id <- rep(NA_character_, nrow(k)); cs_tab <- data.table()
    if (length(cs)) {
      pur <- as.data.table(fit$sets$purity)
      cs_tab <- rbindlist(lapply(seq_along(cs), function(j) data.table(
        locus_id = loc$job_id, cs = names(cs)[j], cs_size = length(cs[[j]]), coverage = fit$sets$coverage[j],
        min_abs_corr = pur$min.abs.corr[j], mean_abs_corr = pur$mean.abs.corr[j], median_abs_corr = pur$median.abs.corr[j],
        SNPs = paste(k$SNP[cs[[j]]], collapse = ","))))
      for (j in seq_along(cs)) cs_id[cs[[j]]] <- ifelse(is.na(cs_id[cs[[j]]]), names(cs)[j], paste(cs_id[cs[[j]]], names(cs)[j], sep = ","))
    }
    snp <- k[, .(phenotype = loc$phenotype, GWAS_method = loc$method, chr, position, SNP, effect_allele, other_allele,
                 beta = beta_h, SE, z, P, log10P, PIP = fit$pip, credible_set = cs_id, locus_id = loc$job_id, flipped, ambiguous)]
    fwrite(snp, file.path(out_dir, paste0(prefix, "snp_pip.tsv")), sep = "\t")
    fwrite(cs_tab, file.path(out_dir, paste0(prefix, "credible_sets.tsv")), sep = "\t")
    length(cs)
  }
  status_row <- function(r, erv) {
    f <- r$fit
    if (inherits(f, "error")) {
      msg <- conditionMessage(f)
      return(data.table(fit_status = "error", failure_category = classify(msg), message = msg, converged = NA, niter = NA,
                        n_cs = NA_integer_, sigma2 = NA_real_, elbo = NA_real_, seconds = r$secs, warnings = r$warn))
    }
    data.table(fit_status = "ok", failure_category = if (isTRUE(f$converged)) "none" else "not_converged", message = "",
               converged = isTRUE(f$converged), niter = f$niter, n_cs = length(f$sets$cs), sigma2 = f$sigma2,
               elbo = tail(f$elbo, 1), seconds = r$secs, warnings = r$warn)
  }

  prim <- run_susie(PAR$estimate_residual_variance)
  st <- data.table(qc, fit_mode = "ERV_TRUE", status_row(prim))
  if (st$fit_status == "ok") { saveRDS(prim$fit, file.path(out_dir, "fit.rds")); write_fit(prim$fit, "") }
  else writeLines(st$message, file.path(out_dir, "FAILED"))
  sens <- NULL
  if (st$fit_status == "error") {                                   # sensitivity only; never enters primary tables
    sens <- run_susie(FALSE)
    s2 <- data.table(qc, fit_mode = "ERV_FALSE_sensitivity", status_row(sens))
    if (s2$fit_status == "ok") write_fit(sens$fit, "sens_ERV_FALSE_")
    fwrite(s2, file.path(out_dir, "sens_ERV_FALSE_fit_summary.tsv"), sep = "\t")
  }
  fwrite(st, file.path(out_dir, "fit_summary.tsv"), sep = "\t")
  writeLines("DONE", file.path(out_dir, "DONE"))
  st
}

manifest <- fread(file.path(W, "loci/manifest.tsv"))
manifest[, pheno_num := as.integer(sub("pheno", "", phenotype))]
out_path <- function(loc) file.path(L, "susie_results", MDIR[[loc$method]], loc$phenotype, loc$locus_id)
ld_path  <- function(loc) file.path(L, "ld", MDIR[[loc$method]], loc$phenotype, loc$locus_id)

args <- commandArgs(TRUE)
ncores <- if (length(args)) as.integer(args[1]) else 20L
if (length(args) > 1 && args[2] == "test") {                         # a few chr17 loci per method
  manifest <- manifest[chr == 17 & phenotype %in% c("pheno1", "pheno2")]
  cat("TEST MODE:", nrow(manifest), "loci\n")
}
groups <- split(manifest, by = c("method", "phenotype"))
res <- mclapply(groups, function(g) {
  gw <- read_gwas(g$method[1], g$pheno_num[1])
  rbindlist(lapply(seq_len(nrow(g)), function(i) {
    loc <- g[i]; od <- out_path(loc)
    if (file.exists(file.path(od, "DONE"))) return(fread(file.path(od, "fit_summary.tsv")))
    if (!file.exists(file.path(ld_path(loc), "summary.tsv"))) return(data.table(locus_id = loc$job_id, fit_status = "ld_missing"))
    tryCatch(fit_locus(loc, gw, ld_path(loc), od),
             error = function(e) data.table(locus_id = loc$job_id, fit_status = "pipeline_error", message = conditionMessage(e)))
  }), fill = TRUE)
}, mc.cores = ncores, mc.preschedule = FALSE)
st <- rbindlist(res, fill = TRUE)
fwrite(st, file.path(L, if (length(args) > 1 && args[2] == "test") "qc/L4_test_status.tsv" else "qc/L4_fit_status_all.tsv"), sep = "\t")
print(st[, .N, by = .(method, fit_status, failure_category)])
