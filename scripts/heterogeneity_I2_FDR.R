#!/usr/bin/env Rscript

# Standard heterogeneity reporting for context-wise correlations:
# - Cochran's Q p-values
# - I^2
# - BH FDR across tested pairs (within an analysis domain)
#
# Outputs:
# - results/heterogeneity_summary.csv
# - doc/bioinformatics/figures/Figure_Heterogeneity_I2.pdf
#
# This script is designed to run from the repo root and uses the local path
# configuration in $HOME/Paths/context.json (via 00-paths.R).

source("00-paths.R")
source(file.path("R", "contextcorr.R"))

out_csv <- file.path("results", "heterogeneity_summary.csv")
out_pdf <- file.path("doc", "bioinformatics", "figures", "Figure_Heterogeneity_I2.pdf")
cache_rds <- file.path("results", "cache", "heterogeneity_I2_FDR.rds")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cache_rds), recursive = TRUE, showWarnings = FALSE)

f1 <- file.path(paths$mirtcga, "matchedNoGBM.Rda")
f2 <- file.path(paths$mirtcga, "matchedNoGBM.rda")
if (file.exists(f1)) {
  load(f1)
} else if (file.exists(f2)) {
  load(f2)
} else {
  stop("Cannot find matchedNoGBM.Rda or matchedNoGBM.rda in ", paths$mirtcga)
}

set.seed(1)

min_n_per_cohort <- 20
min_k <- 10
eps <- 0.05

# Domain: moderate-to-large (for BH and I^2 distributions), but tractable.
G <- 5000
M <- 200

cache_params <- list(
  G = G,
  M = M,
  min_n_per_cohort = min_n_per_cohort,
  min_k = min_k,
  eps = eps,
  matched_file = if (file.exists(f1)) normalizePath(f1) else normalizePath(f2)
)

cache <- NULL
if (file.exists(cache_rds)) {
  tmp <- readRDS(cache_rds)
  if (is.list(tmp) && identical(tmp$params, cache_params)) cache <- tmp
}

if (is.null(cache)) {
  gene_var <- apply(mRNAdata, 1, var)
  mir_var <- apply(miRdata, 1, var)

  genes <- names(sort(gene_var, decreasing = TRUE))[1:G]
  mirs <- names(sort(mir_var, decreasing = TRUE))[1:M]

  mRNA_sub <- mRNAdata[genes, , drop = FALSE]
  miR_sub <- miRdata[mirs, , drop = FALSE]

  rm(gene_var, mir_var)
  gc()

  cohorts <- levels(cancerType)
  idx_by_ct <- lapply(cohorts, function(ct) which(cancerType == ct))
  names(idx_by_ct) <- cohorts

  sum_w <- matrix(0, nrow = G, ncol = M, dimnames = list(genes, mirs))
  sum_wz <- matrix(0, nrow = G, ncol = M, dimnames = list(genes, mirs))
  sum_wz2 <- matrix(0, nrow = G, ncol = M, dimnames = list(genes, mirs))
  k_valid <- matrix(0L, nrow = G, ncol = M, dimnames = list(genes, mirs))

  pos_counts <- matrix(0L, nrow = G, ncol = M, dimnames = list(genes, mirs))
  neg_counts <- matrix(0L, nrow = G, ncol = M, dimnames = list(genes, mirs))

  for (ct in cohorts) {
    idx <- idx_by_ct[[ct]]
    n <- length(idx)
    if (n < min_n_per_cohort) next

    X <- mRNA_sub[, idx, drop = FALSE]
    Y <- miR_sub[, idx, drop = FALSE]

    scaled_mir <- scale(t(Y)) / (ncol(Y) - 1)
    std_rna <- t(scale(t(X)))
    r_ct <- std_rna %*% scaled_mir

    valid <- is.finite(r_ct) & abs(r_ct) < 1
    if (!any(valid)) next

    w <- n - 3
    z <- contextcorr_fisher_z(r_ct)

    k_valid[valid] <- k_valid[valid] + 1L
    sum_w[valid] <- sum_w[valid] + w
    sum_wz[valid] <- sum_wz[valid] + w * z[valid]
    sum_wz2[valid] <- sum_wz2[valid] + w * (z[valid] ^ 2)

    pos <- valid & (r_ct >= eps)
    neg <- valid & (r_ct <= -eps)
    pos_counts[pos] <- pos_counts[pos] + 1L
    neg_counts[neg] <- neg_counts[neg] + 1L

    rm(X, Y, scaled_mir, std_rna, r_ct, valid, z, pos, neg)
    gc()
  }

df <- pmax(k_valid - 1L, 0L)
Q <- sum_wz2 - (sum_wz ^ 2) / pmax(sum_w, 1e-12)
p_het <- rep(NA_real_, length(Q))
ok <- df >= 1L & is.finite(Q) & Q >= 0
p_het[ok] <- stats::pchisq(as.vector(Q)[ok], df = as.vector(df)[ok], lower.tail = FALSE)

I2 <- rep(NA_real_, length(Q))
I2[ok] <- pmax(0, (as.vector(Q)[ok] - as.vector(df)[ok]) / pmax(as.vector(Q)[ok], 1e-12))

zbar <- as.vector(sum_wz / pmax(sum_w, 1e-12))
r_fixed <- contextcorr_inv_fisher_z(zbar)

fdr <- rep(NA_real_, length(p_het))
fdr[ok] <- stats::p.adjust(p_het[ok], method = "BH")

posv <- as.vector(pos_counts)
negv <- as.vector(neg_counts)
kv <- as.vector(k_valid)
mixed <- (posv > 0L) & (negv > 0L)

summarize_group <- function(name, idx) {
  idx <- idx & ok
  n_test <- sum(idx)
  data.frame(
    group = name,
    n_pairs = sum(idx),
    n_tested = n_test,
    frac_fdr_lt_0p05 = if (n_test > 0) mean(fdr[idx] < 0.05, na.rm = TRUE) else NA_real_,
    median_I2 = if (n_test > 0) stats::median(I2[idx], na.rm = TRUE) else NA_real_,
    I2_q25 = if (n_test > 0) stats::quantile(I2[idx], 0.25, na.rm = TRUE) else NA_real_,
    I2_q75 = if (n_test > 0) stats::quantile(I2[idx], 0.75, na.rm = TRUE) else NA_real_,
    frac_I2_gt_0p75 = if (n_test > 0) mean(I2[idx] > 0.75, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}

idx_base <- kv >= min_k
idx_k20 <- kv >= 20

summary_df <- rbind(
  summarize_group("all (k>=10)", idx_base),
  summarize_group("mixed-sign (k>=10, eps=0.05)", idx_base & mixed),
  summarize_group("consistent-sign (k>=10, eps=0.05)", idx_base & !mixed),
  summarize_group("all (k>=20)", idx_k20),
  summarize_group("mixed-sign (k>=20, eps=0.05)", idx_k20 & mixed),
  summarize_group("consistent-sign (k>=20, eps=0.05)", idx_k20 & !mixed)
)

summary_df$min_n_per_cohort <- min_n_per_cohort
summary_df$min_k <- min_k
summary_df$eps <- eps
summary_df$G <- G
summary_df$M <- M
summary_df$domain_pairs <- G * M

utils::write.csv(summary_df, out_csv, row.names = FALSE)

# Figure: I^2 distribution + I^2 vs |r_fixed| (sampled for visibility)
set.seed(1)
idx_plot <- which(ok & idx_base)
if (length(idx_plot) > 60000) idx_plot <- sample(idx_plot, 60000)

pdf(out_pdf, width = 8.5, height = 4.0)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

hist(I2[idx_plot],
  breaks = 50,
  main = expression(paste("A  Heterogeneity (", I^2, ") distribution")),
  xlab = expression(I^2),
  col = "grey80",
  border = "white"
)

plot(abs(r_fixed[idx_plot]), I2[idx_plot],
  pch = 16,
  col = rgb(0, 0, 0, 0.07),
  xlab = expression(paste("|", r[fixed], "| (within-cohort pooled)")),
  ylab = expression(I^2),
  main = expression(paste("B  ", I^2, " vs effect size"))
)
abline(h = 0.75, col = "#d95f02", lwd = 2, lty = 2)

dev.off()

  saveRDS(
    list(
      params = cache_params,
      summary_df = summary_df,
      ok = ok,
      I2 = I2,
      r_fixed = r_fixed,
      k_valid = as.vector(k_valid),
      mixed = mixed
    ),
    cache_rds
  )
} else {
  summary_df <- cache$summary_df
  ok <- cache$ok
  I2 <- cache$I2
  r_fixed <- cache$r_fixed
  kv <- cache$k_valid
  mixed <- cache$mixed

  # Re-write CSV/figure from cache (idempotent).
  utils::write.csv(summary_df, out_csv, row.names = FALSE)

  set.seed(1)
  idx_base <- kv >= min_k
  idx_plot <- which(ok & idx_base)
  if (length(idx_plot) > 60000) idx_plot <- sample(idx_plot, 60000)

  pdf(out_pdf, width = 8.5, height = 4.0)
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  hist(I2[idx_plot],
    breaks = 50,
    main = expression(paste("A  Heterogeneity (", I^2, ") distribution")),
    xlab = expression(I^2),
    col = "grey80",
    border = "white"
  )
  plot(abs(r_fixed[idx_plot]), I2[idx_plot],
    pch = 16,
    col = rgb(0, 0, 0, 0.07),
    xlab = expression(paste("|", r[fixed], "| (within-cohort pooled)")),
    ylab = expression(I^2),
    main = expression(paste("B  ", I^2, " vs effect size"))
  )
  abline(h = 0.75, col = "#d95f02", lwd = 2, lty = 2)
  dev.off()
}

message("Wrote:")
message("- ", out_csv)
message("- ", out_pdf)
