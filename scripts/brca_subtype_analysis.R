#!/usr/bin/env Rscript

# BRCA intrinsic subtype analysis: mirrors robustness_signflip.R and
# heterogeneity_I2_FDR.R but uses BRCA molecular subtypes (PAM50 or IHC) as the
# context label instead of cancer type.
#
# Lynne Davidson predicts that HER2/HR+/TNBC will show the same pattern of
# pooled-vs-within-context sign disagreement that the paper documents across
# cancer types.  This script computes the headline statistics for both context
# definitions and writes a CSV summary plus per-pair caches for downstream use.
#
# Outputs:
# - results/brca_subtype_summary.csv    (aggregated headlines)
# - results/cache/brca_subtype_PAM50.rds
# - results/cache/brca_subtype_IHC.rds

suppressPackageStartupMessages({
  source("00-paths.R")
  source(file.path("R", "contextcorr.R"))
})

out_csv <- file.path("results", "brca_subtype_summary.csv")
cache_dir <- file.path("results", "cache")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

f1 <- file.path(paths$mirtcga, "matchedNoGBM.Rda")
f2 <- file.path(paths$mirtcga, "matchedNoGBM.rda")
matched_file <- if (file.exists(f1)) f1 else f2
load(matched_file)

clin_path <- file.path(paths$clean, "BRCA_clinicalMatrix.tsv")
clin <- read.delim(clin_path, stringsAsFactors = FALSE, check.names = FALSE)

# ---------------------------------------------------------------------------
# 1. Build BRCA tumor subset and attach subtype labels
# ---------------------------------------------------------------------------
brca_idx <- which(cancerType == "BRCA" & sampleType == "tumor")
brca_ids_dot <- colnames(mRNAdata)[brca_idx]
brca_ids15 <- gsub("\\.", "-", substr(brca_ids_dot, 1, 15))

m <- match(brca_ids15, clin$sampleID)
stopifnot(all(!is.na(m)))

pam50_raw <- clin$PAM50Call_RNAseq[m]
er_raw <- clin$ER_Status_nature2012[m]
pr_raw <- clin$PR_Status_nature2012[m]
her2_raw <- clin$HER2_Final_Status_nature2012[m]

# PAM50 cleanup: drop empty / NA, keep five named subtypes (Normal-like is small
# but we keep it as long as n >= 10).
pam50 <- ifelse(pam50_raw %in% c("Basal", "Her2", "LumA", "LumB", "Normal"),
                pam50_raw, NA_character_)
pam50_tab <- table(pam50, useNA = "no")
keep_pam50 <- names(pam50_tab)[pam50_tab >= 10]
pam50 <- ifelse(pam50 %in% keep_pam50, pam50, NA_character_)

# IHC cleanup: only Positive / Negative kept; build HR+ / HER2+ / TNBC labels.
er <- ifelse(er_raw %in% c("Positive", "Negative"), er_raw, NA_character_)
pr <- ifelse(pr_raw %in% c("Positive", "Negative"), pr_raw, NA_character_)
h2 <- ifelse(her2_raw %in% c("Positive", "Negative"), her2_raw, NA_character_)
ihc <- rep(NA_character_, length(er))
ok_ihc <- !is.na(er) & !is.na(pr) & !is.na(h2)
ihc[ok_ihc & h2 == "Positive"] <- "HER2+"
ihc[ok_ihc & h2 == "Negative" & (er == "Positive" | pr == "Positive")] <- "HR+"
ihc[ok_ihc & h2 == "Negative" & er == "Negative" & pr == "Negative"] <- "TNBC"

# ---------------------------------------------------------------------------
# 2. Subset expression matrices to BRCA tumors with a context label
# ---------------------------------------------------------------------------
mRNA_brca <- mRNAdata[, brca_idx, drop = FALSE]
miR_brca <- miRdata[, brca_idx, drop = FALSE]

# Restrict to top G genes / top M miRs by variance computed within BRCA tumors.
G <- 5000
M <- 200
gene_var <- apply(mRNA_brca, 1, var)
mir_var <- apply(miR_brca, 1, var)
genes <- names(sort(gene_var, decreasing = TRUE))[1:G]
mirs <- names(sort(mir_var, decreasing = TRUE))[1:M]
mRNA_sub <- mRNA_brca[genes, , drop = FALSE]
miR_sub <- miR_brca[mirs, , drop = FALSE]

rm(gene_var, mir_var, mRNA_brca, miR_brca)
gc()

# ---------------------------------------------------------------------------
# 3. Core stratified-correlation computation for a single context label
# ---------------------------------------------------------------------------
analyze_context <- function(label_vec, label_name,
                            min_n_per_ctx = 20, min_k = 3,
                            eps_grid = c(0, 0.05, 0.10),
                            t_grid = c(0, 0.05, 0.10, 0.20, 0.30)) {
  cache_rds <- file.path(cache_dir, sprintf("brca_subtype_%s.rds", label_name))
  cache_params <- list(
    label = label_name, G = G, M = M,
    min_n_per_ctx = min_n_per_ctx, min_k = min_k,
    eps_grid = eps_grid, t_grid = t_grid,
    matched_file = normalizePath(matched_file),
    n_total = sum(!is.na(label_vec))
  )
  if (file.exists(cache_rds)) {
    cached <- readRDS(cache_rds)
    if (is.list(cached) && identical(cached$params, cache_params)) {
      message("Using cache for ", label_name)
      return(cached)
    }
  }

  keep <- which(!is.na(label_vec))
  X <- mRNA_sub[, keep, drop = FALSE]
  Y <- miR_sub[, keep, drop = FALSE]
  ctx <- factor(label_vec[keep])

  # Pooled (global) correlation using BRCA tumors with this label.
  scaled_mir_g <- scale(t(Y)) / (ncol(Y) - 1)
  std_rna_g <- t(scale(t(X)))
  r_global <- std_rna_g %*% scaled_mir_g
  rm(scaled_mir_g, std_rna_g); gc()

  ctx_levels <- levels(ctx)

  # Per-context correlations stored as a 3-D array
  cor_array <- array(NA_real_, dim = c(G, M, length(ctx_levels)),
                     dimnames = list(genes, mirs, ctx_levels))
  n_by_ctx <- integer(length(ctx_levels))
  names(n_by_ctx) <- ctx_levels

  # Heterogeneity accumulators (Fisher-z meta-analysis)
  sum_w <- matrix(0, nrow = G, ncol = M)
  sum_wz <- matrix(0, nrow = G, ncol = M)
  sum_wz2 <- matrix(0, nrow = G, ncol = M)
  k_valid <- matrix(0L, nrow = G, ncol = M)

  # Sign counts per epsilon
  pos_counts <- lapply(eps_grid, function(e) matrix(0L, nrow = G, ncol = M))
  neg_counts <- lapply(eps_grid, function(e) matrix(0L, nrow = G, ncol = M))
  max_abs <- matrix(0, nrow = G, ncol = M)

  for (lev in ctx_levels) {
    idx <- which(ctx == lev)
    n_lev <- length(idx)
    n_by_ctx[lev] <- n_lev
    if (n_lev < min_n_per_ctx) next

    Xs <- X[, idx, drop = FALSE]
    Ys <- Y[, idx, drop = FALSE]
    sm <- scale(t(Ys)) / (ncol(Ys) - 1)
    sr <- t(scale(t(Xs)))
    r_ct <- sr %*% sm

    valid <- is.finite(r_ct) & abs(r_ct) < 1
    cor_array[, , lev][valid] <- r_ct[valid]

    z <- contextcorr_fisher_z(r_ct)
    w <- n_lev - 3
    k_valid[valid] <- k_valid[valid] + 1L
    sum_w[valid] <- sum_w[valid] + w
    sum_wz[valid] <- sum_wz[valid] + w * z[valid]
    sum_wz2[valid] <- sum_wz2[valid] + w * (z[valid]^2)
    max_abs[valid] <- pmax(max_abs[valid], abs(r_ct[valid]))

    for (i in seq_along(eps_grid)) {
      e <- eps_grid[i]
      pos <- valid & (r_ct >= e)
      neg <- valid & (r_ct <= -e)
      pos_counts[[i]][pos] <- pos_counts[[i]][pos] + 1L
      neg_counts[[i]][neg] <- neg_counts[[i]][neg] + 1L
    }

    rm(Xs, Ys, sm, sr, r_ct, z, valid, pos, neg)
    gc()
  }

  rm(X, Y); gc()

  # Heterogeneity Q / I^2 / fixed-effect r
  df_het <- pmax(k_valid - 1L, 0L)
  Q <- sum_wz2 - (sum_wz^2) / pmax(sum_w, 1e-12)
  ok_het <- df_het >= 1L & is.finite(Q) & Q >= 0
  p_het <- rep(NA_real_, length(Q))
  p_het[ok_het] <- stats::pchisq(as.vector(Q)[ok_het],
                                  df = as.vector(df_het)[ok_het],
                                  lower.tail = FALSE)
  I2 <- rep(NA_real_, length(Q))
  I2[ok_het] <- pmax(0,
    (as.vector(Q)[ok_het] - as.vector(df_het)[ok_het]) /
      pmax(as.vector(Q)[ok_het], 1e-12))
  zbar <- as.vector(sum_wz / pmax(sum_w, 1e-12))
  r_fixed <- contextcorr_inv_fisher_z(zbar)
  fdr <- rep(NA_real_, length(p_het))
  fdr[ok_het] <- stats::p.adjust(p_het[ok_het], method = "BH")

  # Sign-mixedness and Simpson reversal vs effect size
  rg_vec <- as.vector(r_global)
  df_rates <- data.frame()
  for (i in seq_along(eps_grid)) {
    e <- eps_grid[i]
    pos <- pos_counts[[i]]
    neg <- neg_counts[[i]]
    used <- (pos + neg) > 0L
    maj <- matrix(0L, nrow = G, ncol = M)
    maj[used] <- sign(pos[used] - neg[used])
    mixed <- (pos > 0L) & (neg > 0L)

    for (t in t_grid) {
      elig_mixed <- (k_valid >= min_k) & (max_abs >= t)
      frac_mixed <- if (sum(elig_mixed) > 0)
        sum(mixed & elig_mixed) / sum(elig_mixed) else NA_real_

      elig_simpson <- (k_valid >= min_k) & (abs(matrix(rg_vec, nrow = G)) >= t) &
                      (maj != 0L) & (sign(matrix(rg_vec, nrow = G)) != 0)
      simpson <- elig_simpson & (sign(matrix(rg_vec, nrow = G)) != maj)
      frac_simpson <- if (sum(elig_simpson) > 0)
        sum(simpson) / sum(elig_simpson) else NA_real_

      df_rates <- rbind(df_rates, data.frame(
        context = label_name, eps = e, threshold = t,
        eligible_mixed = sum(elig_mixed),
        frac_mixed = frac_mixed,
        eligible_simpson = sum(elig_simpson),
        frac_simpson = frac_simpson,
        stringsAsFactors = FALSE
      ))
    }
  }

  # I^2 summary (k>=min_k path)
  kv <- as.vector(k_valid)
  idx_base <- kv >= min_k & ok_het
  het_summary <- data.frame(
    context = label_name,
    n_samples = sum(!is.na(label_vec)),
    n_levels = length(ctx_levels),
    levels = paste(ctx_levels, collapse = "/"),
    n_pairs_total = G * M,
    n_pairs_tested = sum(idx_base),
    median_I2 = stats::median(I2[idx_base], na.rm = TRUE),
    I2_q25 = stats::quantile(I2[idx_base], 0.25, na.rm = TRUE),
    I2_q75 = stats::quantile(I2[idx_base], 0.75, na.rm = TRUE),
    frac_I2_gt_0p75 = mean(I2[idx_base] > 0.75, na.rm = TRUE),
    frac_fdr_lt_0p05 = mean(fdr[idx_base] < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # Simpson reversal table — flag candidate exemplars
  # For figure: pairs with |r_global| >= 0.2 AND maj sign disagrees AND large
  # |max within-context r|.
  pos05 <- pos_counts[[which(eps_grid == 0.05)]]
  neg05 <- neg_counts[[which(eps_grid == 0.05)]]
  used05 <- (pos05 + neg05) > 0L
  maj05 <- matrix(0L, nrow = G, ncol = M)
  maj05[used05] <- sign(pos05[used05] - neg05[used05])
  rg_mat <- matrix(rg_vec, nrow = G)
  reversal_mask <- (k_valid >= min_k) & (abs(rg_mat) >= 0.2) &
                   (maj05 != 0L) & (sign(rg_mat) != 0) & (sign(rg_mat) != maj05)

  out <- list(
    params = cache_params,
    label = label_name,
    ctx_levels = ctx_levels,
    n_by_ctx = n_by_ctx,
    cor_array = cor_array,
    r_global = r_global,
    r_fixed = matrix(r_fixed, nrow = G, dimnames = list(genes, mirs)),
    I2 = matrix(I2, nrow = G, dimnames = list(genes, mirs)),
    p_het = matrix(p_het, nrow = G, dimnames = list(genes, mirs)),
    fdr = matrix(fdr, nrow = G, dimnames = list(genes, mirs)),
    k_valid = k_valid,
    max_abs = max_abs,
    pos_counts = pos_counts,
    neg_counts = neg_counts,
    eps_grid = eps_grid,
    t_grid = t_grid,
    rates = df_rates,
    het_summary = het_summary,
    reversal_mask = reversal_mask
  )

  saveRDS(out, cache_rds)
  out
}

# ---------------------------------------------------------------------------
# 4. Run for both contexts and write the headline table
# ---------------------------------------------------------------------------
res_pam50 <- analyze_context(pam50, "PAM50")
res_ihc <- analyze_context(ihc, "IHC")

rates_combined <- rbind(res_pam50$rates, res_ihc$rates)
het_combined <- rbind(res_pam50$het_summary, res_ihc$het_summary)

utils::write.csv(rates_combined,
  file.path("results", "brca_subtype_rates.csv"), row.names = FALSE)
utils::write.csv(het_combined,
  file.path("results", "brca_subtype_heterogeneity.csv"), row.names = FALSE)

# Single combined headline table
headline <- data.frame(
  context = c("PAM50", "IHC"),
  n_samples = c(res_pam50$het_summary$n_samples, res_ihc$het_summary$n_samples),
  levels = c(res_pam50$het_summary$levels, res_ihc$het_summary$levels),
  median_I2 = c(res_pam50$het_summary$median_I2, res_ihc$het_summary$median_I2),
  frac_I2_gt_0p75 = c(res_pam50$het_summary$frac_I2_gt_0p75, res_ihc$het_summary$frac_I2_gt_0p75),
  frac_fdr_lt_0p05 = c(res_pam50$het_summary$frac_fdr_lt_0p05, res_ihc$het_summary$frac_fdr_lt_0p05),
  mixed_eps0p05_t0 = c(
    subset(res_pam50$rates, eps == 0.05 & threshold == 0)$frac_mixed,
    subset(res_ihc$rates,   eps == 0.05 & threshold == 0)$frac_mixed),
  mixed_eps0p05_t0p2 = c(
    subset(res_pam50$rates, eps == 0.05 & threshold == 0.2)$frac_mixed,
    subset(res_ihc$rates,   eps == 0.05 & threshold == 0.2)$frac_mixed),
  simpson_eps0p05_t0p2 = c(
    subset(res_pam50$rates, eps == 0.05 & threshold == 0.2)$frac_simpson,
    subset(res_ihc$rates,   eps == 0.05 & threshold == 0.2)$frac_simpson),
  stringsAsFactors = FALSE
)
utils::write.csv(headline, out_csv, row.names = FALSE)

message("Wrote: ", out_csv)
message("Wrote: results/brca_subtype_rates.csv")
message("Wrote: results/brca_subtype_heterogeneity.csv")
