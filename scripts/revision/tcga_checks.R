#!/usr/bin/env Rscript
# =============================================================================
# scripts/revision/tcga_checks.R
#
# Sample-type and Fisher-variance checks (TCGA miRNA-mRNA, 10^6 high-variance domain).
#
#   Part A   Sample type. The matched TCGA matrix mixes primary tumors (code 01), recurrent
#                  (02) and additional-primary (05) tumors, metastases (06, 07) and adjacent
#                  normal tissue (11). Recompute every Table 1 quantity on
#                    A1  all samples (reproduces the published values),
#                    A2  primary tumors only (code 01),
#                    A3  all tumor samples, adjacent normal excluded,
#                    A4  primary tumors only, high-variance domain re-selected on primary tumors,
#                  plus per-cohort sign stability of r_c between A1 and A2.
#   Part B   Fisher-z variance 1/(N_c - 3). Checks on the observed marginals:
#                    B1  zero inflation / tails of the analyzed features within cohorts,
#                    B2  within-cohort permutation null (rho = 0 with the observed marginals):
#                        calibration of Var(z_c), Cochran's Q and I^2 on all 10^6 pairs,
#                    B3  within-cohort nonparametric bootstrap of z_c at the observed correlation
#                        (random pair subset): Var_boot(z_c) * (N_c - 3), and Q / I^2 recomputed
#                        with bootstrap-variance weights,
#                    B4  rank-based (Spearman) within-cohort heterogeneity with the Fieller-Hartley-
#                        Pearson variance 1.06 / (N_c - 3).
#
# Metric engine: compute_metrics() and helpers are taken verbatim from
# scripts/revision/uniform_table.R (parsed and evaluated, not copied), so definitions match Table 1.
# Run from the repo root:  Rscript scripts/revision/tcga_checks.R
# Outputs: results/revision/sampletype_fisher_checks/
# =============================================================================

t_start <- Sys.time()
suppressPackageStartupMessages({
  source("00-paths.R")
  source(file.path("R", "contextcorr.R"))
  library(parallel)
})

N_CORES <- max(1L, min(as.integer(Sys.getenv("CHECKS_CORES", "16")), parallel::detectCores()))
EPS <- 0.05; THR <- 0.2; K_MIN <- 2L; BOOT_CHUNKS <- 20L
B_BOOT <- as.integer(Sys.getenv("CHECKS_B", "1000"))
SEED <- 20261111L
out_dir <- file.path("results", "revision", "sampletype_fisher_checks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- pull the function definitions from uniform_table.R ----------------------
ex <- parse(file.path("scripts", "revision", "uniform_table.R"))
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))) eval(e, envir = globalenv())
}
stopifnot(exists("compute_metrics"), exists("cor_rows_product"), exists("residualize_by_group"),
          exists("rowcor_pairs"), exists("p_from_r"))

# --- data ---------------------------------------------------------------------
f_matched <- file.path(paths$mirtcga, "matchedNoGBM.rda")
load(f_matched)  # mRNAdata, miRdata, cancerType, sampleType
code <- substr(colnames(mRNAdata), 14, 15)
stopifnot(all(code %in% c("01", "02", "05", "06", "07", "11")))
cancerType <- as.character(cancerType)

comp <- as.data.frame.matrix(table(cancerType, code))
comp$cohort <- rownames(comp)
comp$n_all <- rowSums(comp[, c("01", "02", "05", "06", "07", "11")])
comp <- comp[, c("cohort", "n_all", "01", "02", "05", "06", "07", "11")]
names(comp) <- c("cohort", "n_all", "n_primary_01", "n_recurrent_02", "n_additional_primary_05",
                 "n_metastatic_06", "n_metastatic_07", "n_adjacent_normal_11")
write.csv(comp, file.path(out_dir, "tcga_sampletype_by_cohort.csv"), row.names = FALSE)

G <- 5000; M <- 200
top_by_var <- function(A, n) names(sort(apply(A, 1, var), decreasing = TRUE))[seq_len(n)]
genes_all <- top_by_var(mRNAdata, G); mirs_all <- top_by_var(miRdata, M)

run_variant <- function(keep, genes, mirs, label, boot_seed, return_R = FALSE) {
  t0 <- Sys.time()
  X <- mRNAdata[genes, keep, drop = FALSE]; Y <- miRdata[mirs, keep, drop = FALSE]
  ct <- cancerType[keep]
  tab <- table(ct)
  retained <- names(tab)[tab >= 20]
  pop <- which(ct %in% retained)
  ctx <- droplevels(factor(ct[pop], levels = retained))
  rg <- as.vector(cor_rows_product(X[, pop], Y[, pop]))
  R <- matrix(NA_real_, nrow = G * M, ncol = length(retained), dimnames = list(NULL, retained))
  for (cc in retained) {
    idx <- pop[ctx == cc]
    R[, cc] <- as.vector(cor_rows_product(X[, idx, drop = FALSE], Y[, idx, drop = FALSE]))
  }
  rr <- as.vector(cor_rows_product(residualize_by_group(X[, pop], ctx), residualize_by_group(Y[, pop], ctx)))
  n_ctx <- as.integer(tab[retained]); names(n_ctx) <- retained
  res <- compute_metrics(rg, rr, R, n_ctx, boot_B = B_BOOT, boot_seed = boot_seed)
  message(sprintf("[%s] %d samples, %d cohorts, reversal %.2f%%, %.1f s", label, length(pop), length(retained),
                  res$metrics[["m5_pct_simpson_reversal_ge_0.2"]], as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out <- list(label = label, n_samples = length(pop), retained = retained, dropped = setdiff(names(tab), retained),
              n_ctx = n_ctx, res = res, rg = rg, rr = rr)
  if (return_R) out$R <- R
  out
}

# =============================================================================
# Part A: sample type
# =============================================================================
A1 <- run_variant(seq_len(ncol(mRNAdata)), genes_all, mirs_all, "all_samples", SEED + 1L, return_R = TRUE)
A2 <- run_variant(which(code == "01"), genes_all, mirs_all, "primary_tumor_only", SEED + 11L, return_R = TRUE)
A3 <- run_variant(which(code != "11"), genes_all, mirs_all, "tumor_samples_no_adjacent_normal", SEED + 21L)
genes_p <- top_by_var(mRNAdata[, code == "01"], G); mirs_p <- top_by_var(miRdata[, code == "01"], M)
A4 <- run_variant(which(code == "01"), genes_p, mirs_p, "primary_tumor_only_domain_reselected", SEED + 31L)
domain_overlap <- c(genes = length(intersect(genes_all, genes_p)), mirs = length(intersect(mirs_all, mirs_p)))

variants <- list(A1, A2, A3, A4)
met <- do.call(rbind, lapply(variants, function(v) {
  m <- v$res$metrics
  ci <- v$res$ci
  data.frame(variant = v$label, n_samples = v$n_samples, n_cohorts = length(v$retained),
             cohorts_dropped = paste(v$dropped, collapse = ";"),
             t(m), check.names = FALSE,
             simpson_ci_lo = ci$ci95_lo_pct[ci$rate == "simpson_reversal_ge_0.2"],
             simpson_ci_hi = ci$ci95_hi_pct[ci$rate == "simpson_reversal_ge_0.2"],
             resid_flip_all_ci_lo = ci$ci95_lo_pct[ci$rate == "resid_flip_all"],
             resid_flip_all_ci_hi = ci$ci95_hi_pct[ci$rate == "resid_flip_all"],
             nominal_flip_all_ci_lo = ci$ci95_lo_pct[ci$rate == "nominal_flip_all"],
             nominal_flip_all_ci_hi = ci$ci95_hi_pct[ci$rate == "nominal_flip_all"],
             stringsAsFactors = FALSE)
}))
met$domain_genes_shared_with_main <- c(G, G, G, domain_overlap[["genes"]])
met$domain_mirs_shared_with_main <- c(M, M, M, domain_overlap[["mirs"]])
write.csv(met, file.path(out_dir, "tcga_sampletype_metrics.csv"), row.names = FALSE)

# Pair-level concordance of the Simpson classification between all-sample and primary-only (same pairs).
p1 <- A1$res$pairs; p2 <- A2$res$pairs
# compute_metrics() drops non-analyzed pairs; recover pair indices.
idx_analyzed <- function(v) {
  valid <- is.finite(v$R) & abs(v$R) < 1
  which(is.finite(v$rg) & is.finite(v$rr) & rowSums(valid) >= K_MIN)  # same rule as compute_metrics()
}
i1 <- idx_analyzed(A1); i2 <- idx_analyzed(A2)
stopifnot(length(i1) == length(p1$simpson), length(i2) == length(p2$simpson))
rev1 <- rep(NA, G * M); rev2 <- rep(NA, G * M)
rev1[i1] <- ifelse(p1$elig_simpson_strong, p1$simpson, NA)
rev2[i2] <- ifelse(p2$elig_simpson_strong, p2$simpson, NA)
both <- !is.na(rev1) & !is.na(rev2)
pair_conc <- data.frame(
  n_pairs_eligible_in_both = sum(both),
  n_reversed_all = sum(rev1[both]), n_reversed_primary = sum(rev2[both]),
  n_reversed_in_both = sum(rev1[both] & rev2[both]),
  pct_of_all_sample_reversals_still_reversed = 100 * sum(rev1[both] & rev2[both]) / sum(rev1[both]),
  pct_classification_agreement = 100 * mean(rev1[both] == rev2[both]),
  reversal_rate_all_on_common = 100 * mean(rev1[both]), reversal_rate_primary_on_common = 100 * mean(rev2[both]))
write.csv(pair_conc, file.path(out_dir, "tcga_sampletype_pair_concordance.csv"), row.names = FALSE)

# Per-cohort sign stability of r_c (eps = 0.05 sign) between all samples and primary tumors only.
sgn <- function(r) ifelse(is.na(r) | !is.finite(r) | abs(r) >= 1, NA, ifelse(r >= EPS, 1L, ifelse(r <= -EPS, -1L, 0L)))
common <- intersect(A1$retained, A2$retained)
coh <- do.call(rbind, lapply(common, function(cc) {
  a <- A1$R[, cc]; b <- A2$R[, cc]
  sa <- sgn(a); sb <- sgn(b); ok <- !is.na(sa) & !is.na(sb)
  data.frame(cohort = cc, n_all = A1$n_ctx[[cc]], n_primary = A2$n_ctx[[cc]],
             n_removed = A1$n_ctx[[cc]] - A2$n_ctx[[cc]],
             pct_sign_category_changed = 100 * mean(sa[ok] != sb[ok]),
             pct_opposite_sign = 100 * mean(sa[ok] * sb[ok] == -1),
             pct_opposite_sign_among_abs_rc_ge_0.2 = { s <- ok & abs(a) >= THR; if (any(s)) 100 * mean(sa[s] * sb[s] == -1) else NA_real_ },
             cor_rc_all_vs_primary = if (all(a[ok] == b[ok])) 1 else stats::cor(a[ok], b[ok]),
             median_abs_diff_rc = stats::median(abs(a[ok] - b[ok])),
             stringsAsFactors = FALSE)
}))
coh <- merge(coh, comp[, c("cohort", "n_adjacent_normal_11", "n_metastatic_06", "n_metastatic_07", "n_recurrent_02", "n_additional_primary_05")], by = "cohort")
coh <- coh[order(-coh$n_removed), ]
write.csv(coh, file.path(out_dir, "tcga_sampletype_rc_stability_by_cohort.csv"), row.names = FALSE)
A2$R <- NULL; gc(verbose = FALSE)
message("Part A done")

# =============================================================================
# Part B: Fisher-z variance checks on the main (all-sample) analysis
# =============================================================================
X <- mRNAdata[genes_all, , drop = FALSE]; Y <- miRdata[mirs_all, , drop = FALSE]
cohorts <- A1$retained
n_ctx <- A1$n_ctx
stopifnot(identical(sort(cohorts), sort(unique(cancerType))))
mir_floor <- 0; mrna_floor <- log2(10)
stopifnot(min(Y) >= mir_floor - 1e-12, min(X) >= mrna_floor - 1e-9)

# ---- B1: marginals within cohorts -------------------------------------------
exkurt <- function(A) {  # row-wise excess kurtosis
  m <- rowMeans(A); C <- A - m; v <- rowMeans(C^2)
  ifelse(v > 0, rowMeans(C^4) / v^2 - 3, NA_real_)
}
zero_mir <- sapply(cohorts, function(cc) rowMeans(Y[, cancerType == cc, drop = FALSE] <= mir_floor + 1e-12))  # M x k
zero_mrna <- sapply(cohorts, function(cc) rowMeans(X[, cancerType == cc, drop = FALSE] <= mrna_floor + 1e-9))  # G x k
kurt_mir <- sapply(cohorts, function(cc) exkurt(Y[, cancerType == cc, drop = FALSE]))
kurt_mrna <- sapply(cohorts, function(cc) exkurt(X[, cancerType == cc, drop = FALSE]))
q <- function(x) stats::quantile(x, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE, names = FALSE)
marg <- rbind(
  data.frame(quantity = "miRNA zero fraction within cohort (feature-cohort cells)", t(q(zero_mir))),
  data.frame(quantity = "mRNA zero-count fraction within cohort (feature-cohort cells)", t(q(zero_mrna))),
  data.frame(quantity = "miRNA excess kurtosis within cohort (nonzero-variance cells)", t(q(kurt_mir))),
  data.frame(quantity = "mRNA excess kurtosis within cohort (nonzero-variance cells)", t(q(kurt_mrna)))
)
names(marg)[2:6] <- c("p10", "p25", "p50", "p75", "p90")
marg$pct_cells_zero_fraction_ge_0.5 <- c(100 * mean(zero_mir >= 0.5), 100 * mean(zero_mrna >= 0.5), NA, NA)
marg$pct_cells_zero_fraction_ge_0.2 <- c(100 * mean(zero_mir >= 0.2), 100 * mean(zero_mrna >= 0.2), NA, NA)
marg$pct_cells_abs_excess_kurtosis_gt_2 <- c(NA, NA, 100 * mean(abs(kurt_mir) > 2, na.rm = TRUE), 100 * mean(abs(kurt_mrna) > 2, na.rm = TRUE))
write.csv(marg, file.path(out_dir, "fisher_B1_marginals.csv"), row.names = FALSE)

# Zero-fraction stratum of each pair-cohort cell = max of the two features' zero fractions.
pair_g <- rep(seq_len(G), times = M); pair_m <- rep(seq_len(M), each = G)  # as.vector(G x M) order
zf_cell <- function(cj) pmax(zero_mrna[pair_g, cj], zero_mir[pair_m, cj])
zbreaks <- c(-Inf, 0.05, 0.2, 0.5, Inf); zlabels <- c("<5%", "5-20%", "20-50%", ">=50%")

# ---- B2: within-cohort permutation null on all 10^6 pairs ---------------------
het_stats <- function(R, n_ctx, w_mat = NULL) {
  valid <- is.finite(R) & abs(R) < 1
  Rv <- R; Rv[!valid] <- 0
  z <- contextcorr_fisher_z(Rv); z[!valid] <- 0
  W <- if (is.null(w_mat)) matrix(rep(as.numeric(n_ctx) - 3, each = nrow(R)), nrow = nrow(R)) else w_mat
  W[!valid | !is.finite(W)] <- 0
  k <- rowSums(W > 0)
  sw <- rowSums(W); swz <- rowSums(W * z); swz2 <- rowSums(W * z^2)
  Q <- pmax(0, swz2 - swz^2 / sw); df <- k - 1
  ok <- k >= K_MIN
  p <- stats::pchisq(Q, df = df, lower.tail = FALSE)
  I2 <- pmax(0, (Q - df) / pmax(Q, 1e-12))
  list(Q = Q[ok], df = df[ok], p = p[ok], I2 = I2[ok], ok = ok)
}
het_summary <- function(h, label) {
  data.frame(analysis = label, n_pairs = length(h$I2),
             pct_Q_p_lt_0.05 = 100 * mean(h$p < 0.05),
             pct_Q_BH_lt_0.05 = 100 * mean(stats::p.adjust(h$p, "BH") < 0.05),
             n_Q_BH_lt_0.05 = sum(stats::p.adjust(h$p, "BH") < 0.05),
             min_Q_p = min(h$p),
             median_I2 = stats::median(h$I2), I2_q25 = stats::quantile(h$I2, 0.25, names = FALSE),
             I2_q75 = stats::quantile(h$I2, 0.75, names = FALSE), pct_I2_gt_0.75 = 100 * mean(h$I2 > 0.75),
             mean_Q_over_df = mean(h$Q / h$df), stringsAsFactors = FALSE)
}

N_PERM <- as.integer(Sys.getenv("CHECKS_NPERM", "5"))
perm_rows <- list(); perm_var <- list()
for (b in seq_len(N_PERM)) {
  set.seed(SEED + 100L * b)
  Rp <- matrix(NA_real_, nrow = G * M, ncol = length(cohorts), dimnames = list(NULL, cohorts))
  for (cc in cohorts) {
    idx <- which(cancerType == cc)
    Yc <- Y[, idx, drop = FALSE]
    Yc <- t(apply(Yc, 1, function(v) v[sample.int(length(v))]))  # independent permutation per miRNA
    Rp[, cc] <- as.vector(cor_rows_product(X[, idx, drop = FALSE], Yc))
  }
  h <- het_stats(Rp, n_ctx)
  perm_rows[[b]] <- het_summary(h, sprintf("permutation_%d", b))
  # Var(z_c) * (N_c - 3) under the permutation null, by cohort and by zero-fraction stratum
  for (j in seq_along(cohorts)) {
    r <- Rp[, j]; ok <- is.finite(r) & abs(r) < 1
    z <- contextcorr_fisher_z(r[ok])
    strat <- cut(zf_cell(j)[ok], zbreaks, labels = zlabels, right = FALSE)
    # Normalized second moment (N_c - 3) E[z^2] across pairs, and the centered version (N_c - 3) Var(z)
    # across pairs. Under permutation E[r] = 0 exactly but E[atanh(r)] need not be 0 for skewed margins,
    # so both are reported. These pool different pairs (one permutation each per replicate); they are
    # aggregate summaries, not per-pair conditional variances.
    vs <- tapply(z, strat, function(x) mean(x^2) * (n_ctx[[j]] - 3))
    vc <- tapply(z, strat, function(x) stats::var(x) * (n_ctx[[j]] - 3))
    ms <- tapply(z, strat, mean)
    ns <- tapply(z, strat, length)
    perm_var[[length(perm_var) + 1]] <- data.frame(perm = b, cohort = cohorts[j], n = n_ctx[[j]],
      stratum = c("all", names(vs)), n_cells = c(sum(ok), as.integer(ns)),
      var_z_times_n_minus_3 = c(mean(z^2) * (n_ctx[[j]] - 3), as.numeric(vs)),
      centered_var_z_times_n_minus_3 = c(stats::var(z) * (n_ctx[[j]] - 3), as.numeric(vc)),
      mean_z_times_sqrt_n_minus_3 = c(mean(z) * sqrt(n_ctx[[j]] - 3), as.numeric(ms) * sqrt(n_ctx[[j]] - 3)),
      tail_pct_abs_z_gt_1.96_sd = c(100 * mean(abs(z) * sqrt(n_ctx[[j]] - 3) > 1.96),
                                    as.numeric(tapply(abs(z) * sqrt(n_ctx[[j]] - 3) > 1.96, strat, mean)) * 100),
      stringsAsFactors = FALSE)
  }
  rm(Rp); gc(verbose = FALSE)
  message(sprintf("permutation %d done", b))
}
perm_df <- do.call(rbind, perm_rows)
perm_var_df <- do.call(rbind, perm_var)
write.csv(perm_df, file.path(out_dir, "fisher_B2_permutation_heterogeneity.csv"), row.names = FALSE)
write.csv(perm_var_df, file.path(out_dir, "fisher_B2_permutation_var_by_cohort.csv"), row.names = FALSE)
pv <- perm_var_df[perm_var_df$stratum != "all", ]
pv_strat <- do.call(rbind, lapply(split(pv, pv$stratum), function(d) {
  d <- d[is.finite(d$var_z_times_n_minus_3) & d$n_cells > 0, ]
  data.frame(stratum = d$stratum[1], n_pair_cohort_cells = sum(d$n_cells) / N_PERM,
             var_ratio_weighted = sum(d$var_z_times_n_minus_3 * d$n_cells) / sum(d$n_cells),
             centered_var_ratio_weighted = sum(d$centered_var_z_times_n_minus_3 * d$n_cells) / sum(d$n_cells),
             max_abs_mean_z_scaled = max(abs(d$mean_z_times_sqrt_n_minus_3)),
             tail_pct_weighted = sum(d$tail_pct_abs_z_gt_1.96_sd * d$n_cells) / sum(d$n_cells), stringsAsFactors = FALSE)
}))
pv_all <- perm_var_df[perm_var_df$stratum == "all", ]
pv_strat <- rbind(pv_strat, data.frame(stratum = "all", n_pair_cohort_cells = sum(pv_all$n_cells) / N_PERM,
  var_ratio_weighted = sum(pv_all$var_z_times_n_minus_3 * pv_all$n_cells) / sum(pv_all$n_cells),
  centered_var_ratio_weighted = sum(pv_all$centered_var_z_times_n_minus_3 * pv_all$n_cells) / sum(pv_all$n_cells),
  max_abs_mean_z_scaled = max(abs(pv_all$mean_z_times_sqrt_n_minus_3)),
  tail_pct_weighted = sum(pv_all$tail_pct_abs_z_gt_1.96_sd * pv_all$n_cells) / sum(pv_all$n_cells)))
write.csv(pv_strat, file.path(out_dir, "fisher_B2_permutation_var_by_zero_stratum.csv"), row.names = FALSE)
message("B2 done")

# ---- B3: within-cohort bootstrap of z_c at the observed correlation ----------------
N_BOOT_PAIRS <- as.integer(Sys.getenv("CHECKS_BOOT_PAIRS", "5000"))
B_Z <- as.integer(Sys.getenv("CHECKS_BZ", "400"))
set.seed(SEED + 500L)
pp <- sort(sample.int(G * M, N_BOOT_PAIRS))
ig <- pair_g[pp]; im <- pair_m[pp]
R_obs <- A1$R[pp, , drop = FALSE]
boot_one_cohort <- function(j) {
  idx <- which(cancerType == cohorts[j]); n <- length(idx)
  Xc <- X[, idx, drop = FALSE]; Yc <- Y[, idx, drop = FALSE]
  set.seed(SEED + 1000L + j)
  Z <- matrix(NA_real_, nrow = N_BOOT_PAIRS, ncol = B_Z)
  XY <- rbind(Xc[ig, , drop = FALSE], Yc[im, , drop = FALSE])
  n_clip <- 0
  for (b in seq_len(B_Z)) {
    s <- sample.int(n, n, replace = TRUE)
    r <- rowcor_pairs(XY[, s, drop = FALSE], XY[, s, drop = FALSE], seq_len(N_BOOT_PAIRS), N_BOOT_PAIRS + seq_len(N_BOOT_PAIRS))
    r[abs(r) >= 1] <- NA
    n_clip <- n_clip + sum(abs(r) > 0.999999, na.rm = TRUE)  # values clipped inside contextcorr_fisher_z()
    Z[, b] <- contextcorr_fisher_z(r)
  }
  nfin <- rowSums(is.finite(Z))
  v <- apply(Z, 1, stats::var, na.rm = TRUE)
  v[nfin < 0.9 * B_Z] <- NA  # set aside when more than 10% of resamples give an undefined correlation
  list(v = v, n_failed_resamples = sum(B_Z - nfin), n_cells_set_aside = sum(nfin < 0.9 * B_Z), n_clip = n_clip)
}
boot_out <- parallel::mclapply(seq_along(cohorts), boot_one_cohort, mc.cores = min(length(cohorts), N_CORES))
Vb <- do.call(cbind, lapply(boot_out, `[[`, "v"))
colnames(Vb) <- cohorts
boot_fail <- data.frame(cohort = cohorts,
  n_failed_resamples = sapply(boot_out, `[[`, "n_failed_resamples"),
  n_pair_cells_set_aside = sapply(boot_out, `[[`, "n_cells_set_aside"),
  n_resample_correlations_clipped = sapply(boot_out, `[[`, "n_clip"))
write.csv(boot_fail, file.path(out_dir, "fisher_B3_bootstrap_failures.csv"), row.names = FALSE)
ratio <- sweep(Vb, 2, as.numeric(n_ctx) - 3, `*`)          # Var_boot(z) / (1 / (N_c - 3))
valid_obs <- is.finite(R_obs) & abs(R_obs) < 1
ratio[!valid_obs] <- NA
zf_sub <- sapply(seq_along(cohorts), function(j) pmax(zero_mrna[ig, j], zero_mir[im, j]))
rc_abs <- abs(R_obs)
strat_z <- cut(zf_sub, zbreaks, labels = zlabels, right = FALSE)
strat_r <- cut(rc_abs, c(-Inf, 0.1, 0.3, 0.5, Inf), labels = c("|r_c|<0.1", "0.1-0.3", "0.3-0.5", ">=0.5"), right = FALSE)
rs <- function(x, lab, by) {
  do.call(rbind, lapply(split(x, by), function(v) { v <- v[is.finite(v)]
    data.frame(stratifier = lab, stratum = NA, n_cells = length(v), median_ratio = stats::median(v),
               q25 = stats::quantile(v, 0.25, names = FALSE), q75 = stats::quantile(v, 0.75, names = FALSE),
               mean_ratio = mean(v), pct_ratio_gt_1.5 = 100 * mean(v > 1.5), pct_ratio_gt_2 = 100 * mean(v > 2)) }))
}
ratio_tab <- rbind(
  { d <- rs(as.vector(ratio), "all", rep("all", length(ratio))); d$stratum <- "all"; d },
  { d <- rs(as.vector(ratio), "zero fraction (max of pair)", as.vector(strat_z)); d$stratum <- rownames(d); d },
  { d <- rs(as.vector(ratio), "|r_c|", as.vector(strat_r)); d$stratum <- rownames(d); d })
rownames(ratio_tab) <- NULL
write.csv(ratio_tab, file.path(out_dir, "fisher_B3_bootstrap_variance_ratio.csv"), row.names = FALSE)
by_coh <- data.frame(cohort = cohorts, n = as.integer(n_ctx), median_ratio = apply(ratio, 2, stats::median, na.rm = TRUE),
                     q25 = apply(ratio, 2, stats::quantile, 0.25, na.rm = TRUE), q75 = apply(ratio, 2, stats::quantile, 0.75, na.rm = TRUE))
write.csv(by_coh, file.path(out_dir, "fisher_B3_bootstrap_variance_ratio_by_cohort.csv"), row.names = FALSE)

# One common mask: a pair-cohort cell enters all three weightings only if its observed r_c is valid
# AND its bootstrap variance is defined, so the three analyses use identical pair-context sets.
R_cmp <- R_obs
R_cmp[!is.finite(Vb)] <- NA
h_fisher <- het_stats(R_cmp, n_ctx)
Wb <- 1 / Vb; Wb[!is.finite(Wb)] <- 0
h_boot <- het_stats(R_cmp, n_ctx, w_mat = Wb)
stopifnot(identical(h_fisher$ok, h_boot$ok), all(h_fisher$df == h_boot$df))
# Third variant: Fisher weights scaled by each cohort's median bootstrap ratio (smooth, low-noise weights).
Ws <- matrix(rep((as.numeric(n_ctx) - 3) / by_coh$median_ratio, each = N_BOOT_PAIRS), nrow = N_BOOT_PAIRS)
h_scaled <- het_stats(R_cmp, n_ctx, w_mat = Ws)
stopifnot(identical(h_fisher$ok, h_scaled$ok))
b3 <- rbind(het_summary(h_fisher, "Fisher weights N_c - 3 (as in the paper), same pairs"),
            het_summary(h_boot, "pair-cohort bootstrap-variance weights 1 / Var_boot(z_c)"),
            het_summary(h_scaled, "Fisher weights divided by cohort median bootstrap ratio"))
b3$n_bootstrap_pairs <- N_BOOT_PAIRS; b3$B <- B_Z
write.csv(b3, file.path(out_dir, "fisher_B3_heterogeneity_bootstrap_weights.csv"), row.names = FALSE)
message("B3 done")

# ---- B4: rank-based heterogeneity on the full domain ---------------------------
rank_rows <- function(A) t(apply(A, 1, rank, ties.method = "average"))
Rs <- matrix(NA_real_, nrow = G * M, ncol = length(cohorts), dimnames = list(NULL, cohorts))
Xrk <- X; Yrk <- Y
for (cc in cohorts) {
  idx <- which(cancerType == cc)
  Xrk[, idx] <- rank_rows(X[, idx, drop = FALSE]); Yrk[, idx] <- rank_rows(Y[, idx, drop = FALSE])
  Rs[, cc] <- as.vector(cor_rows_product(Xrk[, idx, drop = FALSE], Yrk[, idx, drop = FALSE]))
}
rho_global <- as.vector(cor_rows_product(rank_rows(X), rank_rows(Y)))
ctx_all <- factor(cancerType, levels = cohorts)
rho_resid <- as.vector(cor_rows_product(residualize_by_group(Xrk, ctx_all), residualize_by_group(Yrk, ctx_all)))
# Spearman Fisher-z variance 1.06 / (N - 3): pass effective sizes so that N_eff - 3 = (N - 3) / 1.06.
n_eff <- (as.numeric(n_ctx) - 3) / 1.06 + 3; names(n_eff) <- cohorts
res_sp <- compute_metrics(rho_global, rho_resid, Rs, n_eff, boot_B = B_BOOT, boot_seed = SEED + 41L)
res_sp_fisherw <- compute_metrics(rho_global, rho_resid, Rs, as.numeric(n_ctx), boot_B = 0L)
keep_cols <- c("n_pairs_analyzed", "m2_pct_mixed_sign_eps0.05", "m2_pct_mixed_sign_eps0.05_maxabs_rc_ge_0.2",
               "m3_median_I2", "m3_I2_q25", "m3_I2_q75", "m3_pct_Q_fdr_lt_0.05", "m4_median_abs_r_RE",
               "m4_pct_sign_rRE_ne_rglobal_ge_0.2", "m5_pct_simpson_reversal_ge_0.2", "n_pairs_eligible_simpson_ge_0.2",
               "m5_pct_resid_flip_all")
b4 <- rbind(
  data.frame(analysis = "Pearson, Fisher weights N_c - 3 (Table 1)", t(A1$res$metrics[keep_cols]), check.names = FALSE),
  data.frame(analysis = "Spearman, variance 1.06 / (N_c - 3)", t(res_sp$metrics[keep_cols]), check.names = FALSE),
  data.frame(analysis = "Spearman, variance 1 / (N_c - 3)", t(res_sp_fisherw$metrics[keep_cols]), check.names = FALSE))
b4$simpson_ci_lo <- c(A1$res$ci$ci95_lo_pct[1], res_sp$ci$ci95_lo_pct[1], NA)
b4$simpson_ci_hi <- c(A1$res$ci$ci95_hi_pct[1], res_sp$ci$ci95_hi_pct[1], NA)
b4$note <- "m5_pct_resid_flip_all for Spearman rows uses cohort-centered within-cohort ranks"
write.csv(b4, file.path(out_dir, "fisher_B4_spearman_heterogeneity.csv"), row.names = FALSE)
message("B4 done")

writeLines(c(sprintf("total seconds: %.1f", as.numeric(difftime(Sys.time(), t_start, units = "secs"))),
             sprintf("B_BOOT = %d, N_PERM = %d, N_BOOT_PAIRS = %d, B_Z = %d, SEED = %d", B_BOOT, N_PERM, N_BOOT_PAIRS, B_Z, SEED),
             capture.output(sessionInfo())), file.path(out_dir, "run_info.txt"))
message("done")
