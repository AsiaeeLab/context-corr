#!/usr/bin/env Rscript
# =============================================================================
# scripts/revision/uniform_table.R
#
# Report five correlation summaries for every setting, a preprocessing table,
# and pair-level bootstrap CIs.
#
# Settings (context in parentheses):
#   1. TCGA miRNA-mRNA (cancer cohort)          2. GTEx mRNA-mRNA (SMTSD tissue)
#   3. PBMC 3k scRNA gene-gene (coarse cell type) 4. CITE-seq protein-RNA (joint WNN cluster)
#   5. TCGA BRCA miRNA-mRNA (PAM50 subtype)     6. TCGA BRCA miRNA-mRNA (IHC receptor group)
#
# Uniform conventions (identical for all six settings):
#   * Analysis population = samples/cells whose context has at least n_min members
#     (n_min is the per-setting minimum used in the paper: 20 for bulk, 200 cells
#     for single-cell). r_global, r_c and r_resid are ALL computed on this population,
#     so the between-vs-within diagnostic compares like with like.
#   * r_c valid if finite and |r_c| < 1; pairs need >= 2 valid contexts.
#   * eps = 0.05 sign tolerance: s_c = +1 (r_c >= eps), -1 (r_c <= -eps), else 0.
#   * mixed-sign = both +1 and -1 present. Simpson reversal = sign(r_global) != sign(sum_c s_c),
#     ties and r_global == 0 excluded, evaluated among |r_global| >= 0.2.
#   * Fisher z with weights N_c - 3; Cochran Q; I2 = max(0, (Q - (k-1))/Q); BH FDR within setting.
#   * r_RE = DerSimonian-Laird random-effects pooled correlation (same algebra as
#     random_effects_pool() in R/contextcorr.R, vectorized over pairs).
#   * r_resid = Pearson correlation after subtracting context means from both variables.
#   * Pair-level percentile bootstrap (B = 1000) for the reversal rate at |r_global| >= 0.2
#     and the residualized flip rates.
#
# Alternative population conventions are also computed for sensitivity checks
# in reproduction_checks.csv; the main reported results are in table1_metrics.csv.
#
# Run from the repo root:   Rscript scripts/revision/uniform_table.R
# Optional env var UNIFORM_TABLE_B overrides the bootstrap B (for quick tests).
# Outputs: results/revision/uniform_table/. Nothing else is written.
# =============================================================================

t_start <- Sys.time()
timings <- list()
tic <- function() Sys.time()
toc <- function(label, t0) {
  timings[[label]] <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("[%s] %.1f s", label, timings[[label]]))
}

suppressPackageStartupMessages({
  source("00-paths.R")
  source(file.path("R", "contextcorr.R"))
  library(Matrix)
  library(parallel)
  library(SeuratObject)  # only to read the two single-cell .rds inputs (Seurat itself not needed)
})

# Forked parallelism only (mclapply); results are seeded per block so they do not depend on N_CORES.
N_CORES <- max(1L, min(as.integer(Sys.getenv("UNIFORM_TABLE_CORES", "20")), parallel::detectCores()))

out_dir <- file.path("results", "revision", "uniform_table")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

B_BOOT <- as.integer(Sys.getenv("UNIFORM_TABLE_B", "1000"))
SEED_BOOT <- 20260903L
EPS <- 0.05
THR <- 0.2
K_MIN <- 2L

SETTINGS <- c("TCGA", "GTEx", "PBMC", "CITEseq", "BRCA_PAM50", "BRCA_IHC")
SETTING_LABELS <- c(
  TCGA = "TCGA miRNA-mRNA", GTEx = "GTEx mRNA-mRNA", PBMC = "PBMC scRNA gene-gene",
  CITEseq = "CITE-seq protein-RNA", BRCA_PAM50 = "BRCA PAM50 miRNA-mRNA", BRCA_IHC = "BRCA IHC miRNA-mRNA"
)
CONTEXT_LABELS <- c(
  TCGA = "cancer cohort", GTEx = "tissue (SMTSD)", PBMC = "coarse cell type",
  CITEseq = "joint WNN cluster", BRCA_PAM50 = "PAM50 subtype", BRCA_IHC = "IHC receptor group"
)

# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------

# Pearson correlations between every row of X (G x n) and every row of Y (M x n):
# G x M matrix. Zero-variance rows give NaN (treated as missing downstream).
cor_rows_product <- function(X, Y) {
  n <- ncol(X)
  sy <- suppressWarnings(scale(t(Y))) / (n - 1)
  sx <- suppressWarnings(t(scale(t(X))))
  sx %*% sy
}

# Subtract per-group row means (group = context label per column).
residualize_by_group <- function(X, g) {
  g <- as.factor(g)
  for (lv in levels(g)) {
    idx <- which(g == lv)
    if (length(idx) == 0) next
    X[, idx] <- X[, idx, drop = FALSE] - rowMeans(X[, idx, drop = FALSE])
  }
  X
}

# Pearson correlation for specific row pairs (ix from X, iy from Y) across columns.
rowcor_pairs <- function(X, Y, ix, iy, chunk = 20000L) {
  zx <- suppressWarnings(t(scale(t(X))))
  zy <- if (identical(X, Y)) zx else suppressWarnings(t(scale(t(Y))))
  denom <- ncol(X) - 1
  out <- numeric(length(ix))
  for (a in seq(1L, length(ix), by = chunk)) {
    b <- min(length(ix), a + chunk - 1L)
    ii <- a:b
    out[ii] <- rowSums(zx[ix[ii], , drop = FALSE] * zy[iy[ii], , drop = FALSE]) / denom
  }
  out[!is.finite(out)] <- NA_real_
  out
}

# Pair-level percentile bootstrap CI (in percent) for a rate = mean(indicator).
# Replicates are split into BOOT_CHUNKS blocks, each seeded deterministically from `seed`,
# so the result depends only on (B, seed, BOOT_CHUNKS) and not on the number of cores used.
BOOT_CHUNKS <- 20L
boot_pct_ci <- function(ind, B, seed) {
  ind <- ind[!is.na(ind)]
  n <- length(ind)
  if (n == 0 || B <= 0) return(c(lo = NA_real_, hi = NA_real_, n = n))
  x <- as.numeric(ind)
  blocks <- split(seq_len(B), cut(seq_len(B), BOOT_CHUNKS, labels = FALSE))
  est <- unlist(parallel::mclapply(seq_along(blocks), function(bi) {
    set.seed(seed + 1000L * bi)
    vapply(blocks[[bi]], function(i) mean(x[sample.int(n, n, replace = TRUE)]), numeric(1))
  }, mc.cores = min(BOOT_CHUNKS, N_CORES)), use.names = FALSE)
  q <- stats::quantile(est, c(0.025, 0.975), names = FALSE)
  c(lo = 100 * q[1], hi = 100 * q[2], n = n)
}

# Two-sided p-value of a Pearson r on nu residual degrees of freedom (t-test on r);
# identical to p_from_r() in scripts/revision/multidonor.R.
p_from_r <- function(r, nu) {
  p <- rep(NA_real_, length(r))
  ok <- is.finite(r) & abs(r) < 1 & nu >= 1
  tt <- r[ok] * sqrt(nu / (1 - r[ok]^2))
  p[ok] <- 2 * stats::pt(-abs(tt), df = nu)
  p[is.finite(r) & abs(r) >= 1 & nu >= 1] <- 0
  p
}

row_max_abs <- function(A) {
  m <- rep(0, nrow(A))
  for (j in seq_len(ncol(A))) m <- pmax(m, A[, j])
  m
}

# -----------------------------------------------------------------------------
# The uniform metric engine
#   r_global, r_resid : numeric vectors (length P)
#   R                 : P x k matrix of within-context r_c (NA where not estimable)
#   n_ctx             : named integer vector, context sizes for the k columns of R
# -----------------------------------------------------------------------------
compute_metrics <- function(r_global, r_resid, R, n_ctx, eps = EPS, thr = THR, k_min = K_MIN,
                            boot_B = 0L, boot_seed = 1L) {
  stopifnot(length(r_global) == nrow(R), length(r_resid) == nrow(R), ncol(R) == length(n_ctx))
  P0 <- nrow(R)
  valid <- is.finite(R) & abs(R) < 1
  k <- rowSums(valid)
  analyzed <- is.finite(r_global) & is.finite(r_resid) & k >= k_min
  drop <- c(
    n_pairs_input = P0,
    n_dropped_nonfinite_r_global = sum(!is.finite(r_global)),
    n_dropped_nonfinite_r_resid = sum(is.finite(r_global) & !is.finite(r_resid)),
    n_dropped_k_lt_kmin = sum(is.finite(r_global) & is.finite(r_resid) & k < k_min)
  )
  R <- R[analyzed, , drop = FALSE]
  valid <- valid[analyzed, , drop = FALSE]
  k <- k[analyzed]
  rg <- r_global[analyzed]
  rr <- r_resid[analyzed]
  P <- length(rg)
  Rv <- R
  Rv[!valid] <- 0

  # (1) pooled r_global
  abs_rg <- abs(rg)
  strong <- abs_rg >= thr

  # (2) within-context distribution
  s1 <- rowSums(Rv)
  s2 <- rowSums(Rv^2)
  sd_rc <- sqrt(pmax(0, (s2 - s1^2 / k) / (k - 1)))
  pos <- rowSums(valid & (Rv >= eps))
  neg <- rowSums(valid & (Rv <= -eps))
  mixed <- pos > 0 & neg > 0
  max_abs <- row_max_abs(abs(Rv))
  strong_ctx <- max_abs >= thr

  # (3) heterogeneity: Fisher z, weights N_c - 3, Cochran Q, I2, BH FDR
  w <- as.numeric(n_ctx) - 3
  W <- matrix(rep(w, each = P), nrow = P) * valid
  z <- contextcorr_fisher_z(Rv)
  z[!valid] <- 0
  sum_w <- rowSums(W)
  sum_wz <- rowSums(W * z)
  sum_wz2 <- rowSums(W * z^2)
  Q <- pmax(0, sum_wz2 - sum_wz^2 / sum_w)
  df <- k - 1
  p_het <- stats::pchisq(Q, df = df, lower.tail = FALSE)
  I2 <- pmax(0, (Q - df) / pmax(Q, 1e-12))
  fdr <- stats::p.adjust(p_het, method = "BH")
  z_fixed <- sum_wz / sum_w
  r_fixed <- contextcorr_inv_fisher_z(z_fixed)

  # (4) random-effects (DerSimonian-Laird), vectorized version of random_effects_pool()
  sum_w2 <- rowSums(W^2)
  C <- sum_w - sum_w2 / sum_w
  tau2 <- pmax(0, (Q - df) / C)
  WRE <- 1 / (1 / W + tau2)   # W == 0 (invalid context) -> weight 0
  zRE <- rowSums(WRE * z) / rowSums(WRE)
  rRE <- contextcorr_inv_fisher_z(zRE)
  rm(W, WRE, z)
  re_disagree <- sign(rRE) != sign(rg) & rg != 0 & rRE != 0

  # (5) between-vs-within diagnostics
  flip <- sign(rr) != sign(rg) & rg != 0 & rr != 0
  # Nominally significant flip = residualized flip in which r_global and r_resid
  # are each individually significant at two-sided p < 0.05 (t-test on r). The analysis population
  # is the union of the retained contexts, so N = sum(n_ctx); r_resid removes k context means,
  # giving nu = N - k - 1 (as in multidonor.R); r_global has nu = N - 2.
  # NOTE: N_pop and k_ctx describe the uniform analysis population (all retained contexts). For the
  # legacy paper-convention reproduction calls (res_paper, res_pbmc_paper, res_cite_paper), r_global is
  # pooled over a larger set, so their nominal-flip metrics are not meaningful and are not reported.
  N_pop <- sum(as.numeric(n_ctx))
  k_ctx <- length(n_ctx)
  p_glob <- p_from_r(rg, N_pop - 2)
  p_res <- p_from_r(rr, N_pop - k_ctx - 1)
  both_sig <- p_glob < 0.05 & p_res < 0.05
  nom_flip <- flip & both_sig
  maj <- sign(pos - neg)
  elig_simp_all <- maj != 0 & sign(rg) != 0
  simpson <- sign(rg) != maj
  elig_simp_strong <- elig_simp_all & strong

  metrics <- c(
    n_pairs_analyzed = P,
    n_contexts_retained = length(n_ctx),
    m1_median_abs_r_global = stats::median(abs_rg),
    m1_pct_abs_r_global_ge_0.2 = 100 * mean(strong),
    n_pairs_abs_r_global_ge_0.2 = sum(strong),
    m2_median_sd_rc = stats::median(sd_rc),
    m2_pct_mixed_sign_eps0.05 = 100 * mean(mixed),
    m2_pct_mixed_sign_eps0.05_maxabs_rc_ge_0.2 = 100 * mean(mixed[strong_ctx]),
    n_pairs_maxabs_rc_ge_0.2 = sum(strong_ctx),
    m2_median_k_valid_contexts = stats::median(k),
    m3_median_I2 = stats::median(I2),
    m3_I2_q25 = stats::quantile(I2, 0.25, names = FALSE),
    m3_I2_q75 = stats::quantile(I2, 0.75, names = FALSE),
    m3_pct_Q_fdr_lt_0.05 = 100 * mean(fdr < 0.05),
    m4_median_abs_r_RE = stats::median(abs(rRE)),
    m4_pct_sign_rRE_ne_rglobal_ge_0.2 = 100 * mean(re_disagree[strong]),
    m4_pct_sign_rRE_ne_rglobal_all = 100 * mean(re_disagree),
    m4_median_abs_r_fixed = stats::median(abs(r_fixed)),
    m5_pct_resid_flip_all = 100 * mean(flip),
    m5_pct_resid_flip_ge_0.2 = 100 * mean(flip[strong]),
    m5_pct_simpson_reversal_ge_0.2 = 100 * mean(simpson[elig_simp_strong]),
    n_pairs_eligible_simpson_ge_0.2 = sum(elig_simp_strong),
    m5_pct_simpson_reversal_all = 100 * mean(simpson[elig_simp_all]),
    n_pairs_eligible_simpson_all = sum(elig_simp_all),
    m5_pct_nominal_flip_all = 100 * mean(nom_flip),
    m5_pct_nominal_flip_ge_0.2 = 100 * mean(nom_flip[strong]),
    n_pairs_both_significant = sum(both_sig),
    m5_pct_nominal_flip_among_both_significant = 100 * mean(nom_flip[both_sig]),
    n_pop_nominal_test = N_pop,
    crit_abs_r_global_p05 = { nu <- N_pop - 2; tq <- stats::qt(0.975, nu); tq / sqrt(nu + tq^2) },
    crit_abs_r_resid_p05 = { nu <- N_pop - k_ctx - 1; tq <- stats::qt(0.975, nu); tq / sqrt(nu + tq^2) },
    n_nominal_flips = sum(nom_flip),
    n_resid_flips = sum(flip),
    median_abs_r_global_nominal_flip = if (any(nom_flip)) stats::median(abs(rg[nom_flip])) else NA_real_,
    median_abs_r_resid_nominal_flip = if (any(nom_flip)) stats::median(abs(rr[nom_flip])) else NA_real_,
    pct_nominal_flips_abs_r_resid_ge_0.1 = if (any(nom_flip)) 100 * mean(abs(rr[nom_flip]) >= 0.1) else NA_real_
  )

  ci <- NULL
  if (boot_B > 0) {
    ci_rev <- boot_pct_ci(simpson[elig_simp_strong], boot_B, boot_seed + 1L)
    ci_flip_all <- boot_pct_ci(flip, boot_B, boot_seed + 2L)
    ci_flip_strong <- boot_pct_ci(flip[strong], boot_B, boot_seed + 3L)
    ci_nom_all <- boot_pct_ci(nom_flip, boot_B, boot_seed + 4L)
    ci_nom_strong <- boot_pct_ci(nom_flip[strong], boot_B, boot_seed + 5L)
    ci <- data.frame(
      rate = c("simpson_reversal_ge_0.2", "resid_flip_all", "resid_flip_ge_0.2", "nominal_flip_all", "nominal_flip_ge_0.2"),
      estimate_pct = c(metrics[["m5_pct_simpson_reversal_ge_0.2"]],
                       metrics[["m5_pct_resid_flip_all"]], metrics[["m5_pct_resid_flip_ge_0.2"]],
                       metrics[["m5_pct_nominal_flip_all"]], metrics[["m5_pct_nominal_flip_ge_0.2"]]),
      ci95_lo_pct = c(ci_rev[["lo"]], ci_flip_all[["lo"]], ci_flip_strong[["lo"]], ci_nom_all[["lo"]], ci_nom_strong[["lo"]]),
      ci95_hi_pct = c(ci_rev[["hi"]], ci_flip_all[["hi"]], ci_flip_strong[["hi"]], ci_nom_all[["hi"]], ci_nom_strong[["hi"]]),
      n_pairs = c(ci_rev[["n"]], ci_flip_all[["n"]], ci_flip_strong[["n"]], ci_nom_all[["n"]], ci_nom_strong[["n"]]),
      B = boot_B,
      seed = boot_seed + 1:5,
      method = "pair-level percentile bootstrap (resample eligible pairs with replacement)",
      stringsAsFactors = FALSE
    )
  }

  list(metrics = metrics, drop = drop, ci = ci,
       pairs = list(I2 = I2, k = k, mixed = mixed, flip = flip, nom_flip = nom_flip, simpson = simpson,
                    elig_simpson_strong = elig_simp_strong, strong = strong, rRE = rRE))
}

results <- list()
cis <- list()
ctx_sizes <- list()
repro <- list()
table2 <- list()
checks <- character()
note_check <- function(...) {
  s <- sprintf(...)
  message("CHECK: ", s)
  checks <<- c(checks, s)
}
add_repro <- function(setting, quantity, published, source, paper_conv, uniform, tol, note = "") {
  repro[[length(repro) + 1]] <<- data.frame(
    setting = setting, quantity = quantity, published_value = published, published_source = source,
    recomputed_paper_convention = paper_conv, recomputed_uniform = uniform,
    paper_convention_reproduces = if (is.na(paper_conv)) NA else abs(paper_conv - published) <= tol,
    uniform_reproduces = if (is.na(uniform)) NA else abs(uniform - published) <= tol,
    rounding_tolerance = tol, note = note, stringsAsFactors = FALSE
  )
}
add_ctx_sizes <- function(setting, tab, retained) {
  ctx_sizes[[setting]] <<- data.frame(setting = setting, context = names(tab), n = as.integer(tab),
                                      retained = names(tab) %in% retained, stringsAsFactors = FALSE)
}

# =============================================================================
# 1. TCGA (cancer cohort)  --  full recomputation on the paper's 10^6 domain
# =============================================================================
t0 <- tic()
message("Loading TCGA matched data ...")
f_matched <- file.path(paths$mirtcga, "matchedNoGBM.rda")
if (!file.exists(f_matched)) f_matched <- file.path(paths$mirtcga, "matchedNoGBM.Rda")
load(f_matched)  # mRNAdata, miRdata, cancerType, sampleType
n_tcga <- ncol(mRNAdata)
tcga_sampletype <- table(sampleType)
G <- 5000; M <- 200
gene_var <- apply(mRNAdata, 1, var)
mir_var <- apply(miRdata, 1, var)
genes_tcga <- names(sort(gene_var, decreasing = TRUE))[1:G]
mirs_tcga <- names(sort(mir_var, decreasing = TRUE))[1:M]
X <- mRNAdata[genes_tcga, , drop = FALSE]
Y <- miRdata[mirs_tcga, , drop = FALSE]
rm(gene_var, mir_var)

cohort_tab <- table(cancerType)
n_min_tcga <- 20L
retained_tcga <- names(cohort_tab)[cohort_tab >= n_min_tcga]
add_ctx_sizes("TCGA", cohort_tab, retained_tcga)
pop <- which(as.character(cancerType) %in% retained_tcga)
ctx <- droplevels(factor(cancerType[pop], levels = retained_tcga))

r_global_tcga <- as.vector(cor_rows_product(X[, pop], Y[, pop]))
R_tcga <- matrix(NA_real_, nrow = G * M, ncol = length(retained_tcga), dimnames = list(NULL, retained_tcga))
for (ct in retained_tcga) {
  idx <- which(ctx == ct)
  R_tcga[, ct] <- as.vector(cor_rows_product(X[, pop[idx], drop = FALSE], Y[, pop[idx], drop = FALSE]))
}
r_resid_tcga <- as.vector(cor_rows_product(residualize_by_group(X[, pop], ctx), residualize_by_group(Y[, pop], ctx)))
n_ctx_tcga <- as.integer(table(ctx)[retained_tcga]); names(n_ctx_tcga) <- retained_tcga

res_tcga <- compute_metrics(r_global_tcga, r_resid_tcga, R_tcga, n_ctx_tcga, boot_B = B_BOOT, boot_seed = SEED_BOOT + 10L)
results[["TCGA"]] <- res_tcga
cis[["TCGA"]] <- cbind(setting = "TCGA", res_tcga$ci)

# Cross-check against the paper's cached pipeline outputs (same domain, same definitions).
het_cache <- readRDS(file.path("results", "cache", "heterogeneity_I2_FDR.rds"))
i2_diff <- max(abs(het_cache$I2[het_cache$ok] - res_tcga$pairs$I2), na.rm = TRUE)
note_check("TCGA: max |I2(recomputed) - I2(results/cache/heterogeneity_I2_FDR.rds)| = %.2e over %d pairs", i2_diff, length(res_tcga$pairs$I2))
rob <- read.csv(file.path("results", "robustness_signflip_summary.csv"))
rob_row <- rob[rob$eps == 0.05 & rob$threshold == 0.2, ]
het_row <- read.csv(file.path("results", "heterogeneity_summary.csv"))[1, ]
add_repro("TCGA", "Simpson reversal % at |r_global|>=0.2 (eps 0.05)", 13.3, "main.tex line 106; results/robustness_signflip_summary.csv frac_simpson (eps 0.05, t 0.2)",
          res_tcga$metrics[["m5_pct_simpson_reversal_ge_0.2"]], res_tcga$metrics[["m5_pct_simpson_reversal_ge_0.2"]], 0.05,
          sprintf("cached pipeline value %.4f%% on %d eligible pairs; recomputed on %d eligible pairs", 100 * rob_row$frac_simpson, rob_row$eligible_simpson, res_tcga$metrics[["n_pairs_eligible_simpson_ge_0.2"]]))
add_repro("TCGA", "median I2", 0.86, "main.tex line 108; results/heterogeneity_summary.csv median_I2",
          res_tcga$metrics[["m3_median_I2"]], res_tcga$metrics[["m3_median_I2"]], 0.005, sprintf("cached %.4f", het_row$median_I2))
add_repro("TCGA", "% pairs with Q significant at FDR<0.05", 99.5, "main.tex line 108; results/heterogeneity_summary.csv frac_fdr_lt_0p05",
          res_tcga$metrics[["m3_pct_Q_fdr_lt_0.05"]], res_tcga$metrics[["m3_pct_Q_fdr_lt_0.05"]], 0.05, sprintf("cached %.4f%%", 100 * het_row$frac_fdr_lt_0p05))
add_repro("TCGA", "% mixed-sign (eps 0.05, all pairs)", 99.5, "Fig 3A; results/robustness_signflip_summary.csv frac_mixed (eps 0.05, t 0)",
          res_tcga$metrics[["m2_pct_mixed_sign_eps0.05"]], res_tcga$metrics[["m2_pct_mixed_sign_eps0.05"]], 0.05,
          sprintf("cached %.4f%%", 100 * rob[rob$eps == 0.05 & rob$threshold == 0, "frac_mixed"]))

table2[["TCGA"]] <- list(
  data_source = "TCGA via FireBrowse standardized data (release 2016-01-28): matched mRNA-seq (RSEM, upper-quartile normalized) and miRNA-seq (RPM); GBM and LAML excluded",
  sample_unit = "bulk sample",
  n_total_in_dataset = n_tcga,
  n_samples_analyzed = length(pop),
  n_features_x = G, n_features_y = M,
  feature_filter_rule = sprintf("top %s mRNAs (of %s) and top %s miRNAs (of %s) by variance across all %s samples", format(G, big.mark = ","), format(nrow(mRNAdata), big.mark = ","), M, format(nrow(miRdata), big.mark = ","), format(n_tcga, big.mark = ",")),
  transform = "mRNA: log2(10 + RSEM normalized count); miRNA: log2(1 + RPM)",
  n_contexts_total = length(cohort_tab), n_contexts_retained = length(retained_tcga),
  min_context_size_rule = sprintf(">= %d samples per cohort", n_min_tcga),
  n_pairs_analyzed = res_tcga$metrics[["n_pairs_analyzed"]],
  pair_selection = sprintf("all %s x %s = %s pairs of the high-variance domain", format(G, big.mark = ","), M, format(G * M, big.mark = ",")),
  notes = sprintf("sampleType in matchedNoGBM.rda: %s tumor, %s normal, %s metastasis (all used, as in the paper)", tcga_sampletype[["tumor"]], tcga_sampletype[["normal"]], tcga_sampletype[["metastasis"]])
)
rm(R_tcga, res_tcga, het_cache); gc(verbose = FALSE)
toc("TCGA", t0)

# =============================================================================
# 5/6. BRCA subtypes (PAM50, IHC) -- reuse caches; recompute r_global/r_resid on the analysis population
# =============================================================================
t0 <- tic()
message("BRCA subtypes ...")
clin <- read.delim(file.path(paths$clean, "BRCA_clinicalMatrix.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
brca_idx <- which(cancerType == "BRCA" & sampleType == "tumor")
brca_ids15 <- gsub("\\.", "-", substr(colnames(mRNAdata)[brca_idx], 1, 15))
m <- match(brca_ids15, clin$sampleID)
stopifnot(all(!is.na(m)))
pam50_raw <- clin$PAM50Call_RNAseq[m]
pam50 <- ifelse(pam50_raw %in% c("Basal", "Her2", "LumA", "LumB", "Normal"), pam50_raw, NA_character_)
pam50_tab <- table(pam50, useNA = "no")
pam50 <- ifelse(pam50 %in% names(pam50_tab)[pam50_tab >= 10], pam50, NA_character_)
er <- clin$ER_Status_nature2012[m]; pr <- clin$PR_Status_nature2012[m]; h2 <- clin$HER2_Final_Status_nature2012[m]
er <- ifelse(er %in% c("Positive", "Negative"), er, NA); pr <- ifelse(pr %in% c("Positive", "Negative"), pr, NA); h2 <- ifelse(h2 %in% c("Positive", "Negative"), h2, NA)
ihc <- rep(NA_character_, length(er)); ok_ihc <- !is.na(er) & !is.na(pr) & !is.na(h2)
ihc[ok_ihc & h2 == "Positive"] <- "HER2+"
ihc[ok_ihc & h2 == "Negative" & (er == "Positive" | pr == "Positive")] <- "HR+"
ihc[ok_ihc & h2 == "Negative" & er == "Negative" & pr == "Negative"] <- "TNBC"

mRNA_brca <- mRNAdata[, brca_idx, drop = FALSE]
miR_brca <- miRdata[, brca_idx, drop = FALSE]
gv <- apply(mRNA_brca, 1, var); mv <- apply(miR_brca, 1, var)
genes_brca <- names(sort(gv, decreasing = TRUE))[1:G]
mirs_brca <- names(sort(mv, decreasing = TRUE))[1:M]
Xb <- mRNA_brca[genes_brca, , drop = FALSE]
Yb <- miR_brca[mirs_brca, , drop = FALSE]
rm(mRNA_brca, miR_brca, gv, mv, mRNAdata, miRdata); gc(verbose = FALSE)

brca_rates <- read.csv(file.path("results", "brca_subtype_rates.csv"))
brca_het <- read.csv(file.path("results", "brca_subtype_heterogeneity.csv"))

run_brca <- function(label_vec, label_name, setting, seed_off, published) {
  cache <- readRDS(file.path("results", "cache", sprintf("brca_subtype_%s.rds", label_name)))
  stopifnot(identical(rownames(cache$r_fixed), genes_brca), identical(colnames(cache$r_fixed), mirs_brca))
  labeled <- which(!is.na(label_vec))
  lab_tab <- table(label_vec[labeled])
  n_min <- 20L
  retained <- names(lab_tab)[lab_tab >= n_min]
  add_ctx_sizes(setting, lab_tab, retained)
  # Paper convention: r_global pooled over ALL labeled tumors (cache); check label reconstruction.
  rg_paper_re <- as.vector(cor_rows_product(Xb[, labeled], Yb[, labeled]))
  note_check("%s: max |r_global(recomputed, n=%d) - cache$r_global| = %.2e", setting, length(labeled), max(abs(rg_paper_re - as.vector(cache$r_global)), na.rm = TRUE))
  R <- sapply(retained, function(l) as.vector(cache$cor_array[, , l]))
  n_ctx <- as.integer(lab_tab[retained]); names(n_ctx) <- retained
  ctx_all <- factor(label_vec[labeled])
  rr_paper <- as.vector(cor_rows_product(residualize_by_group(Xb[, labeled], ctx_all), residualize_by_group(Yb[, labeled], ctx_all)))
  res_paper <- compute_metrics(as.vector(cache$r_global), rr_paper, R, n_ctx, boot_B = 0L)
  # Uniform convention: analysis population = tumors in retained contexts.
  pop <- labeled[label_vec[labeled] %in% retained]
  ctx <- droplevels(factor(label_vec[pop], levels = retained))
  rg <- as.vector(cor_rows_product(Xb[, pop], Yb[, pop]))
  rr <- as.vector(cor_rows_product(residualize_by_group(Xb[, pop], ctx), residualize_by_group(Yb[, pop], ctx)))
  res <- compute_metrics(rg, rr, R, n_ctx, boot_B = B_BOOT, boot_seed = SEED_BOOT + seed_off)
  results[[setting]] <<- res
  cis[[setting]] <<- cbind(setting = setting, res$ci)
  rr_row <- brca_rates[brca_rates$context == label_name & brca_rates$eps == 0.05 & brca_rates$threshold == 0.2, ]
  hh <- brca_het[brca_het$context == label_name, ]
  add_repro(setting, "Simpson reversal % at |r_global|>=0.2 (eps 0.05)", published$rev, "main.tex line 128; results/brca_subtype_rates.csv frac_simpson (eps 0.05, t 0.2)",
            res_paper$metrics[["m5_pct_simpson_reversal_ge_0.2"]], res$metrics[["m5_pct_simpson_reversal_ge_0.2"]], 0.05,
            sprintf("cached %.4f%% (%d eligible); paper convention pools r_global over all %d labeled tumors, uniform over %d tumors in contexts with >= %d", 100 * rr_row$frac_simpson, rr_row$eligible_simpson, length(labeled), length(pop), n_min))
  add_repro(setting, "% mixed-sign (eps 0.05) among max_c|r_c|>=0.2", published$mixed, "main.tex line 128; results/brca_subtype_rates.csv frac_mixed (eps 0.05, t 0.2)",
            res_paper$metrics[["m2_pct_mixed_sign_eps0.05_maxabs_rc_ge_0.2"]], res$metrics[["m2_pct_mixed_sign_eps0.05_maxabs_rc_ge_0.2"]], 0.05, sprintf("cached %.4f%%", 100 * rr_row$frac_mixed))
  add_repro(setting, "median I2", published$i2, "main.tex line 128; results/brca_subtype_heterogeneity.csv median_I2",
            res_paper$metrics[["m3_median_I2"]], res$metrics[["m3_median_I2"]], 0.005, sprintf("cached %.4f", hh$median_I2))
  table2[[setting]] <<- list(
    data_source = sprintf("TCGA BRCA primary tumors from the matched set (n = %d); %s from UCSC Xena BRCA_clinicalMatrix (%s)", length(brca_idx), if (label_name == "PAM50") "PAM50 call" else "ER/PR/HER2 status", if (label_name == "PAM50") "PAM50Call_RNAseq" else "ER_Status_nature2012, PR_Status_nature2012, HER2_Final_Status_nature2012"),
    sample_unit = "bulk primary tumor",
    n_total_in_dataset = length(labeled),
    n_samples_analyzed = length(pop),
    n_features_x = G, n_features_y = M,
    feature_filter_rule = sprintf("top %s mRNAs and top %s miRNAs by variance within the %d BRCA primary tumors", format(G, big.mark = ","), M, length(brca_idx)),
    transform = "mRNA: log2(10 + RSEM normalized count); miRNA: log2(1 + RPM)",
    n_contexts_total = length(lab_tab), n_contexts_retained = length(retained),
    min_context_size_rule = sprintf(">= %d tumors per subtype", n_min),
    n_pairs_analyzed = res$metrics[["n_pairs_analyzed"]],
    pair_selection = sprintf("all %s x %s = %s pairs", format(G, big.mark = ","), M, format(G * M, big.mark = ",")),
    notes = if (length(retained) < length(lab_tab)) sprintf("excluded from within-context statistics and from the analysis population: %s", paste(sprintf("%s (n = %d)", setdiff(names(lab_tab), retained), lab_tab[setdiff(names(lab_tab), retained)]), collapse = ", ")) else "all labeled tumors analyzed"
  )
  invisible(NULL)
}
run_brca(pam50, "PAM50", "BRCA_PAM50", 50L, list(rev = 5.5, mixed = 50.6, i2 = 0.14))
run_brca(ihc, "IHC", "BRCA_IHC", 60L, list(rev = 0.35, mixed = 28.6, i2 = 0.00))
rm(Xb, Yb, X, Y); gc(verbose = FALSE)
toc("BRCA", t0)

# =============================================================================
# 2. GTEx (SMTSD tissue) -- recomputed (no cache exists for the SMTSD-based numbers)
# =============================================================================
t0 <- tic()
message("GTEx ...")
f_gtex <- file.path(paths$mirtcga, "xena-gtex.Rda")
if (!file.exists(f_gtex)) stop("Missing GTEx input: ", f_gtex)
load(f_gtex)  # gtex$mrna: genes x samples, log2(TPM + 0.001)
Xg <- as.matrix(gtex$mrna); rm(gtex)
attr_file <- Sys.getenv("GTEX_ATTRIBUTES", file.path(paths$mirtcga, "GTEx_SampleAttributes_v8.txt"))
attr <- utils::read.delim(attr_file, stringsAsFactors = FALSE, quote = "")
smtsd <- attr$SMTSD[match(gsub("\\.", "-", colnames(Xg)), attr$SAMPID)]
n_gtex_total <- ncol(Xg); n_gtex_unmatched <- sum(is.na(smtsd))
tis_tab <- sort(table(smtsd), decreasing = TRUE)
min_n_gtex <- 150L
keep_tissues <- names(tis_tab)[tis_tab >= min_n_gtex]
keep_tissues <- keep_tissues[seq_len(min(20, length(keep_tissues)))]
add_ctx_sizes("GTEx", tis_tab, keep_tissues)

# The Rmd's pair sampler (408-gtex-mrna-mrna.Rmd), verbatim.
sample_pairs <- function(G, P) {
  i <- sample.int(G, size = 2 * P, replace = TRUE)
  j <- sample.int(G, size = 2 * P, replace = TRUE)
  ok <- i != j
  i <- i[ok]; j <- j[ok]
  a <- pmin(i, j); b <- pmax(i, j)
  key <- paste(a, b, sep = "_")
  u <- !duplicated(key)
  a <- a[u]; b <- b[u]
  if (length(a) < P) stop("Not enough unique pairs sampled.")
  list(i = a[1:P], j = b[1:P])
}
G_gtex <- 1500L; P_gtex <- 200000L
set.seed(1)
pairs_gtex <- sample_pairs(G_gtex, P_gtex)
ii_g <- pairs_gtex$i; jj_g <- pairs_gtex$j

pop <- which(smtsd %in% keep_tissues)
n_gtex_pop <- length(pop)
n_ctx_gtex <- as.integer(tis_tab[keep_tissues]); names(n_ctx_gtex) <- keep_tissues

# Gene ranking under both candidate selection rules (the published GTEx numbers came
# from a re-run with SMTSD labels whose script was not saved, so we evaluate the grid).
gvar_all <- apply(Xg, 1, var)
gvar_20 <- apply(Xg[, pop, drop = FALSE], 1, var)
genes_by_all <- names(sort(gvar_all, decreasing = TRUE))[1:G_gtex]
genes_by_20 <- names(sort(gvar_20, decreasing = TRUE))[1:G_gtex]
n_genes_gtex <- length(gvar_all)
n_genes_overlap <- length(intersect(genes_by_all, genes_by_20))

# One evaluation of the GTEx pipeline for a given (gene set, population) combination.
gtex_eval <- function(genes, pop_mode) {
  Xs <- Xg[genes, , drop = FALSE]
  Rl <- parallel::mclapply(keep_tissues, function(tt) {
    idx <- which(smtsd == tt)
    rowcor_pairs(Xs[, idx, drop = FALSE], Xs[, idx, drop = FALSE], ii_g, jj_g)
  }, mc.cores = min(length(keep_tissues), N_CORES))
  R <- matrix(unlist(Rl, use.names = FALSE), nrow = P_gtex, dimnames = list(NULL, keep_tissues))
  if (pop_mode == "all_samples") {
    # 408-gtex-mrna-mrna.Rmd path: correlations over every sample; only kept-tissue
    # columns are mean-centered, the rest enter r_resid unresidualized.
    use <- seq_len(ncol(Xs))
    Xr <- Xs
    for (tt in keep_tissues) {
      idx <- which(smtsd == tt)
      Xr[, idx] <- Xs[, idx, drop = FALSE] - rowMeans(Xs[, idx, drop = FALSE])
    }
  } else {
    use <- pop
    Xr <- residualize_by_group(Xs[, pop, drop = FALSE], droplevels(factor(smtsd[pop], levels = keep_tissues)))
  }
  rg <- rowcor_pairs(Xs[, use, drop = FALSE], Xs[, use, drop = FALSE], ii_g, jj_g)
  rr <- rowcor_pairs(Xr, Xr, ii_g, jj_g)
  list(r_global = rg, r_resid = rr, R = R)
}

# Summaries under the two sign rules, for the convention grid.
gtex_summaries <- function(ev) {
  R <- ev$R; rg <- ev$r_global; rr <- ev$r_resid
  valid <- is.finite(R) & abs(R) < 1
  k <- rowSums(valid)
  Rv <- R; Rv[!valid] <- 0
  maj0 <- rep(NA_integer_, P_gtex)
  maj0[k > 0] <- sign(rowSums(valid & Rv > 0) - rowSums(valid & Rv < 0))[k > 0]
  rev0 <- k > 0 & sign(rg) != maj0 & sign(rg) != 0 & maj0 != 0
  maj5 <- sign(rowSums(valid & Rv >= EPS) - rowSums(valid & Rv <= -EPS))
  el5 <- maj5 != 0 & sign(rg) != 0 & k >= K_MIN
  rev5 <- sign(rg) != maj5
  flip <- sign(rg) != sign(rr) & rg != 0 & rr != 0
  W <- matrix(rep(as.numeric(n_ctx_gtex) - 3, each = P_gtex), nrow = P_gtex) * valid
  z <- contextcorr_fisher_z(Rv); z[!valid] <- 0
  Q <- pmax(0, rowSums(W * z^2) - rowSums(W * z)^2 / rowSums(W))
  df <- pmax(k - 1, 0)
  I2 <- rep(NA_real_, P_gtex); okq <- df >= 1 & is.finite(Q)
  I2[okq] <- pmax(0, (Q[okq] - df[okq]) / pmax(Q[okq], 1e-12))
  c(disagree_eps0_all = 100 * mean(rev0, na.rm = TRUE),
    disagree_eps0_ge_0.2 = 100 * mean(rev0[abs(rg) >= THR], na.rm = TRUE),
    reversal_eps0.05_all = 100 * mean(rev5[el5]),
    reversal_eps0.05_ge_0.2 = 100 * mean(rev5[el5 & abs(rg) >= THR]),
    resid_flip_all = 100 * mean(flip, na.rm = TRUE),
    resid_flip_ge_0.2 = 100 * mean(flip[abs(rg) >= THR], na.rm = TRUE),
    median_I2 = stats::median(I2, na.rm = TRUE))
}

grid_defs <- list(
  list(id = "genes_all_pop_all", genes = genes_by_all, pop = "all_samples",
       label = "genes by variance over all samples; correlations over all samples (literal 408-gtex-mrna-mrna.Rmd path)"),
  list(id = "genes_all_pop_20tissues", genes = genes_by_all, pop = "tissues20",
       label = "genes by variance over all samples; correlations over the 20-tissue population"),
  list(id = "genes_20tissues_pop_all", genes = genes_by_20, pop = "all_samples",
       label = "genes by variance within the 20 tissues; correlations over all samples"),
  list(id = "genes_20tissues_pop_20tissues", genes = genes_by_20, pop = "tissues20",
       label = "genes by variance within the 20 tissues; correlations over the 20-tissue population (UNIFORM convention used in Table 1)")
)
gtex_grid <- list()
ev_uniform <- NULL
for (gd in grid_defs) {
  message("  GTEx convention: ", gd$id)
  ev <- gtex_eval(gd$genes, gd$pop)
  gtex_grid[[gd$id]] <- data.frame(convention = gd$id, description = gd$label,
                                   t(gtex_summaries(ev)), check.names = FALSE, stringsAsFactors = FALSE)
  if (gd$id == "genes_20tissues_pop_20tissues") ev_uniform <- ev else rm(ev)
  gc(verbose = FALSE)
}
gtex_grid_df <- do.call(rbind, gtex_grid); rownames(gtex_grid_df) <- NULL
gtex_grid_df$published_reference <- "main.tex line 114: 21.0% disagree (all), 6.2% (>=0.2), 23.5% residualized flip, median I2 0.71"
write.csv(gtex_grid_df, file.path(out_dir, "gtex_convention_grid.csv"), row.names = FALSE)

genes_gtex <- genes_by_20
res_gtex <- compute_metrics(ev_uniform$r_global, ev_uniform$r_resid, ev_uniform$R, n_ctx_gtex,
                            boot_B = B_BOOT, boot_seed = SEED_BOOT + 20L)
results[["GTEx"]] <- res_gtex
cis[["GTEx"]] <- cbind(setting = "GTEx", res_gtex$ci)

gg <- function(conv, col) gtex_grid_df[[col]][gtex_grid_df$convention == conv]
best_for <- function(col, target) {
  d <- abs(gtex_grid_df[[col]] - target)
  gtex_grid_df$convention[which.min(d)]
}
add_repro("GTEx", "% pooled sign disagreeing with tissue majority (all pairs)", 21.0,
          "main.tex line 114 (eps = 0 majority, ties counted in the denominator)",
          gg("genes_20tissues_pop_20tissues", "disagree_eps0_all"), res_gtex$metrics[["m5_pct_simpson_reversal_all"]], 0.05,
          sprintf("literal Rmd path gives %.2f%%; the 20-tissue population with 20-tissue gene ranking gives %.2f%% (closest convention: %s). Uniform value uses eps = 0.05 with ties excluded.",
                  gg("genes_all_pop_all", "disagree_eps0_all"), gg("genes_20tissues_pop_20tissues", "disagree_eps0_all"), best_for("disagree_eps0_all", 21.0)))
add_repro("GTEx", "% disagreement among |r_global|>=0.2", 6.2,
          "main.tex line 114 (eps = 0 majority)",
          gg("genes_20tissues_pop_20tissues", "disagree_eps0_ge_0.2"), res_gtex$metrics[["m5_pct_simpson_reversal_ge_0.2"]], 0.05,
          sprintf("literal Rmd path gives %.2f%%; eps = 0 on the 20-tissue population gives %.2f%%; the uniform eps = 0.05 rule on the same population gives %.2f%%",
                  gg("genes_all_pop_all", "disagree_eps0_ge_0.2"), gg("genes_20tissues_pop_20tissues", "disagree_eps0_ge_0.2"), res_gtex$metrics[["m5_pct_simpson_reversal_ge_0.2"]]))
add_repro("GTEx", "% residualized sign flip (all pairs)", 23.5, "main.tex line 114",
          gg("genes_20tissues_pop_20tissues", "resid_flip_all"), res_gtex$metrics[["m5_pct_resid_flip_all"]], 0.05,
          sprintf("literal Rmd path gives %.2f%%; %.2f%% under genes-by-all-samples with the 20-tissue population (closest convention: %s)",
                  gg("genes_all_pop_all", "resid_flip_all"), gg("genes_all_pop_20tissues", "resid_flip_all"), best_for("resid_flip_all", 23.5)))
add_repro("GTEx", "median I2", 0.71, "main.tex line 114",
          gg("genes_20tissues_pop_20tissues", "median_I2"), res_gtex$metrics[["m3_median_I2"]], 0.005,
          sprintf("literal Rmd path (genes ranked over all samples) gives %.4f; ranking genes within the 20 tissues gives %.4f",
                  gg("genes_all_pop_all", "median_I2"), gg("genes_20tissues_pop_20tissues", "median_I2")))

table2[["GTEx"]] <- list(
  data_source = "GTEx v8 via UCSC Xena TCGA-TARGET-GTEx expression matrix (log2(TPM + 0.001)); tissue = SMTSD field of GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt joined on sample ID",
  sample_unit = "bulk tissue sample",
  n_total_in_dataset = n_gtex_total,
  n_samples_analyzed = n_gtex_pop,
  n_features_x = G_gtex, n_features_y = G_gtex,
  feature_filter_rule = sprintf("%s most variable genes (of %s) among the %s samples of the 20 retained tissues", format(G_gtex, big.mark = ","), format(n_genes_gtex, big.mark = ","), format(n_gtex_pop, big.mark = ",")),
  transform = "log2(TPM + 0.001) as distributed by Xena",
  n_contexts_total = length(tis_tab), n_contexts_retained = length(keep_tissues),
  min_context_size_rule = sprintf(">= %d samples per tissue, 20 largest tissues kept", min_n_gtex),
  n_pairs_analyzed = res_gtex$metrics[["n_pairs_analyzed"]],
  pair_selection = sprintf("%s unordered gene pairs sampled uniformly without replacement from %s possible (seed 1)", format(P_gtex, big.mark = ","), format(choose(G_gtex, 2), big.mark = ",")),
  notes = sprintf("%d of %d Xena samples have no SMTSD match; %d tissues have >= %d samples, the 21st (%s, n = %d) is dropped by the 20-tissue cap; %d of the 1,500 genes are shared between the two candidate variance rankings (see gtex_convention_grid.csv)", n_gtex_unmatched, n_gtex_total, sum(tis_tab >= min_n_gtex), min_n_gtex, names(tis_tab)[21], tis_tab[21], n_genes_overlap)
)
rm(Xg, gvar_all, gvar_20, ev_uniform, res_gtex); gc(verbose = FALSE)
toc("GTEx", t0)

# =============================================================================
# 3/4. Single-cell settings -- reuse cached r_c / labels; recompute r_global on the analysis population
# =============================================================================
get_counts <- function(obj, assay) {
  tryCatch(SeuratObject::LayerData(obj[[assay]], layer = "counts"),
           error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = "counts"))
}
# Seurat LogNormalize: log1p(count / total_count_in_cell * 1e4)
lognorm <- function(cnt, scale_factor = 1e4) {
  tot <- Matrix::colSums(cnt)
  out <- cnt
  out@x <- log1p(out@x / rep.int(tot, diff(out@p)) * scale_factor)
  out
}
# Seurat QC used by sc_preprocess(): >= 200 detected features, <= 20% mitochondrial ("^MT-")
qc_keep <- function(cnt, min_features = 200, max_mt = 20) {
  nfeat <- Matrix::colSums(cnt > 0)
  mt <- grepl("^MT-", rownames(cnt))
  pmt <- 100 * Matrix::colSums(cnt[mt, , drop = FALSE]) / Matrix::colSums(cnt)
  nfeat >= min_features & pmt <= max_mt
}
# Seurat CLR with margin = 2: per cell across features
clr_cell <- function(x) log1p(x / exp(sum(log1p(x[x > 0]), na.rm = TRUE) / length(x)))

t0 <- tic()
message("PBMC scRNA ...")
seu <- readRDS(file.path(paths$scratch, "sc_inputs", "pbmc_scrna.rds"))
cnt <- get_counts(seu, "RNA"); n_pbmc_raw <- ncol(cnt); n_genes_pbmc_raw <- nrow(cnt)
keep <- qc_keep(cnt)
ln <- lognorm(cnt[, keep]); rm(seu, cnt)
ar <- readRDS(file.path("results", "sc", "sc_100_pbmc_coarse_celltypes", "analysis_ready.rds"))
inp <- readRDS(file.path("results", "sc", "sc_100_pbmc_coarse_celltypes", "inputs.rds"))
feats <- inp$feature_x
stopifnot(all(feats %in% rownames(ln)), length(ar$context) == ncol(ln))
E <- as.matrix(ln[feats, ]); rm(ln)
ctx_pbmc <- ar$context
stopifnot(identical(as.integer(table(ctx_pbmc)[names(inp$context_table)]), as.integer(inp$context_table)))
ix <- match(ar$pairs$feature_x, feats); iy <- match(ar$pairs$feature_y, feats)
stopifnot(all(ix == ar$pairs$idx_x), all(iy == ar$pairs$idx_y))
n_min_sc <- 200L
pbmc_tab <- table(ctx_pbmc)
retained_pbmc <- ar$context_levels[as.integer(pbmc_tab[ar$context_levels]) >= n_min_sc]
add_ctx_sizes("PBMC", pbmc_tab, retained_pbmc)
# Paper-convention check: r_global over all QC-passing cells must equal the cache.
rg_all_pbmc <- rowcor_pairs(E, E, ix, iy)
note_check("PBMC: max |r_global(recomputed, all %d cells) - cached| = %.2e", ncol(E), max(abs(rg_all_pbmc - ar$summary_pairs$r_global), na.rm = TRUE))
pop <- which(ctx_pbmc %in% retained_pbmc)
ctx <- droplevels(factor(ctx_pbmc[pop], levels = retained_pbmc))
Ek <- E[, pop]
rg_pbmc <- rowcor_pairs(Ek, Ek, ix, iy)
Er <- residualize_by_group(Ek, ctx)
rr_pbmc <- rowcor_pairs(Er, Er, ix, iy)
note_check("PBMC: max |r_resid(recomputed, %d cells in %d contexts) - cached r_resid| = %.2e", length(pop), length(retained_pbmc), max(abs(rr_pbmc - ar$summary_pairs$r_resid), na.rm = TRUE))
R_pbmc <- ar$r_by_context[, retained_pbmc, drop = FALSE]
n_ctx_pbmc <- as.integer(pbmc_tab[retained_pbmc]); names(n_ctx_pbmc) <- retained_pbmc
res_pbmc <- compute_metrics(rg_pbmc, rr_pbmc, R_pbmc, n_ctx_pbmc, boot_B = B_BOOT, boot_seed = SEED_BOOT + 30L)
res_pbmc_paper <- compute_metrics(ar$summary_pairs$r_global, ar$summary_pairs$r_resid, R_pbmc, n_ctx_pbmc, boot_B = 0L)
results[["PBMC"]] <- res_pbmc
cis[["PBMC"]] <- cbind(setting = "PBMC", res_pbmc$ci)
sp <- ar$summary_pairs; okf <- is.finite(sp$r_global) & is.finite(sp$r_resid)
paper_flip_pbmc <- 100 * mean(sign(sp$r_global[okf]) != sign(sp$r_resid[okf]) & sp$r_global[okf] != 0 & sp$r_resid[okf] != 0)
curve_pbmc <- read.csv(file.path("results", "sc", "sc_100_pbmc_coarse_celltypes", "reversal_rate_curve.csv"))
paper_rev_pbmc <- 100 * curve_pbmc[curve_pbmc$eps == 0.05 & abs(curve_pbmc$abs_r_global_threshold - 0.2) < 1e-9, "reversal_rate"]
add_repro("PBMC", "% residualized sign flip (all pairs)", 13.1, "main.tex line 120; scripts/sc_make_paper_figures_v2.R (cached r_global on all cells, cached r_resid)",
          paper_flip_pbmc, res_pbmc$metrics[["m5_pct_resid_flip_all"]], 0.05,
          sprintf("paper convention: r_global over all %d cells, r_resid over the %d cells in contexts with >= %d cells; uniform: both over the %d cells", ncol(E), length(pop), n_min_sc, length(pop)))
add_repro("PBMC", "Simpson reversal % at |r_global|>=0.2 (eps 0.05)", 10.3, "supplement.tex line 82; results/sc/sc_100_pbmc_coarse_celltypes/reversal_rate_curve.csv (eps 0.05, thr 0.2)",
          res_pbmc_paper$metrics[["m5_pct_simpson_reversal_ge_0.2"]], res_pbmc$metrics[["m5_pct_simpson_reversal_ge_0.2"]], 0.05,
          sprintf("cached curve value %.4f%%; paper convention recomputed here %.4f%% (%d eligible)", paper_rev_pbmc, res_pbmc_paper$metrics[["m5_pct_simpson_reversal_ge_0.2"]], res_pbmc_paper$metrics[["n_pairs_eligible_simpson_ge_0.2"]]))
table2[["PBMC"]] <- list(
  data_source = sprintf("10x Genomics PBMC 3k (single healthy donor; Cell Ranger 1.1.0, hg19): %s cells x %s genes; local copy %s", format(n_pbmc_raw, big.mark = ","), format(n_genes_pbmc_raw, big.mark = ","), "sc_inputs/pbmc_scrna.rds"),
  sample_unit = "cell",
  n_total_in_dataset = ncol(E),
  n_samples_analyzed = length(pop),
  n_features_x = length(feats), n_features_y = length(feats),
  feature_filter_rule = sprintf("%s genes: top 2,000 highly variable genes (Seurat vst) plus %d marker genes, detected in >= 1%% of cells, capped at 1,800", format(length(feats), big.mark = ","), inp$extra$feature_selection$marker_panel_size),
  transform = "QC: >= 200 detected genes and <= 20% mitochondrial UMIs; Seurat LogNormalize (UMI / cell total x 1e4, log1p)",
  n_contexts_total = length(pbmc_tab), n_contexts_retained = length(retained_pbmc),
  min_context_size_rule = sprintf(">= %d cells per cell type", n_min_sc),
  n_pairs_analyzed = res_pbmc$metrics[["n_pairs_analyzed"]],
  pair_selection = sprintf("%s gene pairs sampled without replacement from %s unordered pairs (seed 100)", format(nrow(ar$pairs), big.mark = ","), format(choose(length(feats), 2), big.mark = ",")),
  notes = sprintf("coarse cell types by marker-gene scoring assigned before correlation analysis; labels with < %d cells excluded from the analysis population: %s", n_min_sc, paste(sprintf("%s (n = %d)", setdiff(names(pbmc_tab), retained_pbmc), pbmc_tab[setdiff(names(pbmc_tab), retained_pbmc)]), collapse = ", "))
)
rm(E, Ek, Er, ar, R_pbmc, res_pbmc); gc(verbose = FALSE)
toc("PBMC", t0)

t0 <- tic()
message("CITE-seq ...")
seu <- readRDS(file.path(paths$scratch, "sc_inputs", "pbmc_citeseq.rds"))
rna <- get_counts(seu, "RNA"); adt <- get_counts(seu, "ADT")
n_cite_raw <- ncol(rna); n_genes_cite_raw <- nrow(rna); n_adt_raw <- nrow(adt)
keep <- qc_keep(rna)
lnr <- lognorm(rna[, keep])
adt_clr <- apply(as.matrix(adt[, keep]), 2, clr_cell); dimnames(adt_clr) <- dimnames(adt[, keep])
rm(seu, rna, adt)
ac <- readRDS(file.path("results", "sc", "sc_200_citeseq_protein_rna_contexts", "analysis_ready.rds"))
inpc <- readRDS(file.path("results", "sc", "sc_200_citeseq_protein_rna_contexts", "inputs.rds"))
rf <- inpc$feature_sets$rna_features; af <- inpc$feature_sets$adt_features
stopifnot(all(rf %in% rownames(lnr)), all(af %in% rownames(adt_clr)))
ER <- as.matrix(lnr[rf, ]); EA <- adt_clr[af, , drop = FALSE]; rm(lnr, adt_clr)
ctx_cite <- ac$context_defs$Joint_RNA_ADT
stopifnot(length(ctx_cite) == ncol(ER))
ix <- match(ac$pairs$feature_x, af); iy <- match(ac$pairs$feature_y, rf)
stopifnot(all(ix == ac$pairs$idx_x), all(iy == ac$pairs$idx_y))
st <- ac$per_definition_stats$Joint_RNA_ADT
spj <- ac$summary_pairs[ac$summary_pairs$context_definition == "Joint_RNA_ADT", ]
stopifnot(nrow(spj) == nrow(ac$pairs))
cite_tab <- table(ctx_cite)
retained_cite <- st$context_levels[as.integer(cite_tab[st$context_levels]) >= n_min_sc]
add_ctx_sizes("CITEseq", cite_tab, retained_cite)
rg_all_cite <- rowcor_pairs(EA, ER, ix, iy)
note_check("CITE-seq: max |r_global(recomputed, all %d cells) - cached| = %.2e", ncol(ER), max(abs(rg_all_cite - spj$r_global), na.rm = TRUE))
pop <- which(ctx_cite %in% retained_cite)
ctx <- droplevels(factor(ctx_cite[pop], levels = retained_cite))
rg_cite <- rowcor_pairs(EA[, pop, drop = FALSE], ER[, pop], ix, iy)
rr_cite <- rowcor_pairs(residualize_by_group(EA[, pop, drop = FALSE], ctx), residualize_by_group(ER[, pop], ctx), ix, iy)
note_check("CITE-seq: max |r_resid(recomputed, %d cells in %d clusters) - cached r_resid| = %.2e", length(pop), length(retained_cite), max(abs(rr_cite - spj$r_resid), na.rm = TRUE))
R_cite <- st$r_by_context[, retained_cite, drop = FALSE]
n_ctx_cite <- as.integer(cite_tab[retained_cite]); names(n_ctx_cite) <- retained_cite
res_cite <- compute_metrics(rg_cite, rr_cite, R_cite, n_ctx_cite, boot_B = B_BOOT, boot_seed = SEED_BOOT + 40L)
res_cite_paper <- compute_metrics(spj$r_global, spj$r_resid, R_cite, n_ctx_cite, boot_B = 0L)
results[["CITEseq"]] <- res_cite
cis[["CITEseq"]] <- cbind(setting = "CITEseq", res_cite$ci)
okf <- is.finite(spj$r_global) & is.finite(spj$r_resid)
paper_flip_cite <- 100 * mean(sign(spj$r_global[okf]) != sign(spj$r_resid[okf]))
add_repro("CITEseq", "% residualized sign flip (all pairs)", 37.9, "main.tex line 122; scripts/gb_fig2C_citeseq_globalvsresid.R (cached r_global on all cells, cached r_resid)",
          paper_flip_cite, res_cite$metrics[["m5_pct_resid_flip_all"]], 0.05,
          sprintf("paper convention: r_global over all %d cells, r_resid over the %d cells in the %d clusters with >= %d cells; uniform: both over the %d cells", ncol(ER), length(pop), length(retained_cite), n_min_sc, length(pop)))
table2[["CITEseq"]] <- list(
  data_source = sprintf("10x Genomics PBMC 10k with TotalSeq-B antibodies (pbmc_10k_protein_v3, Cell Ranger v3 filtered_feature_bc_matrix): %s cells x %s genes + %d antibody-derived tags; local copy %s", format(n_cite_raw, big.mark = ","), format(n_genes_cite_raw, big.mark = ","), n_adt_raw, "sc_inputs/pbmc_citeseq.rds"),
  sample_unit = "cell",
  n_total_in_dataset = ncol(ER),
  n_samples_analyzed = length(pop),
  n_features_x = length(af), n_features_y = length(rf),
  feature_filter_rule = sprintf("%d ADTs (%d minus %d isotype controls) x %s genes (top 1,500 highly variable genes plus marker panel, detected in >= 1%% of cells, capped at 1,300)", length(af), n_adt_raw, n_adt_raw - length(af), format(length(rf), big.mark = ",")),
  transform = "QC on RNA: >= 200 detected genes and <= 20% mitochondrial UMIs; RNA Seurat LogNormalize; ADT centered log-ratio per cell (Seurat CLR, margin 2)",
  n_contexts_total = length(cite_tab), n_contexts_retained = length(retained_cite),
  min_context_size_rule = sprintf(">= %d cells per cluster", n_min_sc),
  n_pairs_analyzed = res_cite$metrics[["n_pairs_analyzed"]],
  pair_selection = sprintf("all %d x %s = %s protein-RNA pairs", length(af), format(length(rf), big.mark = ","), format(length(af) * length(rf), big.mark = ",")),
  notes = sprintf("contexts = Seurat weighted-nearest-neighbor (WNN) joint RNA+ADT Louvain clusters (resolution 0.6): %d clusters; the %d clusters with < %d cells (%s cells) are excluded from the analysis population", length(cite_tab), length(cite_tab) - length(retained_cite), n_min_sc, format(ncol(ER) - length(pop), big.mark = ","))
)
rm(ER, EA, ac, R_cite, res_cite); gc(verbose = FALSE)
toc("CITEseq", t0)

# =============================================================================
# Assemble outputs
# =============================================================================
t0 <- tic()
metric_names <- names(results[["TCGA"]]$metrics)
tab1 <- data.frame(setting = SETTINGS, setting_label = SETTING_LABELS[SETTINGS], context = CONTEXT_LABELS[SETTINGS], stringsAsFactors = FALSE)
for (mn in metric_names) tab1[[mn]] <- sapply(SETTINGS, function(s) results[[s]]$metrics[[mn]])
tab1$n_samples_analyzed <- sapply(SETTINGS, function(s) table2[[s]]$n_samples_analyzed)
tab1$n_pairs_input <- sapply(SETTINGS, function(s) results[[s]]$drop[["n_pairs_input"]])
tab1$n_pairs_dropped <- tab1$n_pairs_input - tab1$n_pairs_analyzed
ci_df <- do.call(rbind, cis)
for (rt in unique(ci_df$rate)) {
  tab1[[paste0(rt, "_ci95_lo")]] <- sapply(SETTINGS, function(s) ci_df$ci95_lo_pct[ci_df$setting == s & ci_df$rate == rt])
  tab1[[paste0(rt, "_ci95_hi")]] <- sapply(SETTINGS, function(s) ci_df$ci95_hi_pct[ci_df$setting == s & ci_df$rate == rt])
}
tab1$eps <- EPS; tab1$abs_r_global_threshold <- THR; tab1$min_valid_contexts_per_pair <- K_MIN
rownames(tab1) <- NULL
write.csv(tab1, file.path(out_dir, "table1_metrics.csv"), row.names = FALSE)

metric_meta <- data.frame(
  column = c("m1_median_abs_r_global", "m1_pct_abs_r_global_ge_0.2",
             "m2_median_sd_rc", "m2_pct_mixed_sign_eps0.05", "m2_pct_mixed_sign_eps0.05_maxabs_rc_ge_0.2", "m2_median_k_valid_contexts",
             "m3_median_I2", "m3_I2_q25", "m3_I2_q75", "m3_pct_Q_fdr_lt_0.05",
             "m4_median_abs_r_RE", "m4_pct_sign_rRE_ne_rglobal_ge_0.2", "m4_pct_sign_rRE_ne_rglobal_all", "m4_median_abs_r_fixed",
             "m5_pct_resid_flip_all", "m5_pct_resid_flip_ge_0.2", "m5_pct_simpson_reversal_ge_0.2", "m5_pct_simpson_reversal_all",
             "n_pairs_analyzed", "n_pairs_abs_r_global_ge_0.2", "n_pairs_maxabs_rc_ge_0.2", "n_pairs_eligible_simpson_ge_0.2", "n_pairs_eligible_simpson_all", "n_contexts_retained", "n_samples_analyzed"),
  group = c(rep("1. Pooled correlation", 2), rep("2. Within-context distribution", 4), rep("3. Heterogeneity", 4), rep("4. Random-effects estimate", 4), rep("5. Between-vs-within diagnostics", 4), rep("Denominators", 7)),
  description = c("median |r_global| over analyzed pairs", "% of analyzed pairs with |r_global| >= 0.2",
                  "median over pairs of SD of r_c across valid contexts", "% mixed-sign (eps = 0.05) among analyzed pairs", "% mixed-sign (eps = 0.05) among pairs with max_c |r_c| >= 0.2", "median number of valid contexts per pair",
                  "median I2", "I2 first quartile", "I2 third quartile", "% of pairs with Cochran Q significant at BH FDR < 0.05",
                  "median |r_RE| (DerSimonian-Laird)", "% with sign(r_RE) != sign(r_global) among |r_global| >= 0.2", "% with sign(r_RE) != sign(r_global) among all analyzed pairs", "median |r_fixed| (fixed-effects Fisher-z pool)",
                  "% residualized sign flip among all analyzed pairs", "% residualized sign flip among |r_global| >= 0.2", "% Simpson reversal at |r_global| >= 0.2 (eps 0.05, ties excluded)", "% Simpson reversal among all eligible pairs (eps 0.05, ties excluded)",
                  "pairs analyzed (finite r_global, finite r_resid, >= 2 valid contexts)", "pairs with |r_global| >= 0.2", "pairs with max_c |r_c| >= 0.2", "pairs eligible for the reversal rate at |r_global| >= 0.2 (non-tie majority)", "pairs eligible for the all-pairs reversal rate", "contexts retained (>= minimum size)", "samples or cells in the analysis population"),
  denominator_column = c("n_pairs_analyzed", "n_pairs_analyzed", "n_pairs_analyzed", "n_pairs_analyzed", "n_pairs_maxabs_rc_ge_0.2", "n_pairs_analyzed",
                         "n_pairs_analyzed", "n_pairs_analyzed", "n_pairs_analyzed", "n_pairs_analyzed",
                         "n_pairs_analyzed", "n_pairs_abs_r_global_ge_0.2", "n_pairs_analyzed", "n_pairs_analyzed",
                         "n_pairs_analyzed", "n_pairs_abs_r_global_ge_0.2", "n_pairs_eligible_simpson_ge_0.2", "n_pairs_eligible_simpson_all",
                         rep("", 7)),
  stringsAsFactors = FALSE
)
metric_meta <- rbind(metric_meta, data.frame(
  column = c("m5_pct_nominal_flip_all", "m5_pct_nominal_flip_ge_0.2", "n_pairs_both_significant", "m5_pct_nominal_flip_among_both_significant"),
  group = c(rep("5. Between-vs-within diagnostics", 2), "Denominators", "5. Between-vs-within diagnostics"),
  description = c("% nominally significant flip (r_global and r_resid each two-sided p < 0.05, opposite signs) among all analyzed pairs",
                  "% nominally significant flip among |r_global| >= 0.2",
                  "pairs with r_global and r_resid both two-sided p < 0.05",
                  "% opposite signs among pairs with r_global and r_resid both p < 0.05"),
  denominator_column = c("n_pairs_analyzed", "n_pairs_abs_r_global_ge_0.2", "", "n_pairs_both_significant"),
  stringsAsFactors = FALSE))
long <- do.call(rbind, lapply(seq_len(nrow(metric_meta)), function(i) {
  data.frame(setting = SETTINGS, metric_group = metric_meta$group[i], metric = metric_meta$column[i],
             description = metric_meta$description[i], value = tab1[[metric_meta$column[i]]],
             denominator = if (metric_meta$denominator_column[i] != "") tab1[[metric_meta$denominator_column[i]]] else NA,
             csv = "table1_metrics.csv", csv_column = metric_meta$column[i], csv_row_setting = SETTINGS, stringsAsFactors = FALSE)
}))
write.csv(long, file.path(out_dir, "table1_metrics_long.csv"), row.names = FALSE)
write.csv(ci_df, file.path(out_dir, "bootstrap_cis_by_setting.csv"), row.names = FALSE)

ctx_df <- do.call(rbind, ctx_sizes); rownames(ctx_df) <- NULL
write.csv(ctx_df, file.path(out_dir, "context_sizes.csv"), row.names = FALSE)

ctx_stats <- function(s) {
  d <- ctx_df[ctx_df$setting == s & ctx_df$retained, ]
  c(min = min(d$n), median = stats::median(d$n), max = max(d$n))
}
tab2 <- do.call(rbind, lapply(SETTINGS, function(s) {
  t2 <- table2[[s]]; cs <- ctx_stats(s)
  data.frame(setting = s, setting_label = SETTING_LABELS[[s]], context_definition = CONTEXT_LABELS[[s]],
             data_source = t2$data_source, sample_unit = t2$sample_unit,
             n_total_in_dataset = t2$n_total_in_dataset, n_samples_analyzed = t2$n_samples_analyzed,
             n_features_x = t2$n_features_x, n_features_y = t2$n_features_y, feature_filter_rule = t2$feature_filter_rule,
             transform = t2$transform, n_contexts_total = t2$n_contexts_total, n_contexts_retained = t2$n_contexts_retained,
             min_context_size_rule = t2$min_context_size_rule,
             context_size_min = cs[["min"]], context_size_median = cs[["median"]], context_size_max = cs[["max"]],
             n_pairs_analyzed = t2$n_pairs_analyzed, pair_selection = t2$pair_selection,
             correlation_type = "Pearson", sign_tolerance_eps = EPS, min_valid_contexts_per_pair = K_MIN,
             notes = t2$notes, stringsAsFactors = FALSE)
}))
rownames(tab2) <- NULL
write.csv(tab2, file.path(out_dir, "table2_preprocessing.csv"), row.names = FALSE)

repro_df <- do.call(rbind, repro); rownames(repro_df) <- NULL
write.csv(repro_df, file.path(out_dir, "reproduction_checks.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# LaTeX tables (numbers read back from the CSVs just written)
# -----------------------------------------------------------------------------
t1 <- read.csv(file.path(out_dir, "table1_metrics.csv"), stringsAsFactors = FALSE, check.names = FALSE)
t2 <- read.csv(file.path(out_dir, "table2_preprocessing.csv"), stringsAsFactors = FALSE, check.names = FALSE)
cisv <- read.csv(file.path(out_dir, "bootstrap_cis_by_setting.csv"), stringsAsFactors = FALSE)
fpct <- function(x) sprintf("%.1f", x)
fnum <- function(x) sprintf("%.2f", x)
fint <- function(x) format(as.integer(round(x)), big.mark = ",", trim = TRUE)
tex_escape <- function(s) {
  s <- gsub("\\\\", "\\\\textbackslash{}", s)
  s <- gsub("([&%#_{}])", "\\\\\\1", s)
  s <- gsub(">=", "$\\\\geq$", s, fixed = FALSE)
  s <- gsub("<=", "$\\\\leq$", s, fixed = FALSE)
  s <- gsub(" x ", " $\\\\times$ ", s)
  s <- gsub("<", "$<$", s, fixed = TRUE)
  s
}
col_of <- function(col, s) t1[[col]][t1$setting == s]
ci_cell <- function(rate, s) {
  r <- cisv[cisv$setting == s & cisv$rate == rate, ]
  sprintf("%s [%s, %s]", fpct(r$estimate_pct), fpct(r$ci95_lo_pct), fpct(r$ci95_hi_pct))
}
head_labels <- c(TCGA = "TCGA", GTEx = "GTEx", PBMC = "PBMC", CITEseq = "CITE-seq", BRCA_PAM50 = "BRCA PAM50", BRCA_IHC = "BRCA IHC")
row_tex <- function(label, cells) paste0(label, " & ", paste(cells, collapse = " & "), " \\\\")
grp_tex <- function(g) sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", g)
L <- c(
  "% Generated by scripts/revision/uniform_table.R from results/revision/uniform_table/table1_metrics.csv",
  "% and bootstrap_cis_by_setting.csv. Requires \\usepackage{booktabs}. Bracketed values are pair-level",
  sprintf("%% percentile bootstrap 95%% CIs (B = %d). eps = %.2f; threshold on |r_global| = %.1f.", B_BOOT, EPS, THR),
  "\\begin{tabular}{l*{6}{r}}",
  "\\toprule",
  row_tex("Metric", head_labels[SETTINGS]),
  row_tex("Context", tex_escape(CONTEXT_LABELS[SETTINGS])),
  "\\midrule",
  grp_tex("1. Pooled correlation $r_{\\mathrm{global}}$"),
  row_tex("Median $|r_{\\mathrm{global}}|$", sapply(SETTINGS, function(s) fnum(col_of("m1_median_abs_r_global", s)))),
  row_tex("\\% pairs with $|r_{\\mathrm{global}}| \\geq 0.2$", sapply(SETTINGS, function(s) fpct(col_of("m1_pct_abs_r_global_ge_0.2", s)))),
  "\\midrule",
  grp_tex("2. Within-context distribution of $r_c$ ($\\varepsilon = 0.05$)"),
  row_tex("Median SD of $r_c$ across contexts", sapply(SETTINGS, function(s) fnum(col_of("m2_median_sd_rc", s)))),
  row_tex("\\% mixed-sign, all pairs", sapply(SETTINGS, function(s) fpct(col_of("m2_pct_mixed_sign_eps0.05", s)))),
  row_tex("\\% mixed-sign, $\\max_c |r_c| \\geq 0.2$", sapply(SETTINGS, function(s) fpct(col_of("m2_pct_mixed_sign_eps0.05_maxabs_rc_ge_0.2", s)))),
  "\\midrule",
  grp_tex("3. Heterogeneity (Cochran's $Q$, $I^2$)"),
  row_tex("Median $I^2$", sapply(SETTINGS, function(s) fnum(col_of("m3_median_I2", s)))),
  row_tex("$I^2$ IQR", sapply(SETTINGS, function(s) sprintf("%s--%s", fnum(col_of("m3_I2_q25", s)), fnum(col_of("m3_I2_q75", s))))),
  row_tex("\\% pairs with $Q$ significant (FDR $< 0.05$)", sapply(SETTINGS, function(s) fpct(col_of("m3_pct_Q_fdr_lt_0.05", s)))),
  "\\midrule",
  grp_tex("4. Random-effects estimate $r_{\\mathrm{RE}}$"),
  row_tex("Median $|r_{\\mathrm{RE}}|$", sapply(SETTINGS, function(s) fnum(col_of("m4_median_abs_r_RE", s)))),
  row_tex("\\% sign$(r_{\\mathrm{RE}}) \\neq$ sign$(r_{\\mathrm{global}})$, $|r_{\\mathrm{global}}| \\geq 0.2$", sapply(SETTINGS, function(s) fpct(col_of("m4_pct_sign_rRE_ne_rglobal_ge_0.2", s)))),
  "\\midrule",
  grp_tex("5. Between-vs-within diagnostics"),
  row_tex("\\% residualized sign flip, all pairs", sapply(SETTINGS, function(s) ci_cell("resid_flip_all", s))),
  row_tex("\\% residualized sign flip, $|r_{\\mathrm{global}}| \\geq 0.2$", sapply(SETTINGS, function(s) ci_cell("resid_flip_ge_0.2", s))),
  row_tex("\\% Simpson reversal, $|r_{\\mathrm{global}}| \\geq 0.2$", sapply(SETTINGS, function(s) ci_cell("simpson_reversal_ge_0.2", s))),
  row_tex("\\% nominally significant flip, all pairs", sapply(SETTINGS, function(s) ci_cell("nominal_flip_all", s))),
  row_tex("\\% nominally significant flip, $|r_{\\mathrm{global}}| \\geq 0.2$", sapply(SETTINGS, function(s) ci_cell("nominal_flip_ge_0.2", s))),
  "\\midrule",
  grp_tex("Denominators"),
  row_tex("Samples or cells analyzed", sapply(SETTINGS, function(s) fint(col_of("n_samples_analyzed", s)))),
  row_tex("Contexts retained", sapply(SETTINGS, function(s) fint(col_of("n_contexts_retained", s)))),
  row_tex("Pairs analyzed", sapply(SETTINGS, function(s) fint(col_of("n_pairs_analyzed", s)))),
  row_tex("Pairs with $|r_{\\mathrm{global}}| \\geq 0.2$", sapply(SETTINGS, function(s) fint(col_of("n_pairs_abs_r_global_ge_0.2", s)))),
  row_tex("Pairs eligible for reversal at $\\geq 0.2$", sapply(SETTINGS, function(s) fint(col_of("n_pairs_eligible_simpson_ge_0.2", s)))),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(L, file.path(out_dir, "table1.tex"))

t2_of <- function(col, s) t2[[col]][t2$setting == s]
short_source <- c(
  TCGA = "TCGA (FireBrowse, 2016-01-28 release)", GTEx = "GTEx v8 (UCSC Xena matrix; SMTSD annotations)",
  PBMC = "10x PBMC 3k (single donor)", CITEseq = "10x PBMC 10k TotalSeq-B (pbmc\\_10k\\_protein\\_v3)",
  BRCA_PAM50 = "TCGA BRCA primary tumors; Xena clinical matrix", BRCA_IHC = "TCGA BRCA primary tumors; Xena clinical matrix"
)
short_features <- c(
  TCGA = "top 5,000 mRNAs $\\times$ top 200 miRNAs by variance (all samples)",
  GTEx = "1,500 most variable genes (within the 20 tissues)",
  PBMC = "1,655 genes: 2,000 HVGs + markers, detected in $\\geq 1\\%$ of cells, cap 1,800",
  CITEseq = "14 ADTs (isotype controls removed) $\\times$ 1,055 genes (1,500 HVGs + markers, detected $\\geq 1\\%$, cap 1,300)",
  BRCA_PAM50 = "top 5,000 mRNAs $\\times$ top 200 miRNAs by variance within BRCA tumors",
  BRCA_IHC = "top 5,000 mRNAs $\\times$ top 200 miRNAs by variance within BRCA tumors"
)
short_transform <- c(
  TCGA = "mRNA $\\log_2(10 + x)$; miRNA $\\log_2(1 + x)$", GTEx = "$\\log_2(\\mathrm{TPM} + 0.001)$",
  PBMC = "QC ($\\geq 200$ genes, $\\leq 20\\%$ mito); LogNormalize", CITEseq = "QC as PBMC; RNA LogNormalize; ADT CLR",
  BRCA_PAM50 = "as TCGA", BRCA_IHC = "as TCGA"
)
short_pairs <- c(
  TCGA = "all $10^6$", GTEx = "200,000 sampled (seed 1)", PBMC = "200,000 sampled (seed 100)",
  CITEseq = "all 14,770", BRCA_PAM50 = "all $10^6$", BRCA_IHC = "all $10^6$"
)
L2 <- c(
  "% Generated by scripts/revision/uniform_table.R from results/revision/uniform_table/table2_preprocessing.csv.",
  "% Requires \\usepackage{booktabs,array}. Full-length descriptions are in the CSV (columns data_source,",
  "% feature_filter_rule, transform, pair_selection, notes).",
  "\\begin{tabular}{>{\\raggedright\\arraybackslash}p{3.0cm}*{6}{>{\\raggedright\\arraybackslash}p{2.3cm}}}",
  "\\toprule",
  row_tex("", head_labels[SETTINGS]),
  "\\midrule",
  row_tex("Context", tex_escape(CONTEXT_LABELS[SETTINGS])),
  row_tex("Data source", short_source[SETTINGS]),
  row_tex("Unit", sapply(SETTINGS, function(s) tex_escape(t2_of("sample_unit", s)))),
  row_tex("N in dataset (after QC / with label)", sapply(SETTINGS, function(s) fint(t2_of("n_total_in_dataset", s)))),
  row_tex("N analyzed (in retained contexts)", sapply(SETTINGS, function(s) fint(t2_of("n_samples_analyzed", s)))),
  row_tex("Features after filtering", short_features[SETTINGS]),
  row_tex("Transform", short_transform[SETTINGS]),
  row_tex("Contexts (total / retained)", sapply(SETTINGS, function(s) sprintf("%s / %s", fint(t2_of("n_contexts_total", s)), fint(t2_of("n_contexts_retained", s))))),
  row_tex("Minimum context size", sapply(SETTINGS, function(s) tex_escape(t2_of("min_context_size_rule", s)))),
  row_tex("Context size min / median / max", sapply(SETTINGS, function(s) sprintf("%s / %s / %s", fint(t2_of("context_size_min", s)), fint(t2_of("context_size_median", s)), fint(t2_of("context_size_max", s))))),
  row_tex("Pairs analyzed", sapply(SETTINGS, function(s) fint(t2_of("n_pairs_analyzed", s)))),
  row_tex("Pair selection", short_pairs[SETTINGS]),
  row_tex("Correlation", sapply(SETTINGS, function(s) t2_of("correlation_type", s))),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(L2, file.path(out_dir, "table2.tex"))

# -----------------------------------------------------------------------------
# Runtime and session information
# -----------------------------------------------------------------------------
toc("assemble_outputs", t0)
total_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
rt <- data.frame(step = names(timings), seconds = round(unlist(timings), 1), stringsAsFactors = FALSE)
rt <- rbind(rt, data.frame(step = "total", seconds = round(total_secs, 1)))
write.csv(rt, file.path(out_dir, "runtime.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
