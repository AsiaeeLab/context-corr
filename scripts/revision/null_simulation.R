#!/usr/bin/env Rscript

# Null calibration of mixed-sign, Simpson-reversal, and I^2 statistics under
# NO true within-context heterogeneity, plus a downsampling experiment on the
# real TCGA data. Context sizes and imbalance are varied to quantify the
# contribution of estimation noise to each statistic.
#
# Run from the repo root:  Rscript scripts/revision/null_simulation.R
# Quick smoke test (small sizes, scratch output dir):
#   NULLSIM_QUICK=1 NULLSIM_OUT=/some/dir Rscript scripts/revision/null_simulation.R
#
# Conventions are copied from the paper's scripts (scripts/robustness_signflip.R,
# scripts/heterogeneity_I2_FDR.R, R/contextcorr.R):
#   - high-variance domain: top 5,000 mRNAs x top 200 miRNAs by global variance
#   - a within-context r_c is valid when finite and |r_c| < 1
#   - a pair is evaluable when at least min_k = 10 contexts are valid
#   - sign tolerance eps: s_c = +1 if r_c >= eps, -1 if r_c <= -eps, else 0
#   - mixed-sign: both +1 and -1 present among contexts
#   - Simpson reversal: sign(r_global) != sign(sum_c s_c), ties and r_global = 0
#     excluded from the denominator; reported among |r_global| >= 0.2
#   - Q with weights N_c - 3 on Fisher z; I^2 = max(0, (Q - (k - 1)) / Q)
#
# Outputs (results/revision/null_simulation/):
#   null_rates.csv                one row per simulated scenario
#   null_pairs.rds                per-pair statistics for every scenario
#   meanshift_empirical.csv       empirical mean-shift distribution summaries
#   meanshift_sampled_pairs.rds   the 20,000 sampled real pairs' cohort means
#   tcga_reference.csv            observed TCGA rates recomputed on the full
#                                 1e6-pair domain with the same code path
#   downsampling_tcga.csv         mean and 2.5/97.5% range over replicates per n
#   downsampling_tcga_replicates.csv   one row per (n, replicate)
#   Figure_NullCalibration.pdf
#   runtime.txt

t_start <- Sys.time()
section_times <- list()
tic <- function() Sys.time()
toc <- function(label, t0) {
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  section_times[[label]] <<- dt
  message(sprintf("[%s] %.1f s", label, dt))
}

source("00-paths.R")
source(file.path("R", "contextcorr.R"))

QUICK <- identical(Sys.getenv("NULLSIM_QUICK"), "1")
out_dir <- Sys.getenv("NULLSIM_OUT", unset = file.path("results", "revision", "null_simulation"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_cores <- max(1L, min(20L, parallel::detectCores() - 2L))

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
G <- 5000                       # top mRNAs by global variance (paper domain)
M <- 200                        # top miRNAs by global variance (paper domain)
min_n_per_cohort <- 20          # paper convention (no TCGA cohort is below 45)
min_k <- 10                     # paper convention
eps_main <- 0.05                # sign tolerance for headline rates
eps_grid <- c(0, 0.05, 0.10)    # paper's eps sensitivity grid
t_main <- 0.2                   # |r_global| threshold for reversal rate

rhos <- c(0, 0.1, 0.3)
balanced_ns <- c(20, 36, 50, 100, 300)
pairs_per_scenario <- if (QUICK) 2000L else 20000L
chunk_size <- 2000L
down_ns <- c(20, 36, 50, 100)
down_reps <- if (QUICK) 2L else 20L
down_pairs <- if (QUICK) 20000L else 100000L

seed_pairs_B <- 20260903        # which real pairs supply the mean shifts
seed_down_pairs <- 20260904     # the fixed 100,000-pair TCGA subsample
seed_sim_base <- 10000L         # per-(design, chunk) seeds: base * design + chunk
seed_down_base <- 300000L       # per-(replicate, n) seeds: base + 1000 * rep + n

# ---------------------------------------------------------------------------
# Load matched TCGA data
# ---------------------------------------------------------------------------
t0 <- tic()
f1 <- file.path(paths$mirtcga, "matchedNoGBM.Rda")
f2 <- file.path(paths$mirtcga, "matchedNoGBM.rda")
if (file.exists(f1)) {
  load(f1)
} else if (file.exists(f2)) {
  load(f2)
} else {
  stop("Cannot find matchedNoGBM.Rda or matchedNoGBM.rda in ", paths$mirtcga)
}

cohorts <- levels(cancerType)
idx_by_ct <- lapply(cohorts, function(ct) which(cancerType == ct))
names(idx_by_ct) <- cohorts
N_c <- vapply(idx_by_ct, length, integer(1))
k_ctx <- length(cohorts)
stopifnot(all(N_c >= min_n_per_cohort))

gene_var <- apply(mRNAdata, 1, var)
mir_var <- apply(miRdata, 1, var)
genes <- names(sort(gene_var, decreasing = TRUE))[1:G]
mirs <- names(sort(mir_var, decreasing = TRUE))[1:M]
mRNA_sub <- mRNAdata[genes, , drop = FALSE]
miR_sub <- miRdata[mirs, , drop = FALSE]
rm(mRNAdata, miRdata, gene_var, mir_var)
invisible(gc())
toc("load_data", t0)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# Pearson correlations between every row of A and every row of B, using the
# columns idx (same code path as scripts/robustness_signflip.R).
cor_rows <- function(A, B, idx) {
  A <- A[, idx, drop = FALSE]
  B <- B[, idx, drop = FALSE]
  scaled_B <- scale(t(B)) / (ncol(B) - 1)
  std_A <- t(scale(t(A)))
  std_A %*% scaled_B
}

# Column-wise Pearson correlation between X[idx, ] and Y[idx, ] (one value per
# column = one simulated pair). Used by the simulation.
col_cor <- function(X, Y, idx) {
  Xc <- X[idx, , drop = FALSE]
  Yc <- Y[idx, , drop = FALSE]
  n <- length(idx)
  mx <- colMeans(Xc)
  my <- colMeans(Yc)
  sxx <- colSums(Xc * Xc) - n * mx^2
  syy <- colSums(Yc * Yc) - n * my^2
  sxy <- colSums(Xc * Yc) - n * mx * my
  sxy / sqrt(sxx * syy)
}

# Per-pair statistics from a vector of pooled correlations and a P x k matrix of
# within-context correlations (NA allowed). n_c gives the context sizes used as
# Q weights (N_c - 3).
pair_metrics <- function(r_global, R, n_c, eps_grid = c(0, 0.05, 0.10)) {
  stopifnot(length(r_global) == nrow(R), length(n_c) == ncol(R))
  P <- nrow(R)
  k <- ncol(R)
  valid <- is.finite(R) & abs(R) < 1
  k_valid <- rowSums(valid)

  max_abs <- rep(0, P)
  for (j in seq_len(k)) {
    a <- abs(R[, j])
    a[!valid[, j]] <- 0
    max_abs <- pmax(max_abs, a)
  }

  W <- matrix(n_c - 3, nrow = P, ncol = k, byrow = TRUE)
  W[!valid] <- 0
  Z <- contextcorr_fisher_z(R)
  Z[!valid] <- 0
  sum_w <- rowSums(W)
  sum_wz <- rowSums(W * Z)
  sum_wz2 <- rowSums(W * Z^2)
  Q <- sum_wz2 - sum_wz^2 / pmax(sum_w, 1e-12)
  df <- pmax(k_valid - 1L, 0L)
  ok <- df >= 1L & is.finite(Q) & Q >= 0
  I2 <- rep(NA_real_, P)
  I2[ok] <- pmax(0, (Q[ok] - df[ok]) / pmax(Q[ok], 1e-12))
  p_Q <- rep(NA_real_, P)
  p_Q[ok] <- stats::pchisq(Q[ok], df = df[ok], lower.tail = FALSE)

  out <- data.frame(r_global = as.numeric(r_global), k_valid = k_valid,
                    max_abs = max_abs, Q = Q, df = df, I2 = I2, p_Q = p_Q)
  for (eps in eps_grid) {
    tag <- gsub("\\.", "p", formatC(eps, format = "f", digits = 2))
    out[[paste0("n_pos_eps", tag)]] <- rowSums(valid & (R >= eps))
    out[[paste0("n_neg_eps", tag)]] <- rowSums(valid & (R <= -eps))
  }
  out
}

# Scenario-level rates from a pair_metrics data frame.
rate_summary <- function(pm, eps = 0.05, t = 0.2, min_k = 10) {
  tagf <- function(e) gsub("\\.", "p", formatC(e, format = "f", digits = 2))
  valid <- pm$k_valid >= min_k & is.finite(pm$r_global)
  rg <- pm$r_global
  sg <- sign(rg)

  one_eps <- function(e) {
    pos <- pm[[paste0("n_pos_eps", tagf(e))]]
    neg <- pm[[paste0("n_neg_eps", tagf(e))]]
    mixed <- pos > 0 & neg > 0
    maj <- sign(pos - neg)
    elig_t <- valid & abs(rg) >= t & maj != 0 & sg != 0
    elig_0 <- valid & maj != 0 & sg != 0
    list(
      frac_mixed = mean(mixed[valid]),
      frac_mixed_given_maxabs = if (any(valid & pm$max_abs >= t)) mean(mixed[valid & pm$max_abs >= t]) else NA_real_,
      n_elig_simpson = sum(elig_t),
      frac_simpson = if (sum(elig_t) > 0) mean(sg[elig_t] != maj[elig_t]) else NA_real_,
      n_elig_simpson_t0 = sum(elig_0),
      frac_simpson_t0 = if (sum(elig_0) > 0) mean(sg[elig_0] != maj[elig_0]) else NA_real_,
      frac_tie = mean((maj == 0)[valid])
    )
  }
  m <- one_eps(eps)
  m0 <- one_eps(0)
  m10 <- one_eps(0.10)

  binom_se <- function(p, n) if (is.finite(p) && n > 0) sqrt(p * (1 - p) / n) else NA_real_

  data.frame(
    n_pairs = nrow(pm),
    n_valid = sum(valid),
    frac_mixed_eps0p05 = m$frac_mixed,
    se_frac_mixed_eps0p05 = binom_se(m$frac_mixed, sum(valid)),
    se_frac_simpson_eps0p05_t0p2 = binom_se(m$frac_simpson, m$n_elig_simpson),
    frac_maxabs_ge0p2 = mean((pm$max_abs >= t)[valid]),
    frac_mixed_eps0p05_given_maxabs_ge0p2 = m$frac_mixed_given_maxabs,
    frac_rglobal_ge0p2 = mean((abs(rg) >= t)[valid]),
    median_abs_rglobal = stats::median(abs(rg[valid])),
    n_elig_simpson_t0p2 = m$n_elig_simpson,
    frac_simpson_eps0p05_t0p2 = m$frac_simpson,
    n_elig_simpson_t0 = m$n_elig_simpson_t0,
    frac_simpson_eps0p05_t0 = m$frac_simpson_t0,
    frac_majority_tie_eps0p05 = m$frac_tie,
    median_I2 = stats::median(pm$I2[valid], na.rm = TRUE),
    mean_I2 = mean(pm$I2[valid], na.rm = TRUE),
    I2_q25 = as.numeric(stats::quantile(pm$I2[valid], 0.25, na.rm = TRUE)),
    I2_q75 = as.numeric(stats::quantile(pm$I2[valid], 0.75, na.rm = TRUE)),
    frac_I2_gt0p75 = mean(pm$I2[valid] > 0.75, na.rm = TRUE),
    frac_Qp_lt0p05 = mean(pm$p_Q[valid] < 0.05, na.rm = TRUE),
    frac_mixed_eps0 = m0$frac_mixed,
    frac_simpson_eps0_t0p2 = m0$frac_simpson,
    frac_mixed_eps0p10 = m10$frac_mixed,
    frac_simpson_eps0p10_t0p2 = m10$frac_simpson,
    stringsAsFactors = FALSE
  )
}

# Within-cohort and pooled correlations on the real high-variance domain for a
# given list of per-cohort sample indices; returns the P x k matrix R and the
# pooled vector r_global restricted to the pairs in `pair_idx` (column-major
# indices into the G x M matrix).
tcga_pair_correlations <- function(idx_list, pair_idx) {
  R <- matrix(NA_real_, nrow = length(pair_idx), ncol = length(idx_list))
  for (j in seq_along(idx_list)) {
    idx <- idx_list[[j]]
    if (length(idx) < min_n_per_cohort) next
    r_ct <- cor_rows(mRNA_sub, miR_sub, idx)
    R[, j] <- r_ct[pair_idx]
  }
  all_idx <- sort(unlist(idx_list, use.names = FALSE))
  r_glob <- cor_rows(mRNA_sub, miR_sub, all_idx)[pair_idx]
  list(R = R, r_global = r_glob, n_c = vapply(idx_list, length, integer(1)))
}

# ---------------------------------------------------------------------------
# Part 1: empirical between-cohort mean-shift distribution (condition B)
# ---------------------------------------------------------------------------
t0 <- tic()
E <- rbind(mRNA_sub, miR_sub)                       # (G + M) x N
feature_type <- c(rep("mRNA", G), rep("miRNA", M))
N_total <- ncol(E)
grand_mean <- rowMeans(E)
cohort_mean <- matrix(NA_real_, nrow = nrow(E), ncol = k_ctx, dimnames = list(rownames(E), cohorts))
within_ss <- matrix(NA_real_, nrow = nrow(E), ncol = k_ctx, dimnames = list(rownames(E), cohorts))
for (j in seq_len(k_ctx)) {
  Ec <- E[, idx_by_ct[[j]], drop = FALSE]
  mc <- rowMeans(Ec)
  cohort_mean[, j] <- mc
  within_ss[, j] <- rowSums((Ec - mc)^2)
}
pooled_within_sd <- sqrt(rowSums(within_ss) / (N_total - k_ctx))   # residual SD after cohort means
total_ss <- rowSums((E - grand_mean)^2)
between_ss <- as.numeric((cohort_mean - grand_mean)^2 %*% N_c)
eta2 <- between_ss / total_ss                                       # fraction of variance between cohorts
D <- (cohort_mean - grand_mean) / pooled_within_sd                  # standardized cohort-mean differences
between_sd_std <- sqrt(as.numeric(D^2 %*% N_c) / N_total)           # N_c-weighted SD of D, within-SD units
n_zero_within <- sum(!is.finite(D[, 1]))
if (n_zero_within > 0) stop("Features with zero pooled within-cohort SD: ", n_zero_within)
rm(E, within_ss)
invisible(gc())

# Sample the real pairs whose cohort means will be used in condition B. Same
# pairs for every scenario; pair j of every scenario uses real pair j.
set.seed(seed_pairs_B)
pairB_idx <- sample.int(G * M, pairs_per_scenario)
pairB_gene <- (pairB_idx - 1L) %% G + 1L
pairB_mir <- (pairB_idx - 1L) %/% G + 1L
D_gene <- D[pairB_gene, , drop = FALSE]                  # P x k
D_mir <- D[G + pairB_mir, , drop = FALSE]                # P x k

wcor <- function(A, B, w) {
  w <- w / sum(w)
  ma <- as.numeric(A %*% w)
  mb <- as.numeric(B %*% w)
  Ac <- A - ma
  Bc <- B - mb
  cab <- as.numeric((Ac * Bc) %*% w)
  vaa <- as.numeric((Ac * Ac) %*% w)
  vbb <- as.numeric((Bc * Bc) %*% w)
  list(cor = cab / sqrt(vaa * vbb), cov = cab, va = vaa, vb = vbb)
}
# Implied pooled correlation when rho = 0 and within-context SD = 1:
# r_between = Cov_B / sqrt((1 + Var_B,x)(1 + Var_B,y)), with N_c-weighted moments.
implied_r_between <- function(w) {
  z <- wcor(D_gene, D_mir, w)
  z$cov / sqrt((1 + z$va) * (1 + z$vb))
}
cm_cor_w <- wcor(D_gene, D_mir, N_c)$cor
cm_cor_u <- wcor(D_gene, D_mir, rep(1, k_ctx))$cor
r_between_tcga <- implied_r_between(N_c)
r_between_bal <- implied_r_between(rep(1, k_ctx))

qsum <- function(x, block, feature_type = NA, cohort = NA, n_cohort = NA) {
  x <- x[is.finite(x)]
  q <- stats::quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
  data.frame(block = block, feature_type = feature_type, cohort = cohort, n_cohort = n_cohort,
             n = length(x), mean = mean(x), sd = stats::sd(x), min = min(x),
             q025 = q[1], q25 = q[2], median = q[3], q75 = q[4], q975 = q[5], max = max(x),
             stringsAsFactors = FALSE)
}
ms_rows <- list()
for (ft in c("mRNA", "miRNA")) {
  sel <- feature_type == ft
  for (j in seq_len(k_ctx)) {
    ms_rows[[length(ms_rows) + 1]] <- qsum(D[sel, j], "std_mean_diff_by_cohort", ft, cohorts[j], N_c[j])
  }
  ms_rows[[length(ms_rows) + 1]] <- qsum(as.numeric(D[sel, ]), "std_mean_diff_all_cohorts", ft)
  ms_rows[[length(ms_rows) + 1]] <- qsum(as.numeric(abs(D[sel, ])), "abs_std_mean_diff_all_cohorts", ft)
  ms_rows[[length(ms_rows) + 1]] <- qsum(between_sd_std[sel], "between_cohort_sd_in_within_sd_units", ft)
  ms_rows[[length(ms_rows) + 1]] <- qsum(eta2[sel], "eta2_between_cohort_variance_fraction", ft)
  ms_rows[[length(ms_rows) + 1]] <- qsum(pooled_within_sd[sel], "pooled_within_cohort_sd_log2_units", ft)
}
ms_rows[[length(ms_rows) + 1]] <- qsum(cm_cor_w, "sampled_pair_cohort_mean_cor_Nc_weighted", "mRNA-miRNA pair")
ms_rows[[length(ms_rows) + 1]] <- qsum(cm_cor_u, "sampled_pair_cohort_mean_cor_unweighted", "mRNA-miRNA pair")
ms_rows[[length(ms_rows) + 1]] <- qsum(r_between_tcga, "sampled_pair_implied_r_between_rho0_tcga_sizes", "mRNA-miRNA pair")
ms_rows[[length(ms_rows) + 1]] <- qsum(r_between_bal, "sampled_pair_implied_r_between_rho0_balanced", "mRNA-miRNA pair")
ms_rows[[length(ms_rows) + 1]] <- qsum(as.numeric(abs(r_between_tcga) >= t_main), "sampled_pair_frac_abs_implied_r_between_ge0p2_tcga_sizes", "mRNA-miRNA pair")
ms_rows[[length(ms_rows) + 1]] <- qsum(as.numeric(abs(r_between_bal) >= t_main), "sampled_pair_frac_abs_implied_r_between_ge0p2_balanced", "mRNA-miRNA pair")
ms_rows[[length(ms_rows) + 1]] <- qsum(as.numeric(N_c), "tcga_cohort_sizes", "cohort")
meanshift_df <- do.call(rbind, ms_rows)
utils::write.csv(meanshift_df, file.path(out_dir, "meanshift_empirical.csv"), row.names = FALSE)
saveRDS(list(genes = genes[pairB_gene], mirs = mirs[pairB_mir], D_gene = D_gene, D_mir = D_mir,
             cohorts = cohorts, N_c = N_c, pooled_within_sd_gene = pooled_within_sd[pairB_gene],
             pooled_within_sd_mir = pooled_within_sd[G + pairB_mir],
             cohort_mean_cor_weighted = cm_cor_w, cohort_mean_cor_unweighted = cm_cor_u,
             implied_r_between_tcga_sizes = r_between_tcga, implied_r_between_balanced = r_between_bal,
             seed = seed_pairs_B),
        file.path(out_dir, "meanshift_sampled_pairs.rds"))
toc("meanshift_empirical", t0)

# ---------------------------------------------------------------------------
# Part 2: observed TCGA reference on the full 1e6-pair domain (same code path)
# ---------------------------------------------------------------------------
t0 <- tic()
ref <- tcga_pair_correlations(idx_by_ct, seq_len(G * M))
pm_ref <- pair_metrics(ref$r_global, ref$R, ref$n_c, eps_grid)
ref_summary <- rate_summary(pm_ref, eps = eps_main, t = t_main, min_k = min_k)
ref_long <- data.frame(metric = names(ref_summary), value = as.numeric(ref_summary[1, ]),
                       domain = sprintf("top %d mRNAs x top %d miRNAs, all %d samples, %d cohorts", G, M, N_total, k_ctx),
                       stringsAsFactors = FALSE)
utils::write.csv(ref_long, file.path(out_dir, "tcga_reference.csv"), row.names = FALSE)
message(sprintf("TCGA reference: mixed %.4f | mixed given maxabs>=0.2 %.4f | simpson(|rg|>=0.2) %.4f | frac |rg|>=0.2 %.4f | median I2 %.4f",
                ref_summary$frac_mixed_eps0p05, ref_summary$frac_mixed_eps0p05_given_maxabs_ge0p2,
                ref_summary$frac_simpson_eps0p05_t0p2, ref_summary$frac_rglobal_ge0p2, ref_summary$median_I2))
rm(ref, pm_ref)
invisible(gc())
toc("tcga_reference", t0)

# ---------------------------------------------------------------------------
# Part 3: null simulation
# ---------------------------------------------------------------------------
t0 <- tic()
designs <- c(lapply(balanced_ns, function(n) list(name = "balanced", n_per_context = n, sizes = rep(as.integer(n), k_ctx))),
             list(list(name = "tcga_imbalanced", n_per_context = NA_real_, sizes = as.integer(N_c))))
n_chunks <- as.integer(ceiling(pairs_per_scenario / chunk_size))
tasks <- expand.grid(design = seq_along(designs), chunk = seq_len(n_chunks))

simulate_task <- function(d, j) {
  des <- designs[[d]]
  sizes <- des$sizes
  ctx_of_row <- rep(seq_len(k_ctx), times = sizes)
  idx_list <- split(seq_along(ctx_of_row), ctx_of_row)
  N <- length(ctx_of_row)
  p_from <- (j - 1L) * chunk_size + 1L
  p_to <- min(j * chunk_size, pairs_per_scenario)
  pr <- p_from:p_to
  P <- length(pr)

  set.seed(seed_sim_base * d + j)
  E1 <- matrix(stats::rnorm(N * P), nrow = N, ncol = P)
  E2 <- matrix(stats::rnorm(N * P), nrow = N, ncol = P)
  MUx <- t(D_gene[pr, , drop = FALSE])[ctx_of_row, , drop = FALSE]   # N x P context means
  MUy <- t(D_mir[pr, , drop = FALSE])[ctx_of_row, , drop = FALSE]

  out <- list()
  for (rho in rhos) {
    Y0 <- rho * E1 + sqrt(1 - rho^2) * E2
    for (shift in c("none", "tcga_matched")) {
      if (shift == "none") {
        X <- E1
        Y <- Y0
      } else {
        X <- E1 + MUx
        Y <- Y0 + MUy
      }
      R <- matrix(NA_real_, nrow = P, ncol = k_ctx)
      for (cc in seq_len(k_ctx)) R[, cc] <- col_cor(X, Y, idx_list[[cc]])
      r_glob <- col_cor(X, Y, seq_len(N))
      pm <- pair_metrics(r_glob, R, sizes, eps_grid)
      pm$pair <- pr
      pm$design <- des$name
      pm$n_per_context <- des$n_per_context
      pm$rho <- rho
      pm$mean_shift <- shift
      out[[paste(rho, shift)]] <- pm
    }
  }
  out
}

sim_res <- parallel::mclapply(seq_len(nrow(tasks)), function(i) simulate_task(tasks$design[i], tasks$chunk[i]),
                              mc.cores = n_cores, mc.preschedule = FALSE)
failed <- vapply(sim_res, inherits, logical(1), what = "try-error")
if (any(failed)) stop("Simulation tasks failed: ", paste(which(failed), collapse = ", "))

null_pairs <- do.call(rbind, unlist(sim_res, recursive = FALSE))
rownames(null_pairs) <- NULL
rm(sim_res)
invisible(gc())
saveRDS(null_pairs, file.path(out_dir, "null_pairs.rds"))

scen <- unique(null_pairs[, c("design", "n_per_context", "rho", "mean_shift")])
null_rates <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i) {
  s <- scen[i, ]
  sel <- null_pairs$design == s$design & (is.na(s$n_per_context) | null_pairs$n_per_context == s$n_per_context) &
    null_pairs$rho == s$rho & null_pairs$mean_shift == s$mean_shift
  if (is.na(s$n_per_context)) sel <- sel & is.na(null_pairs$n_per_context)
  rs <- rate_summary(null_pairs[sel, ], eps = eps_main, t = t_main, min_k = min_k)
  sizes <- if (s$design == "balanced") rep(s$n_per_context, k_ctx) else N_c
  cbind(data.frame(scenario = i, design = s$design, mean_shift = s$mean_shift, rho = s$rho,
                   k = k_ctx, n_per_context = s$n_per_context, n_min = min(sizes), n_median = stats::median(sizes),
                   n_max = max(sizes), n_total = sum(sizes), eps = eps_main, threshold = t_main, min_k = min_k,
                   stringsAsFactors = FALSE), rs)
}))
null_rates <- null_rates[order(null_rates$mean_shift, null_rates$rho, null_rates$design != "balanced", null_rates$n_per_context), ]
null_rates$scenario <- seq_len(nrow(null_rates))
rownames(null_rates) <- NULL
utils::write.csv(null_rates, file.path(out_dir, "null_rates.csv"), row.names = FALSE)
toc("null_simulation", t0)

# ---------------------------------------------------------------------------
# Part 4: downsampling on real TCGA
# ---------------------------------------------------------------------------
t0 <- tic()
set.seed(seed_down_pairs)
down_pair_idx <- sort(sample.int(G * M, down_pairs))

full_res <- tcga_pair_correlations(idx_by_ct, down_pair_idx)
full_summary <- rate_summary(pair_metrics(full_res$r_global, full_res$R, full_res$n_c, eps_grid),
                             eps = eps_main, t = t_main, min_k = min_k)
full_row <- cbind(data.frame(n_target = NA_real_, replicate = NA_integer_, n_total = N_total,
                             n_cohorts_subsampled = 0L, stringsAsFactors = FALSE), full_summary)
rm(full_res)

down_grid <- expand.grid(replicate = seq_len(down_reps), n_target = down_ns)
down_one <- function(rep_i, n_target) {
  set.seed(seed_down_base + 1000L * rep_i + as.integer(n_target))
  idx_list <- lapply(idx_by_ct, function(idx) if (length(idx) > n_target) sort(sample(idx, n_target)) else idx)
  res <- tcga_pair_correlations(idx_list, down_pair_idx)
  rs <- rate_summary(pair_metrics(res$r_global, res$R, res$n_c, eps_grid), eps = eps_main, t = t_main, min_k = min_k)
  cbind(data.frame(n_target = n_target, replicate = rep_i, n_total = sum(res$n_c),
                   n_cohorts_subsampled = sum(N_c > n_target), stringsAsFactors = FALSE), rs)
}
down_list <- parallel::mclapply(seq_len(nrow(down_grid)), function(i) down_one(down_grid$replicate[i], down_grid$n_target[i]),
                                mc.cores = n_cores, mc.preschedule = FALSE)
failed <- vapply(down_list, inherits, logical(1), what = "try-error")
if (any(failed)) stop("Downsampling tasks failed: ", paste(which(failed), collapse = ", "))
down_reps_df <- rbind(full_row, do.call(rbind, down_list))
down_reps_df$eps <- eps_main
down_reps_df$threshold <- t_main
down_reps_df$min_k <- min_k
down_reps_df$n_pairs_sampled <- down_pairs
rownames(down_reps_df) <- NULL
utils::write.csv(down_reps_df, file.path(out_dir, "downsampling_tcga_replicates.csv"), row.names = FALSE)

metric_cols <- c("frac_mixed_eps0p05", "frac_maxabs_ge0p2", "frac_mixed_eps0p05_given_maxabs_ge0p2",
                 "frac_rglobal_ge0p2", "median_abs_rglobal", "n_elig_simpson_t0p2", "frac_simpson_eps0p05_t0p2",
                 "frac_simpson_eps0p05_t0", "frac_majority_tie_eps0p05", "median_I2", "frac_I2_gt0p75",
                 "frac_Qp_lt0p05", "n_valid", "frac_mixed_eps0", "frac_simpson_eps0_t0p2",
                 "frac_mixed_eps0p10", "frac_simpson_eps0p10_t0p2")
down_summary <- do.call(rbind, lapply(c(down_ns, NA), function(n) {
  sub <- if (is.na(n)) down_reps_df[is.na(down_reps_df$n_target), ] else down_reps_df[down_reps_df$n_target %in% n, ]
  do.call(rbind, lapply(metric_cols, function(mc) {
    x <- sub[[mc]]
    data.frame(n_target = if (is.na(n)) "all" else as.character(n), n_total_mean = mean(sub$n_total),
               n_cohorts_subsampled = sub$n_cohorts_subsampled[1], n_replicates = nrow(sub), metric = mc,
               mean = mean(x, na.rm = TRUE),
               q025 = if (nrow(sub) > 1) as.numeric(stats::quantile(x, 0.025, na.rm = TRUE)) else NA_real_,
               q975 = if (nrow(sub) > 1) as.numeric(stats::quantile(x, 0.975, na.rm = TRUE)) else NA_real_,
               min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE), stringsAsFactors = FALSE)
  }))
}))
down_summary$eps <- eps_main
down_summary$threshold <- t_main
down_summary$n_pairs_sampled <- down_pairs
rownames(down_summary) <- NULL
utils::write.csv(down_summary, file.path(out_dir, "downsampling_tcga.csv"), row.names = FALSE)
toc("downsampling_tcga", t0)

# ---------------------------------------------------------------------------
# Part 5: figure
# ---------------------------------------------------------------------------
t0 <- tic()
obs_reversal <- ref_summary$frac_simpson_eps0p05_t0p2      # recomputed observed rate (paper: 13.3%)
obs_mixed_hv <- ref_summary$frac_mixed_eps0p05_given_maxabs_ge0p2   # paper: 99.5%
obs_mixed_gw <- 0.948                                       # paper: genome-wide, 23,170,038 pairs, eps = 0
x_labels <- c(as.character(balanced_ns), sprintf("TCGA\n(%d–%d)", min(N_c), max(N_c)))
x_pos <- seq_along(x_labels)
rho_cols <- c("#0072B2", "#E69F00", "#009E73")
min_elig_plot <- 50

get_curve <- function(shift, rho, metric) {
  bal <- null_rates[null_rates$mean_shift == shift & null_rates$rho == rho & null_rates$design == "balanced", ]
  bal <- bal[match(balanced_ns, bal$n_per_context), ]
  imb <- null_rates[null_rates$mean_shift == shift & null_rates$rho == rho & null_rates$design == "tcga_imbalanced", ]
  y <- c(bal[[metric]], imb[[metric]])
  if (metric == "frac_simpson_eps0p05_t0p2") {
    ne <- c(bal$n_elig_simpson_t0p2, imb$n_elig_simpson_t0p2)
    y[ne < min_elig_plot] <- NA
  }
  y
}

ref_lty <- c(reversal = 2, mixed_gw = 4, mixed_hv = 3)
panel <- function(shift, metric, ylab, main, label, ylim, ref_lines, ref_ltys) {
  plot(NA, xlim = c(0.7, length(x_pos) + 0.3), ylim = ylim, xaxt = "n", xlab = "Context size (samples per context)",
       ylab = ylab, main = "", las = 1)
  title(main = main, font.main = 2, cex.main = 1.05)
  mtext(label, side = 3, line = 0.9, adj = 0, at = par("usr")[1] - 0.13 * diff(par("usr")[1:2]), font = 2, cex = 1.25)
  axis(1, at = x_pos, labels = x_labels, padj = 0.4)
  abline(v = length(x_pos) - 0.5, col = "grey80", lty = 3)
  for (i in seq_along(ref_lines)) abline(h = ref_lines[i], col = "grey35", lty = ref_ltys[i], lwd = 1.6)
  for (i in seq_along(rhos)) {
    y <- get_curve(shift, rhos[i], metric)
    nb <- length(balanced_ns)
    lines(x_pos[1:nb], y[1:nb], col = rho_cols[i], lwd = 2.2)
    points(x_pos[1:nb], y[1:nb], col = rho_cols[i], pch = 16, cex = 1.1)
    points(x_pos[nb + 1], y[nb + 1], col = rho_cols[i], pch = 18, cex = 1.7)
  }
}

pdf_file <- file.path(out_dir, "Figure_NullCalibration.pdf")
grDevices::cairo_pdf(pdf_file, width = 10, height = 8.6)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.16))
par(mar = c(4.6, 4.8, 3.2, 1), oma = c(0, 0.6, 0, 0), mgp = c(3, 0.7, 0))
rev_ylab <- expression("Simpson reversal rate among |" * r[global] * "|" >= 0.2)
mix_ylab <- expression("Mixed-sign rate (" * epsilon == 0.05 * ")")
panel("none", "frac_simpson_eps0p05_t0p2", rev_ylab, "Reversal, No Between-Context Mean Shift", "A.", c(0, 0.6),
      obs_reversal, ref_lty["reversal"])
panel("tcga_matched", "frac_simpson_eps0p05_t0p2", rev_ylab, "Reversal, TCGA-Matched Mean Shifts", "B.", c(0, 0.6),
      obs_reversal, ref_lty["reversal"])
panel("none", "frac_mixed_eps0p05", mix_ylab, "Mixed Sign, No Between-Context Mean Shift", "C.", c(0, 1.05),
      c(obs_mixed_gw, obs_mixed_hv), ref_lty[c("mixed_gw", "mixed_hv")])
panel("tcga_matched", "frac_mixed_eps0p05", mix_ylab, "Mixed Sign, TCGA-Matched Mean Shifts", "D.", c(0, 1.05),
      c(obs_mixed_gw, obs_mixed_hv), ref_lty[c("mixed_gw", "mixed_hv")])
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", ncol = 4, bty = "n", cex = 0.95, x.intersp = 0.8, text.width = NA,
       legend = c(expression(rho == 0), expression(rho == 0.1), expression(rho == 0.3), "Balanced contexts",
                  "TCGA-like imbalanced sizes", sprintf("TCGA observed reversal (%.1f%%)", 100 * obs_reversal),
                  sprintf("TCGA genome-wide mixed sign (%.1f%%)", 100 * obs_mixed_gw),
                  sprintf("TCGA high-variance mixed sign (%.1f%%)", 100 * obs_mixed_hv)),
       col = c(rho_cols, "grey20", "grey20", "grey35", "grey35", "grey35"),
       lty = c(1, 1, 1, NA, NA, ref_lty), lwd = c(2.2, 2.2, 2.2, NA, NA, 1.6, 1.6, 1.6),
       pch = c(NA, NA, NA, 16, 18, NA, NA, NA), pt.cex = c(NA, NA, NA, 1.1, 1.7, NA, NA, NA))
dev.off()
toc("figure", t0)

# ---------------------------------------------------------------------------
# Runtime and session info
# ---------------------------------------------------------------------------
total_s <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
rt <- c(sprintf("null_simulation.R finished %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        sprintf("QUICK mode: %s", QUICK),
        sprintf("cores used (mclapply): %d", n_cores),
        sprintf("pairs per scenario: %d; scenarios: %d", pairs_per_scenario, nrow(null_rates)),
        sprintf("downsampling: %d pairs, %d replicates per n, n in {%s}", down_pairs, down_reps, paste(down_ns, collapse = ", ")),
        sprintf("total wall time: %.1f s (%.1f min)", total_s, total_s / 60),
        "section wall times (s):",
        sprintf("  %-22s %8.1f", names(section_times), unlist(section_times)),
        "", capture.output(print(sessionInfo())))
writeLines(rt, file.path(out_dir, "runtime.txt"))
message(sprintf("Done in %.1f min. Outputs in %s", total_s / 60, out_dir))
