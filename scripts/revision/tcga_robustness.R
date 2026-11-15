#!/usr/bin/env Rscript

# TCGA robustness analyses for the BMC Genomics revision.
#
# Run from the repo root:
#   Rscript scripts/revision/tcga_robustness.R
#
# Optional environment variables:
#   TCGA_ROBUSTNESS_CORES     forked workers (default: min(44, detectCores() - 4))
#   TCGA_ROBUSTNESS_B_SAMPLE  sample-level bootstrap replicates (default 200)
#   TCGA_ROBUSTNESS_B_PAIR    pair-level bootstrap replicates (default 1000)
#
# All outputs go to results/revision/tcga_robustness/. Nothing outside that directory is
# modified; existing scripts, caches, results and .tex files are only read. Each CSV is
# written as soon as its analysis finishes so partial results survive an interruption.
#
# Conventions reused from the paper's code (scripts/robustness_signflip.R,
# scripts/heterogeneity_I2_FDR.R, R/contextcorr.R):
#   * high-variance domain: top 5,000 mRNAs x top 200 miRNAs by global variance;
#   * a pair-cohort cell is valid (measurable) when r_c is finite and |r_c| < 1;
#   * sign tolerance eps = 0.05; majority vote over cohorts; ties and r_global = 0 excluded;
#   * cohorts need N >= 20 (all 31 qualify); pairs need >= 10 valid cohorts;
#   * Fisher z with Var = 1/(N_c - 3); Cochran Q with weights N_c - 3; I2 = max(0, (Q-df)/Q);
#   * BH FDR across the 10^6 pairs for the heterogeneity test.
# The genome-wide domain follows 300-correl.Rmd / 311-zetachi.Rmd: all 20,289 mRNAs x the
# 1,142 miRNAs that are not identically zero (23,170,038 pairs).

suppressPackageStartupMessages(library(parallel))

t_start <- Sys.time()
timing <- list()
tic <- function() Sys.time()
toc <- function(label, t0) {
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timing[[label]] <<- el
  message(sprintf("[%s] %.1f s", label, el))
  invisible(el)
}

source("00-paths.R")
source(file.path("R", "contextcorr.R"))

out_dir <- file.path("results", "revision", "tcga_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wcsv <- function(df, name) {
  utils::write.csv(df, file.path(out_dir, name), row.names = FALSE)
  message("  wrote ", name)
}

## ---------------------------------------------------------------------------
## Parameters
## ---------------------------------------------------------------------------
G <- 5000L; M <- 200L
eps <- 0.05
t_strong <- 0.2
min_k <- 10L
min_n_per_cohort <- 20L
t_switch <- 0.3
stab_thresh <- 0.9
z_crit <- 1.96
n_gw_sub <- 100000L

B_pair <- as.integer(Sys.getenv("TCGA_ROBUSTNESS_B_PAIR", "1000"))
B_samp <- as.integer(Sys.getenv("TCGA_ROBUSTNESS_B_SAMPLE", "200"))
seed_main <- 20260903L
seed_boot_base <- 7000000L
seed_pair_base <- 3000000L
n_cores_env <- Sys.getenv("TCGA_ROBUSTNESS_CORES", "")
n_cores <- if (nzchar(n_cores_env)) as.integer(n_cores_env) else max(1L, min(44L, detectCores() - 4L))
message(sprintf("Workers: %d; B_samp = %d; B_pair = %d", n_cores, B_samp, B_pair))

pmap <- function(X, FUN, cores = n_cores) {
  cores <- max(1L, min(cores, length(X)))
  res <- if (cores > 1L) mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE) else lapply(X, FUN)
  bad <- vapply(res, function(x) inherits(x, "try-error"), logical(1))
  if (any(bad)) stop("parallel worker failed: ", paste(unique(unlist(lapply(res[bad], as.character))), collapse = " | "))
  res
}

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------

# Row-standardize a features x samples matrix. Rows with zero variance (e.g. every
# sample at the log2(10) floor) become NaN, so "non-measurable" never depends on
# floating-point round-off.
prep_rows <- function(X) {
  mu <- rowMeans(X)
  Xc <- X - mu
  ss <- sqrt(rowSums(Xc^2))
  zero <- rowSums(X != X[, 1L]) == 0L
  ss[zero] <- 0
  Xs <- Xc / ss
  if (any(zero)) Xs[zero, ] <- NaN
  list(Xs = Xs, mu = mu, ss = ss, zero = zero)
}
cor_rows <- function(X, Y) tcrossprod(prep_rows(X)$Xs, prep_rows(Y)$Xs)

# Pooled (global) r rebuilt from within-cohort sufficient statistics:
#   Sxy = sum_c sum_{i in c}(x_i - xbar_c)(y_i - ybar_c) + sum_c n_c (xbar_c - xbar)(ybar_c - ybar)
pooled_r_from_parts <- function(Sxy_within, MUX, MUY, SSX, SSY, n) {
  N <- sum(n)
  xbar <- as.vector(MUX %*% n) / N
  ybar <- as.vector(MUY %*% n) / N
  DX <- MUX - xbar; DY <- MUY - ybar
  Sxx <- rowSums(SSX^2) + as.vector(DX^2 %*% n)
  Syy <- rowSums(SSY^2) + as.vector(DY^2 %*% n)
  Sxy <- Sxy_within + as.vector(sweep(DX, 2, n, "*") %*% t(DY))
  Sxy / as.vector(outer(sqrt(Sxx), sqrt(Syy)))
}

row_max_na <- function(Mat) { o <- rep(-Inf, nrow(Mat)); for (j in seq_len(ncol(Mat))) o <- pmax(o, Mat[, j], na.rm = TRUE); o[!is.finite(o)] <- NA; o }
row_min_na <- function(Mat) { o <- rep(Inf, nrow(Mat)); for (j in seq_len(ncol(Mat))) o <- pmin(o, Mat[, j], na.rm = TRUE); o[!is.finite(o)] <- NA; o }

cohen_kappa <- function(a, b) {
  lv <- sort(union(unique(a), unique(b)))
  tab <- table(factor(a, levels = lv), factor(b, levels = lv))
  n <- sum(tab); po <- sum(diag(tab)) / n; pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  if (pe >= 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")
fmt_int <- function(x) format(round(as.numeric(x)), big.mark = ",", scientific = FALSE, trim = TRUE)
chunk_seq <- function(n, nch) split(seq_len(n), rep_len(seq_len(min(nch, n)), n))

## ---------------------------------------------------------------------------
## Load data
## ---------------------------------------------------------------------------
t0 <- tic()
f1 <- file.path(paths$mirtcga, "matchedNoGBM.Rda"); f2 <- file.path(paths$mirtcga, "matchedNoGBM.rda")
if (file.exists(f1)) load(f1) else if (file.exists(f2)) load(f2) else stop("Cannot find matchedNoGBM in ", paths$mirtcga)
matched_file <- if (file.exists(f1)) normalizePath(f1) else normalizePath(f2)
stopifnot(identical(colnames(miRdata), colnames(mRNAdata)), length(cancerType) == ncol(miRdata))
if (!exists("sampleType")) sampleType <- factor(rep("unknown", ncol(miRdata)))
cohorts <- levels(cancerType); C <- length(cohorts)
idx_by_ct <- lapply(cohorts, function(ct) which(cancerType == ct)); names(idx_by_ct) <- cohorts
n_by_ct <- vapply(idx_by_ct, length, integer(1)); N_total <- ncol(miRdata)
w_ct <- n_by_ct - 3
toc("load_data", t0)

## ---------------------------------------------------------------------------
## (a) Per-cohort sample sizes  ->  cohort_sizes.csv, cohort_sizes_summary.csv
## ---------------------------------------------------------------------------
t0 <- tic()
st <- table(cancerType, sampleType)
gc_col <- function(nm) if (nm %in% colnames(st)) as.integer(st[, nm]) else rep(0L, nrow(st))
cohort_sizes <- data.frame(row_id = cohorts, cohort = cohorts, N = as.integer(n_by_ct),
                           n_tumor = gc_col("tumor"), n_normal = gc_col("normal"), n_metastasis = gc_col("metastasis"),
                           weight_N_minus_3 = as.integer(n_by_ct) - 3L, frac_of_total = as.numeric(n_by_ct) / N_total,
                           stringsAsFactors = FALSE)
cohort_sizes <- cohort_sizes[order(cohort_sizes$N, cohort_sizes$cohort), ]
cohort_sizes$rank_by_N <- seq_len(nrow(cohort_sizes)); rownames(cohort_sizes) <- NULL
wcsv(cohort_sizes, "cohort_sizes.csv")
lt50 <- cohort_sizes$cohort[cohort_sizes$N < 50]; lt100 <- cohort_sizes$cohort[cohort_sizes$N < 100]
cs_summary <- data.frame(
  row_id = sprintf("S%02d", 1:12),
  quantity = c("n_cohorts", "N_total", "min_N", "min_N_cohort", "max_N", "max_N_cohort", "median_N", "mean_N",
               "n_cohorts_N_lt_50", "cohorts_N_lt_50", "n_cohorts_N_lt_100", "cohorts_N_lt_100"),
  value = c(C, N_total, min(cohort_sizes$N), cohort_sizes$cohort[which.min(cohort_sizes$N)], max(cohort_sizes$N),
            cohort_sizes$cohort[which.max(cohort_sizes$N)], median(cohort_sizes$N), round(mean(cohort_sizes$N), 2),
            length(lt50), paste(lt50, collapse = ";"), length(lt100), paste(lt100, collapse = ";")),
  stringsAsFactors = FALSE)
st_tot <- data.frame(row_id = sprintf("S%02d", 13:(12 + ncol(st))), quantity = paste0("N_total_", colnames(st)),
                     value = as.character(colSums(st)), stringsAsFactors = FALSE)
cs_summary <- rbind(cs_summary, st_tot)
wcsv(cs_summary, "cohort_sizes_summary.csv")
min_cohort <- cs_summary$value[4]
toc("a_cohort_sizes", t0)

## ---------------------------------------------------------------------------
## Per-cohort standardized matrices (shared with forked workers)
## ---------------------------------------------------------------------------
t0 <- tic()
mir_all_zero <- apply(miRdata, 1, function(x) all(x == 0))
Y_gw <- miRdata[!mir_all_zero, , drop = FALSE]
n_genes_gw <- nrow(mRNAdata); n_mirs_gw <- nrow(Y_gw); P_gw <- as.numeric(n_genes_gw) * n_mirs_gw
XS <- vector("list", C); YS <- vector("list", C)
zv_genes_by_ct <- integer(C); zv_mirs_by_ct <- integer(C)
for (j in seq_len(C)) {
  px <- prep_rows(mRNAdata[, idx_by_ct[[j]], drop = FALSE])
  py <- prep_rows(Y_gw[, idx_by_ct[[j]], drop = FALSE])
  XS[[j]] <- px$Xs; YS[[j]] <- py$Xs
  zv_genes_by_ct[j] <- sum(px$zero); zv_mirs_by_ct[j] <- sum(py$zero)
  rm(px, py)
}
gc(verbose = FALSE)
toc("standardize_by_cohort", t0)

## ---------------------------------------------------------------------------
## High-variance domain, full data
## ---------------------------------------------------------------------------
t0 <- tic()
gene_var <- apply(mRNAdata, 1, var); mir_var <- apply(miRdata, 1, var)
genes <- names(sort(gene_var, decreasing = TRUE))[1:G]
mirs <- names(sort(mir_var, decreasing = TRUE))[1:M]
X_hv <- mRNAdata[genes, , drop = FALSE]; Y_hv <- miRdata[mirs, , drop = FALSE]
P_hv <- G * M
pair_gene <- rep(genes, times = M); pair_mir <- rep(mirs, each = G)   # as.vector order of G x M

R <- matrix(NA_real_, P_hv, C, dimnames = list(NULL, cohorts))
MUX <- matrix(0, G, C); MUY <- matrix(0, M, C); SSX <- matrix(0, G, C); SSY <- matrix(0, M, C)
ZVX_hv <- matrix(FALSE, G, C); ZVY_hv <- matrix(FALSE, M, C)
Sxy_within <- numeric(P_hv)
for (j in seq_len(C)) {
  idx <- idx_by_ct[[j]]
  px <- prep_rows(X_hv[, idx, drop = FALSE]); py <- prep_rows(Y_hv[, idx, drop = FALSE])
  r <- as.vector(tcrossprod(px$Xs, py$Xs))
  R[, j] <- r
  W <- r * as.vector(outer(px$ss, py$ss)); W[!is.finite(W)] <- 0
  Sxy_within <- Sxy_within + W
  MUX[, j] <- px$mu; MUY[, j] <- py$mu; SSX[, j] <- px$ss; SSY[, j] <- py$ss
  ZVX_hv[, j] <- px$zero; ZVY_hv[, j] <- py$zero
}
rg <- as.vector(cor_rows(X_hv, Y_hv))
max_diff_pooled <- max(abs(rg - pooled_r_from_parts(Sxy_within, MUX, MUY, SSX, SSY, n_by_ct)))
sg <- sign(rg); rm(Sxy_within)

valid <- is.finite(R) & abs(R) < 1
valid[, n_by_ct < min_n_per_cohort] <- FALSE
n_cells_hv <- length(R)
n_nonmeasurable_hv <- sum(!is.finite(R))
n_degenerate_hv <- sum(is.finite(R) & abs(R) >= 1)
k_valid <- rowSums(valid)
S <- matrix(0L, P_hv, C)
S[valid & R >= eps] <- 1L; S[valid & R <= -eps] <- -1L
n_pos <- rowSums(S == 1L); n_neg <- rowSums(S == -1L)
maj <- sign(n_pos - n_neg); mixed <- n_pos > 0L & n_neg > 0L
Ra <- abs(R); Ra[!valid] <- NA
max_abs <- row_max_na(Ra); max_abs[is.na(max_abs)] <- 0; rm(Ra)

elig_simpson <- k_valid >= min_k & abs(rg) >= t_strong & maj != 0 & sg != 0
reversal <- elig_simpson & sg != maj; concordant <- elig_simpson & sg == maj
n_elig_simpson <- sum(elig_simpson); n_reversal <- sum(reversal)
rate_simpson <- n_reversal / n_elig_simpson
elig_simpson_t0 <- k_valid >= min_k & maj != 0 & sg != 0
rate_simpson_t0 <- sum(elig_simpson_t0 & sg != maj) / sum(elig_simpson_t0)
elig_mixed <- k_valid >= min_k & max_abs >= t_strong
n_elig_mixed <- sum(elig_mixed)
rate_mixed_hv <- mean(mixed[elig_mixed])

Z <- contextcorr_fisher_z(R); Z[!valid] <- NA
sum_w <- numeric(P_hv); sum_w2 <- numeric(P_hv); sum_wz <- numeric(P_hv); sum_wz2 <- numeric(P_hv)
for (j in seq_len(C)) {
  v <- valid[, j]; z <- Z[, j]; z[!v] <- 0
  sum_w <- sum_w + w_ct[j] * v; sum_w2 <- sum_w2 + w_ct[j]^2 * v
  sum_wz <- sum_wz + w_ct[j] * z; sum_wz2 <- sum_wz2 + w_ct[j] * z^2
}
df_het <- k_valid - 1L
Q <- sum_wz2 - sum_wz^2 / pmax(sum_w, 1e-12)
ok_het <- df_het >= 1L & is.finite(Q) & Q >= 0
I2 <- rep(NA_real_, P_hv); I2[ok_het] <- pmax(0, (Q[ok_het] - df_het[ok_het]) / pmax(Q[ok_het], 1e-12))
p_het <- rep(NA_real_, P_hv); p_het[ok_het] <- pchisq(Q[ok_het], df = df_het[ok_het], lower.tail = FALSE)
fdr_het <- rep(NA_real_, P_hv); fdr_het[ok_het] <- p.adjust(p_het[ok_het], method = "BH")
r_fixed <- contextcorr_inv_fisher_z(sum_wz / pmax(sum_w, 1e-12))
median_I2 <- median(I2[ok_het]); I2_q25 <- unname(quantile(I2[ok_het], .25)); I2_q75 <- unname(quantile(I2[ok_het], .75))
frac_fdr <- mean(fdr_het[ok_het] < 0.05)

set.seed(seed_main)
chk <- sample.int(P_hv, 300)
max_diff_I2_helper <- max(abs(vapply(chk, function(p) heterogeneity_Q_I2(R[p, ], n_by_ct)$I2, numeric(1)) - I2[chk]))

sf <- read.csv(file.path("results", "robustness_signflip_summary.csv"))
sf_row <- sf[sf$eps == eps & sf$threshold == t_strong, ]
het_summary <- read.csv(file.path("results", "heterogeneity_summary.csv"))
het_row <- het_summary[het_summary$group == "all (k>=10)", ]
headline <- data.frame(
  row_id = sprintf("H%02d", 1:8),
  quantity = c("simpson_reversal_rate_rg_ge_0.2_eps_0.05", "n_eligible_simpson_rg_ge_0.2",
               "mixed_sign_rate_maxabs_ge_0.2_eps_0.05", "n_eligible_mixed_maxabs_ge_0.2",
               "median_I2", "I2_q25", "I2_q75", "frac_pairs_Q_fdr_lt_0.05"),
  paper_value = c(sf_row$frac_simpson, sf_row$eligible_simpson, sf_row$frac_mixed, sf_row$eligible_mixed,
                  het_row$median_I2, het_row$I2_q25, het_row$I2_q75, het_row$frac_fdr_lt_0p05),
  recomputed_value = c(rate_simpson, n_elig_simpson, rate_mixed_hv, n_elig_mixed, median_I2, I2_q25, I2_q75, frac_fdr),
  paper_source = c(rep("results/robustness_signflip_summary.csv (eps=0.05, threshold=0.2)", 4),
                   rep("results/heterogeneity_summary.csv (group 'all (k>=10)')", 4)),
  stringsAsFactors = FALSE)
headline$abs_diff <- abs(headline$paper_value - headline$recomputed_value)
wcsv(headline, "headline_reproduction.csv")
toc("hv_full_data", t0)

## ---------------------------------------------------------------------------
## Genome-wide pass (parallel over gene chunks)
## ---------------------------------------------------------------------------
t0 <- tic()
gchunks <- chunk_seq(n_genes_gw, n_cores * 3L)
gw_res <- pmap(gchunks, function(rows) {
  nr <- length(rows); np <- nr * n_mirs_gw
  k <- integer(np); nanc <- integer(np); deg <- integer(np)
  p0 <- integer(np); n0 <- integer(np); p5 <- integer(np); n5 <- integer(np)
  for (j in seq_len(C)) {
    r <- as.vector(tcrossprod(XS[[j]][rows, , drop = FALSE], YS[[j]]))
    fin <- is.finite(r); nanc <- nanc + !fin
    v <- fin & abs(r) < 1; deg <- deg + (fin & !v)
    if (n_by_ct[j] >= min_n_per_cohort) {
      k <- k + v
      p0 <- p0 + (v & r > 0); n0 <- n0 + (v & r < 0)
      p5 <- p5 + (v & r >= eps); n5 <- n5 + (v & r <= -eps)
    }
  }
  c(n_cells_nan = sum(nanc), n_cells_deg = sum(deg),
    n_k0 = sum(k == 0L), n_k1 = sum(k == 1L), n_k2 = sum(k >= 2L), n_k10 = sum(k >= min_k),
    n_mixed0 = sum(p0 > 0L & n0 > 0L), n_mixed5 = sum(p5 > 0L & n5 > 0L),
    n_mixed0_k2 = sum((p0 > 0L & n0 > 0L) & k >= 2L), n_mixed5_k2 = sum((p5 > 0L & n5 > 0L) & k >= 2L),
    n_mixed5_k10 = sum((p5 > 0L & n5 > 0L) & k >= min_k), n_pairs = np)
})
gw <- colSums(do.call(rbind, gw_res))
stopifnot(gw[["n_pairs"]] == P_gw)
n_cells_gw <- P_gw * C
n_nonmeasurable_gw <- gw[["n_cells_nan"]]; n_degenerate_gw <- gw[["n_cells_deg"]]
n_pairs_k0_gw <- gw[["n_k0"]]; n_pairs_k1_gw <- gw[["n_k1"]]; n_pairs_k2_gw <- gw[["n_k2"]]; n_pairs_k10_gw <- gw[["n_k10"]]
gw_rate_eps0_all <- gw[["n_mixed0"]] / P_gw
gw_rate_eps0_k2 <- gw[["n_mixed0_k2"]] / n_pairs_k2_gw
gw_rate_eps05_all <- gw[["n_mixed5"]] / P_gw
gw_rate_eps05_k2 <- gw[["n_mixed5_k2"]] / n_pairs_k2_gw
gw_rate_eps05_k10 <- gw[["n_mixed5_k10"]] / n_pairs_k10_gw
rm(gw_res); gc(verbose = FALSE)
toc("genome_wide_pass", t0)

## ---------------------------------------------------------------------------
## (f) Measurability  ->  measurability.csv, measurability_by_cohort.csv
## ---------------------------------------------------------------------------
t0 <- tic()
measurability <- data.frame(
  row_id = sprintf("M%02d", 1:27),
  domain = c(rep("high-variance (5,000 x 200)", 12), rep("genome-wide (20,289 x 1,142)", 15)),
  quantity = c("n_pairs", "n_pair_cohort_cells", "n_nonmeasurable_cells_zero_variance", "frac_nonmeasurable_cells",
               "n_degenerate_cells_abs_r_ge_1", "min_k_valid", "n_pairs_k_valid_eq_31", "n_pairs_k_valid_ge_10",
               "denominator_simpson_rg_ge_0.2", "denominator_mixed_maxabs_ge_0.2", "denominator_I2_and_Q_FDR", "n_pairs_k_valid_lt_31",
               "n_pairs", "n_miRNAs_excluded_all_zero", "n_pair_cohort_cells", "n_nonmeasurable_cells_zero_variance",
               "frac_nonmeasurable_cells", "n_degenerate_cells_abs_r_ge_1", "n_pairs_k_valid_eq_0", "n_pairs_k_valid_eq_1",
               "n_pairs_k_valid_ge_2", "n_pairs_k_valid_ge_10", "mixed_sign_rate_eps0_denominator_all_pairs",
               "mixed_sign_rate_eps0_denominator_k_ge_2", "mixed_sign_rate_eps0.05_denominator_all_pairs",
               "mixed_sign_rate_eps0.05_denominator_k_ge_2", "mixed_sign_rate_eps0.05_denominator_k_ge_10"),
  value = c(P_hv, n_cells_hv, n_nonmeasurable_hv, n_nonmeasurable_hv / n_cells_hv, n_degenerate_hv, min(k_valid),
            sum(k_valid == C), sum(k_valid >= min_k), n_elig_simpson, n_elig_mixed, sum(ok_het), sum(k_valid < C),
            P_gw, sum(mir_all_zero), n_cells_gw, n_nonmeasurable_gw, n_nonmeasurable_gw / n_cells_gw, n_degenerate_gw,
            n_pairs_k0_gw, n_pairs_k1_gw, n_pairs_k2_gw, n_pairs_k10_gw,
            gw_rate_eps0_all, gw_rate_eps0_k2, gw_rate_eps05_all, gw_rate_eps05_k2, gw_rate_eps05_k10),
  note = c("top 5,000 mRNAs x top 200 miRNAs by global variance", "pairs x 31 cohorts",
           "r_c is NaN because the mRNA or the miRNA has zero variance in that cohort", "",
           "finite r_c with |r_c| >= 1 (treated as non-valid, as in scripts/robustness_signflip.R)",
           "smallest number of valid cohorts over the 10^6 pairs", "", "all pairs passing min_k = 10",
           "k >= 10 & |r_global| >= 0.2 & majority != 0 & r_global != 0 (reproduces 13.3%)",
           "k >= 10 & max_c |r_c| >= 0.2 (reproduces 99.5%)",
           "pairs with k >= 2 (reproduces median I2 and the FDR fraction)", "pairs with at least one non-valid cohort",
           "all 20,289 mRNAs x 1,142 miRNAs not identically zero (300-correl.Rmd)", "of 1,146 miRNAs in matchedNoGBM",
           "pairs x 31 cohorts", "", "", "", "no measurable cohort", "exactly one measurable cohort (spread undefined)", "", "",
           "reproduces the paper's 94.8%: eps = 0, denominator ALL 23,170,038 pairs including those with < 2 measurable cohorts (311-zetachi.Rmd, S/TOT)",
           "same numerator, denominator restricted to pairs with >= 2 measurable cohorts",
           "eps = 0.05 (the Methods default tolerance), denominator all pairs",
           "eps = 0.05, denominator pairs with >= 2 measurable cohorts",
           "eps = 0.05, denominator pairs with >= 10 measurable cohorts"),
  stringsAsFactors = FALSE)
wcsv(measurability, "measurability.csv")
meas_by_ct <- data.frame(row_id = cohorts, cohort = cohorts, N = as.integer(n_by_ct),
                         n_zero_variance_mRNAs_of_20289 = zv_genes_by_ct, n_zero_variance_miRNAs_of_1142 = zv_mirs_by_ct,
                         n_zero_variance_hv_mRNAs_of_5000 = colSums(ZVX_hv), n_zero_variance_hv_miRNAs_of_200 = colSums(ZVY_hv),
                         n_nonmeasurable_hv_cells = colSums(!is.finite(R)),
                         n_degenerate_hv_cells = colSums(is.finite(R) & abs(R) >= 1),
                         n_measurable_hv_pairs = colSums(valid), stringsAsFactors = FALSE)
wcsv(meas_by_ct, "measurability_by_cohort.csv")
toc("f_measurability", t0)

## ---------------------------------------------------------------------------
## Per-cohort p-values, CI gate, per-cohort BH  (inputs to (c), (d), (e))
## ---------------------------------------------------------------------------
t0 <- tic()
Pm <- matrix(NA_real_, P_hv, C)
sig_ci <- matrix(FALSE, P_hv, C)
for (j in seq_len(C)) {
  v <- valid[, j]; zj <- Z[v, j]
  Pm[v, j] <- 2 * pnorm(-abs(zj) * sqrt(w_ct[j]))
  sig_ci[v, j] <- abs(zj) > z_crit / sqrt(w_ct[j])
}
sig_p <- valid & !is.na(Pm) & Pm < 0.05
qres <- pmap(as.list(seq_len(C)), function(j) {
  v <- valid[, j]; q <- rep(NA_real_, P_hv)
  q[v] <- p.adjust(Pm[v, j], method = "BH")
  s <- v & !is.na(q) & q < 0.05
  list(sig = s, min_abs_r = if (any(s)) min(abs(R[s, j])) else NA_real_)
})
sig_q <- matrix(FALSE, P_hv, C)
min_abs_r_q <- numeric(C)
for (j in seq_len(C)) { sig_q[, j] <- qres[[j]]$sig; min_abs_r_q[j] <- qres[[j]]$min_abs_r }
rm(qres); gc(verbose = FALSE)
toc("per_cohort_pvalues", t0)

## ---------------------------------------------------------------------------
## (d) Context-switching pairs  ->  context_switching.csv, context_switching_top20.csv
## ---------------------------------------------------------------------------
t0 <- tic()
Rv <- R; Rv[!valid] <- NA
r_max <- row_max_na(Rv); r_min <- row_min_na(Rv)
cand <- !is.na(r_max) & !is.na(r_min) & r_max >= t_switch & r_min <= -t_switch
pos_nom <- rowSums(valid & R >= t_switch & sig_p) > 0
neg_nom <- rowSums(valid & R <= -t_switch & sig_p) > 0
switch_nom <- cand & pos_nom & neg_nom
pos_q <- rowSums(valid & R >= t_switch & sig_q) > 0
neg_q <- rowSums(valid & R <= -t_switch & sig_q) > 0
switch_q <- cand & pos_q & neg_q

cand_idx <- which(cand)
bh_res <- pmap(chunk_seq(length(cand_idx), n_cores * 2L), function(ii) {
  pid <- cand_idx[ii]; m <- length(pid)
  hp <- logical(m); hn <- logical(m)
  jp <- integer(m); jn <- integer(m); qp <- rep(NA_real_, m); qn <- rep(NA_real_, m)
  for (a in seq_len(m)) {
    p <- pid[a]; v <- which(valid[p, ])
    q <- p.adjust(Pm[p, v], method = "BH")
    rr <- R[p, v]
    okp <- q < 0.05 & rr >= t_switch; okn <- q < 0.05 & rr <= -t_switch
    if (any(okp)) { b <- which(okp)[which.max(rr[okp])]; hp[a] <- TRUE; jp[a] <- v[b]; qp[a] <- q[b] }
    if (any(okn)) { b <- which(okn)[which.min(rr[okn])]; hn[a] <- TRUE; jn[a] <- v[b]; qn[a] <- q[b] }
  }
  list(i = ii, hp = hp, hn = hn, jp = jp, jn = jn, qp = qp, qn = qn)
})
hp_all <- logical(length(cand_idx)); hn_all <- logical(length(cand_idx))
jp_all <- integer(length(cand_idx)); jn_all <- integer(length(cand_idx))
qp_all <- rep(NA_real_, length(cand_idx)); qn_all <- rep(NA_real_, length(cand_idx))
for (b in bh_res) { hp_all[b$i] <- b$hp; hn_all[b$i] <- b$hn; jp_all[b$i] <- b$jp; jn_all[b$i] <- b$jn; qp_all[b$i] <- b$qp; qn_all[b$i] <- b$qn }
rm(bh_res)
switch_bh <- logical(P_hv); switch_bh[cand_idx[hp_all & hn_all]] <- TRUE

context_switching <- data.frame(
  row_id = sprintf("C%02d", 1:8),
  criterion = c("r_max >= 0.3 and r_min <= -0.3 (no significance requirement)",
                "both extremes nominally significant (Fisher-z p < 0.05)",
                "both extremes significant after BH within pair (q < 0.05 across that pair's measurable cohorts)",
                "both extremes significant after per-cohort BH across pairs (q < 0.05)",
                "nominal criterion, among baseline Simpson reversals (|r_global| >= 0.2)",
                "within-pair-BH criterion, among baseline Simpson reversals (|r_global| >= 0.2)",
                "nominal criterion, among pairs eligible for the reversal rate (|r_global| >= 0.2)",
                "within-pair-BH criterion, among pairs eligible for the reversal rate (|r_global| >= 0.2)"),
  n_pairs = c(sum(cand), sum(switch_nom), sum(switch_bh), sum(switch_q), sum(switch_nom & reversal),
              sum(switch_bh & reversal), sum(switch_nom & elig_simpson), sum(switch_bh & elig_simpson)),
  denominator = c(rep(P_hv, 4), n_reversal, n_reversal, n_elig_simpson, n_elig_simpson), stringsAsFactors = FALSE)
context_switching$fraction <- context_switching$n_pairs / context_switching$denominator
wcsv(context_switching, "context_switching.csv")

sw <- which(hp_all & hn_all)
if (length(sw) > 0) {
  pidx <- cand_idx[sw]
  top <- data.frame(pair_index = pidx, mRNA = pair_gene[pidx], miRNA = pair_mir[pidx], r_global = rg[pidx],
                    k_valid = k_valid[pidx], I2 = I2[pidx],
                    cohort_pos = cohorts[jp_all[sw]], r_pos = R[cbind(pidx, jp_all[sw])], N_pos = as.integer(n_by_ct[jp_all[sw]]),
                    p_pos = Pm[cbind(pidx, jp_all[sw])], q_pos_within_pair = qp_all[sw],
                    cohort_neg = cohorts[jn_all[sw]], r_neg = R[cbind(pidx, jn_all[sw])], N_neg = as.integer(n_by_ct[jn_all[sw]]),
                    p_neg = Pm[cbind(pidx, jn_all[sw])], q_neg_within_pair = qn_all[sw],
                    n_cohorts_pos_ge_0.3 = rowSums(Rv[pidx, , drop = FALSE] >= t_switch, na.rm = TRUE),
                    n_cohorts_neg_le_m0.3 = rowSums(Rv[pidx, , drop = FALSE] <= -t_switch, na.rm = TRUE),
                    baseline_simpson_reversal = reversal[pidx], stringsAsFactors = FALSE)
  top$spread <- top$r_pos - top$r_neg
  top <- top[order(-top$spread), ]
  top20 <- head(top, 20); top20$rank <- seq_len(nrow(top20))
  top20 <- top20[, c("rank", setdiff(names(top20), "rank"))]
} else top20 <- data.frame(rank = integer(0))
wcsv(top20, "context_switching_top20.csv")
rm(Rv); gc(verbose = FALSE)
toc("d_context_switching", t0)

## ---------------------------------------------------------------------------
## Rule ingredients (ii) (iii) (iv) and the p/q gates
## ---------------------------------------------------------------------------
t0 <- tic()
vote_w <- numeric(P_hv); vote_ci <- numeric(P_hv); vote_p <- numeric(P_hv); vote_q <- numeric(P_hv)
for (j in seq_len(C)) {
  vote_w <- vote_w + w_ct[j] * S[, j]
  sj <- sign(R[, j]); sj[!valid[, j]] <- 0        # sign(NaN) would poison the running sum
  vote_ci <- vote_ci + sj * sig_ci[, j]
  vote_p <- vote_p + S[, j] * sig_p[, j]
  vote_q <- vote_q + S[, j] * sig_q[, j]
}
stopifnot(all(is.finite(vote_w)), all(is.finite(vote_ci)), all(is.finite(vote_p)), all(is.finite(vote_q)))
maj_w <- sign(vote_w); maj_ci <- sign(vote_ci); maj_p <- sign(vote_p); maj_q <- sign(vote_q)

denom_dl <- sum_w - sum_w2 / pmax(sum_w, 1e-12)
tau2 <- pmax(0, (Q - df_het) / pmax(denom_dl, 1e-12))
s_wre <- numeric(P_hv); s_wre_z <- numeric(P_hv)
for (j in seq_len(C)) {
  v <- valid[, j]
  wre <- 1 / (1 / w_ct[j] + tau2[v])
  s_wre[v] <- s_wre[v] + wre; s_wre_z[v] <- s_wre_z[v] + wre * Z[v, j]
}
r_RE <- contextcorr_inv_fisher_z(s_wre_z / s_wre); r_RE[k_valid < 2L] <- NA
max_diff_rRE_helper <- max(abs(vapply(chk, function(p) random_effects_pool(R[p, ], n_by_ct)$r_RE, numeric(1)) - r_RE[chk]))
dir_RE <- sign(r_RE); dir_RE[is.na(dir_RE)] <- 0
toc("rule_ingredients", t0)

## ---------------------------------------------------------------------------
## (e) FDR-gated sensitivity  ->  fdr_gated.csv, fdr_gated_by_cohort.csv
## ---------------------------------------------------------------------------
t0 <- tic()
rev_idx <- which(reversal)
mq_rev <- maj_q[rev_idx]
n_rev_det_q <- sum(mq_rev != 0); n_rev_hold_q <- sum(mq_rev == maj[rev_idx])
n_rev_flip_q <- sum(mq_rev != 0 & mq_rev == sg[rev_idx]); n_rev_und_q <- sum(mq_rev == 0)
n_rev_novote_q <- sum(rowSums(sig_q[rev_idx, , drop = FALSE]) == 0)
mp_rev <- maj_p[rev_idx]; n_rev_hold_p <- sum(mp_rev == maj[rev_idx])
elig_q <- k_valid >= min_k & abs(rg) >= t_strong & maj_q != 0 & sg != 0
rate_q <- sum(elig_q & sg != maj_q) / sum(elig_q)
elig_q0 <- k_valid >= min_k & maj_q != 0 & sg != 0
rate_q0 <- sum(elig_q0 & sg != maj_q) / sum(elig_q0)
elig_p <- k_valid >= min_k & abs(rg) >= t_strong & maj_p != 0 & sg != 0
rate_p <- sum(elig_p & sg != maj_p) / sum(elig_p)
elig_p0 <- k_valid >= min_k & maj_p != 0 & sg != 0
rate_p0 <- sum(elig_p0 & sg != maj_p) / sum(elig_p0)
n_sig_q_cells <- sum(sig_q); n_valid_cells <- sum(valid)
per_cohort_q <- data.frame(row_id = cohorts, cohort = cohorts, N = as.integer(n_by_ct),
                           n_measurable_pairs = colSums(valid), n_q_lt_0.05 = colSums(sig_q),
                           frac_q_lt_0.05 = colSums(sig_q) / colSums(valid), n_p_lt_0.05 = colSums(sig_p),
                           frac_p_lt_0.05 = colSums(sig_p) / colSums(valid),
                           min_abs_r_with_q_lt_0.05 = min_abs_r_q, stringsAsFactors = FALSE)
wcsv(per_cohort_q, "fdr_gated_by_cohort.csv")
fdr_gated <- data.frame(
  row_id = sprintf("F%02d", 1:16),
  quantity = c("n_baseline_reversals", "n_baseline_reversals_with_determinate_q_gated_majority",
               "n_baseline_reversals_majority_holds_q_gated", "frac_majority_holds_among_determinate_q_gated",
               "frac_baseline_reversals_surviving_q_gated", "n_baseline_reversals_flipped_to_global_sign_q_gated",
               "n_baseline_reversals_tie_or_no_vote_q_gated", "n_baseline_reversals_with_zero_q_significant_cohorts",
               "overall_reversal_rate_q_gated_rg_ge_0.2", "n_eligible_q_gated_rg_ge_0.2",
               "overall_reversal_rate_q_gated_all_rg", "overall_reversal_rate_nominal_p_gated_rg_ge_0.2",
               "n_eligible_nominal_p_gated_rg_ge_0.2", "overall_reversal_rate_nominal_p_gated_all_rg",
               "frac_baseline_reversals_surviving_nominal_p_gated", "frac_measurable_cells_with_q_lt_0.05"),
  value = c(n_reversal, n_rev_det_q, n_rev_hold_q, n_rev_hold_q / n_rev_det_q, n_rev_hold_q / n_reversal,
            n_rev_flip_q, n_rev_und_q, n_rev_novote_q, rate_q, sum(elig_q), rate_q0, rate_p, sum(elig_p),
            rate_p0, n_rev_hold_p / n_reversal, n_sig_q_cells / n_valid_cells),
  numerator = c(NA, NA, NA, n_rev_hold_q, n_rev_hold_q, NA, NA, NA, sum(elig_q & sg != maj_q), NA,
                sum(elig_q0 & sg != maj_q), sum(elig_p & sg != maj_p), NA, sum(elig_p0 & sg != maj_p), n_rev_hold_p, n_sig_q_cells),
  denominator = c(NA, NA, NA, n_rev_det_q, n_reversal, NA, NA, NA, sum(elig_q), NA, sum(elig_q0), sum(elig_p),
                  NA, sum(elig_p0), n_reversal, n_valid_cells), stringsAsFactors = FALSE)
wcsv(fdr_gated, "fdr_gated.csv")
toc("e_fdr_gated", t0)

## ---------------------------------------------------------------------------
## Fixed random 100,000-pair subsample of the genome-wide domain
## ---------------------------------------------------------------------------
set.seed(seed_main)
gw_sub_p <- sort(sample.int(P_gw, n_gw_sub))
gw_sub_gi <- as.integer((gw_sub_p - 1) %% n_genes_gw + 1)
gw_sub_mj <- as.integer((gw_sub_p - 1) %/% n_genes_gw + 1)
ugenes <- sort(unique(gw_sub_gi)); gi_map <- match(gw_sub_gi, ugenes)
gw_sub_rates <- function(get_idx) {
  k <- integer(n_gw_sub); p0 <- integer(n_gw_sub); n0 <- integer(n_gw_sub)
  p5 <- integer(n_gw_sub); n5 <- integer(n_gw_sub)
  for (j in seq_len(C)) {
    idx <- get_idx(j)
    if (length(idx) < min_n_per_cohort) next
    Xs <- prep_rows(mRNAdata[ugenes, idx, drop = FALSE])$Xs
    Ys <- prep_rows(Y_gw[, idx, drop = FALSE])$Xs
    r <- numeric(n_gw_sub)
    for (s in seq(1L, n_gw_sub, by = 10000L)) {
      e <- min(s + 9999L, n_gw_sub); ii <- s:e
      r[ii] <- rowSums(Xs[gi_map[ii], , drop = FALSE] * Ys[gw_sub_mj[ii], , drop = FALSE])
    }
    rm(Xs, Ys)
    v <- is.finite(r) & abs(r) < 1
    k <- k + v; p0 <- p0 + (v & r > 0); n0 <- n0 + (v & r < 0)
    p5 <- p5 + (v & r >= eps); n5 <- n5 + (v & r <= -eps)
  }
  m0 <- p0 > 0L & n0 > 0L; m5 <- p5 > 0L & n5 > 0L
  c(gw_eps0_all = mean(m0), gw_eps0_k2 = mean(m0[k >= 2L]), gw_eps05_all = mean(m5),
    gw_eps05_k2 = mean(m5[k >= 2L]), gw_n_k2 = sum(k >= 2L))
}
t0 <- tic()
gw_sub_point <- gw_sub_rates(function(j) idx_by_ct[[j]])
toc("gw_subsample_point", t0)

## ---------------------------------------------------------------------------
## (b) Sample-level bootstrap  ->  bootstrap_sample_replicates.csv
## ---------------------------------------------------------------------------
maj_full <- maj
boot_replicate <- function(b) {
  set.seed(seed_boot_base + b)
  idx_b <- lapply(idx_by_ct, function(idx) idx[sample.int(length(idx), length(idx), replace = TRUE)])
  k_b <- integer(P_hv); npos <- integer(P_hv); nneg <- integer(P_hv); maxabs <- numeric(P_hv)
  s_w <- numeric(P_hv); s_wz <- numeric(P_hv); s_wz2 <- numeric(P_hv); Sxy <- numeric(P_hv)
  MUXb <- matrix(0, G, C); MUYb <- matrix(0, M, C); SSXb <- matrix(0, G, C); SSYb <- matrix(0, M, C)
  for (j in seq_len(C)) {
    idx <- idx_b[[j]]
    px <- prep_rows(X_hv[, idx, drop = FALSE]); py <- prep_rows(Y_hv[, idx, drop = FALSE])
    r <- as.vector(tcrossprod(px$Xs, py$Xs))
    W <- r * as.vector(outer(px$ss, py$ss)); W[!is.finite(W)] <- 0
    Sxy <- Sxy + W
    MUXb[, j] <- px$mu; MUYb[, j] <- py$mu; SSXb[, j] <- px$ss; SSYb[, j] <- py$ss
    if (length(idx) >= min_n_per_cohort) {
      v <- is.finite(r) & abs(r) < 1
      k_b <- k_b + v; npos <- npos + (v & r >= eps); nneg <- nneg + (v & r <= -eps)
      a <- abs(r); a[!v] <- 0; maxabs <- pmax(maxabs, a)
      w <- length(idx) - 3; z <- contextcorr_fisher_z(r); z[!v] <- 0
      s_w <- s_w + w * v; s_wz <- s_wz + w * z; s_wz2 <- s_wz2 + w * z^2
    }
    rm(px, py, r, W)
  }
  rg_b <- pooled_r_from_parts(Sxy, MUXb, MUYb, SSXb, SSYb, n_by_ct)
  sg_b <- sign(rg_b); sg_b[!is.finite(rg_b)] <- 0
  maj_b <- sign(npos - nneg)
  elig <- k_b >= min_k & is.finite(rg_b) & abs(rg_b) >= t_strong & maj_b != 0 & sg_b != 0
  rev_b <- elig & sg_b != maj_b
  elig0 <- k_b >= min_k & maj_b != 0 & sg_b != 0
  mixed_b <- npos > 0L & nneg > 0L
  em <- k_b >= min_k & maxabs >= t_strong
  dfb <- k_b - 1L
  Qb <- s_wz2 - s_wz^2 / pmax(s_w, 1e-12)
  okb <- dfb >= 1L & is.finite(Qb) & Qb >= 0
  I2b <- pmax(0, (Qb[okb] - dfb[okb]) / pmax(Qb[okb], 1e-12))
  pbv <- pchisq(Qb[okb], df = dfb[okb], lower.tail = FALSE)
  gws <- gw_sub_rates(function(j) idx_b[[j]])
  stats <- data.frame(b = b,
    frac_simpson_rg_ge_0.2 = sum(rev_b) / sum(elig), n_elig_simpson = sum(elig),
    frac_simpson_all = sum(elig0 & sg_b != maj_b) / sum(elig0),
    frac_mixed_maxabs_ge_0.2 = mean(mixed_b[em]), n_elig_mixed = sum(em),
    frac_mixed_all_hv = mean(mixed_b[k_b >= min_k]),
    median_I2 = median(I2b), I2_q25 = unname(quantile(I2b, .25)), I2_q75 = unname(quantile(I2b, .75)),
    frac_fdr_lt_0.05 = mean(p.adjust(pbv, method = "BH") < 0.05),
    n_pairs_k_lt_min_k = sum(k_b < min_k),
    gw_sub_mixed_eps0_all = unname(gws["gw_eps0_all"]), gw_sub_mixed_eps0_k2 = unname(gws["gw_eps0_k2"]),
    gw_sub_mixed_eps05_all = unname(gws["gw_eps05_all"]), gw_sub_mixed_eps05_k2 = unname(gws["gw_eps05_k2"]))
  list(stats = stats, maj_b = as.integer(maj_b), rev_b = rev_b, elig_b = elig)
}
t0 <- tic()
res_list <- pmap(chunk_seq(B_samp, n_cores), function(bs) {
  n_same <- integer(P_hv); n_rev <- integer(P_hv); n_elig <- integer(P_hv); stats <- vector("list", length(bs))
  for (i in seq_along(bs)) {
    r <- boot_replicate(bs[i])
    n_same <- n_same + (r$maj_b == maj_full); n_rev <- n_rev + r$rev_b; n_elig <- n_elig + r$elig_b
    stats[[i]] <- r$stats
    rm(r); gc(verbose = FALSE)
  }
  list(stats = do.call(rbind, stats), n_same = n_same, n_rev = n_rev, n_elig = n_elig)
})
boot_stats <- do.call(rbind, lapply(res_list, `[[`, "stats"))
boot_stats <- boot_stats[order(boot_stats$b), ]
n_same_maj <- Reduce(`+`, lapply(res_list, `[[`, "n_same"))
n_rev_boot <- Reduce(`+`, lapply(res_list, `[[`, "n_rev"))
n_elig_boot <- Reduce(`+`, lapply(res_list, `[[`, "n_elig"))
rm(res_list); gc(verbose = FALSE)
wcsv(boot_stats, "bootstrap_sample_replicates.csv")
gw_sub_check <- abs(gw_sub_point[["gw_eps0_all"]] - gw_rate_eps0_all)
toc("b_sample_level_bootstrap", t0)

## ---------------------------------------------------------------------------
## (b) Pair-level bootstrap  ->  bootstrap_cis.csv
## ---------------------------------------------------------------------------
t0 <- tic()
pair_boot <- function(x, B, FUN) {
  n <- length(x)
  unlist(pmap(chunk_seq(B, n_cores), function(bs) vapply(bs, function(i) { set.seed(seed_pair_base + i); FUN(x[sample.int(n, n, replace = TRUE)]) }, numeric(1))), use.names = FALSE)
}
pb_simpson <- pair_boot(reversal[elig_simpson], B_pair, mean)
pb_mixed <- pair_boot(mixed[elig_mixed], B_pair, mean)
pb_I2 <- pair_boot(I2[ok_het], B_pair, median)
set.seed(seed_main + 1L)
pb_gw <- rbinom(B_pair, size = P_gw, prob = gw_rate_eps0_all) / P_gw
pb_gw_eps05 <- rbinom(B_pair, size = P_gw, prob = gw_rate_eps05_all) / P_gw

ci <- function(x) unname(quantile(x, c(.025, .975), type = 7))
mk_ci <- function(id, quantity, domain, point, n_units, type, B, boots, note) {
  q <- ci(boots)
  data.frame(row_id = id, quantity = quantity, domain = domain, point_estimate = point, n_units = n_units,
             bootstrap_type = type, B = B, ci_lower = q[1], ci_upper = q[2],
             ci_lower_basic = 2 * point - q[2], ci_upper_basic = 2 * point - q[1],
             boot_mean = mean(boots), boot_bias = mean(boots) - point, boot_sd = sd(boots),
             note = note, stringsAsFactors = FALSE)
}
sb <- boot_stats
boot_cis <- rbind(
  mk_ci("B01", "simpson_reversal_rate_rg_ge_0.2_eps_0.05", "high-variance 10^6 pairs", rate_simpson, n_elig_simpson, "pair-level", B_pair, pb_simpson, "resample eligible pairs with replacement"),
  mk_ci("B02", "simpson_reversal_rate_rg_ge_0.2_eps_0.05", "high-variance 10^6 pairs", rate_simpson, n_elig_simpson, "sample-level", B_samp, sb$frac_simpson_rg_ge_0.2, "resample tumors within cohort; r_c, r_global and eligibility recomputed per replicate"),
  mk_ci("B03", "mixed_sign_rate_maxabs_ge_0.2_eps_0.05", "high-variance 10^6 pairs", rate_mixed_hv, n_elig_mixed, "pair-level", B_pair, pb_mixed, "resample eligible pairs with replacement"),
  mk_ci("B04", "mixed_sign_rate_maxabs_ge_0.2_eps_0.05", "high-variance 10^6 pairs", rate_mixed_hv, n_elig_mixed, "sample-level", B_samp, sb$frac_mixed_maxabs_ge_0.2, "resample tumors within cohort; eligibility (max_c |r_c| >= 0.2) recomputed per replicate"),
  mk_ci("B05", "mixed_sign_rate_genome_wide_eps_0", "genome-wide 23,170,038 pairs", gw_rate_eps0_all, P_gw, "pair-level", B_pair, pb_gw, "Binomial(n, p_hat)/n, exactly the pair-level bootstrap distribution of a proportion"),
  mk_ci("B06", "mixed_sign_rate_genome_wide_eps_0", "fixed random 100,000-pair subsample of the genome-wide domain", unname(gw_sub_point["gw_eps0_all"]), n_gw_sub, "sample-level", B_samp, sb$gw_sub_mixed_eps0_all, "resample tumors within cohort; the full 23M domain per replicate is infeasible, so a fixed 100,000-pair subsample (seed 20260903) is used; denominator = all subsample pairs (paper convention)"),
  mk_ci("B07", "mixed_sign_rate_genome_wide_eps_0.05", "genome-wide 23,170,038 pairs", gw_rate_eps05_all, P_gw, "pair-level", B_pair, pb_gw_eps05, "Binomial(n, p_hat)/n; eps = 0.05 variant, denominator all pairs"),
  mk_ci("B08", "mixed_sign_rate_genome_wide_eps_0.05", "fixed random 100,000-pair subsample of the genome-wide domain", unname(gw_sub_point["gw_eps05_all"]), n_gw_sub, "sample-level", B_samp, sb$gw_sub_mixed_eps05_all, "eps = 0.05 variant on the same subsample and replicates"),
  mk_ci("B09", "median_I2", "high-variance 10^6 pairs", median_I2, sum(ok_het), "pair-level", B_pair, pb_I2, "resample pairs with replacement, median of I2"),
  mk_ci("B10", "median_I2", "high-variance 10^6 pairs", median_I2, sum(ok_het), "sample-level", B_samp, sb$median_I2, "resample tumors within cohort; Q and I2 recomputed per replicate"),
  mk_ci("B11", "frac_pairs_Q_fdr_lt_0.05", "high-variance 10^6 pairs", frac_fdr, sum(ok_het), "sample-level", B_samp, sb$frac_fdr_lt_0.05, "BH across the 10^6 pairs within each replicate"),
  mk_ci("B12", "simpson_reversal_rate_all_rg_eps_0.05", "high-variance 10^6 pairs", rate_simpson_t0, sum(elig_simpson_t0), "sample-level", B_samp, sb$frac_simpson_all, "no |r_global| threshold (paper Fig. 3B at threshold 0)")
)
wcsv(boot_cis, "bootstrap_cis.csv")
toc("b_pair_level_bootstrap", t0)

## ---------------------------------------------------------------------------
## (c) Alternative reversal rules  ->  reversal_rules.csv, rule_agreement.csv, bootstrap_stability.csv
## ---------------------------------------------------------------------------
t0 <- tic()
stability <- n_same_maj / B_samp
stable <- stability >= stab_thresh
rev_flag_stability <- n_rev_boot / B_samp

D0 <- k_valid >= min_k & abs(rg) >= t_strong & sg != 0
n_D0 <- sum(D0)
label_from_vote <- function(m) {
  lab <- rep(NA_character_, P_hv)
  lab[D0 & m != 0 & sg != m] <- "reversal"
  lab[D0 & m != 0 & sg == m] <- "concordant"
  lab[D0 & m == 0] <- "undetermined"
  lab
}
lab_base <- label_from_vote(maj)
lab_stab <- lab_base; lab_stab[D0 & !stable] <- "undetermined"
rules <- list(
  list(id = "R01", rule = "i_baseline_majority_eps0.05", desc = "sign(sum_c s_c), s_c from eps = 0.05 (paper definition)", lab = lab_base),
  list(id = "R02", rule = "ii_precision_weighted_vote", desc = "sign(sum_c (N_c - 3) s_c), s_c from eps = 0.05", lab = label_from_vote(maj_w)),
  list(id = "R03", rule = "iii_random_effects_direction", desc = "sign(r_RE) from DerSimonian-Laird pooling of Fisher z", lab = label_from_vote(dir_RE)),
  list(id = "R04", rule = "iv_CI_gated_majority", desc = "s_c = sign(r_c) only if |z_c| > 1.96/sqrt(N_c - 3); then majority", lab = label_from_vote(maj_ci)),
  list(id = "R05", rule = "v_bootstrap_stable_majority", desc = sprintf("baseline label kept only if the majority sign is unchanged in >= %.0f%% of %d sample-level bootstrap replicates; else undetermined", 100 * stab_thresh, B_samp), lab = lab_stab),
  list(id = "R06", rule = "nominal_p_gated_majority", desc = "s_c (eps = 0.05) counts only if the Fisher-z p_c < 0.05 (two-sided); then majority", lab = label_from_vote(maj_p)),
  list(id = "R07", rule = "BH_q_gated_majority", desc = "s_c (eps = 0.05) counts only if the per-cohort BH q_c < 0.05 (BH across all measurable high-variance pairs within the cohort); then majority", lab = label_from_vote(maj_q))
)
reversal_rules <- do.call(rbind, lapply(rules, function(rl) {
  lab <- rl$lab[D0]; nr <- sum(lab == "reversal"); nc <- sum(lab == "concordant"); nu <- sum(lab == "undetermined")
  data.frame(row_id = rl$id, rule = rl$rule, description = rl$desc, n_domain = n_D0, n_determinate = nr + nc,
             n_reversal = nr, n_concordant = nc, n_undetermined = nu,
             frac_reversal_of_determinate = nr / (nr + nc), frac_reversal_of_domain = nr / n_D0,
             frac_reversal_of_baseline_eligible = nr / n_elig_simpson, stringsAsFactors = FALSE)
}))
wcsv(reversal_rules, "reversal_rules.csv")

base_lab <- lab_base[D0]
rule_agreement <- do.call(rbind, lapply(rules, function(rl) {
  lab <- rl$lab[D0]
  both <- base_lab != "undetermined" & lab != "undetermined"
  br <- base_lab == "reversal"
  data.frame(row_id = sub("^R", "A", rl$id), rule = rl$rule, n_domain = n_D0,
             frac_identical_3class = mean(lab == base_lab), kappa_3class = cohen_kappa(base_lab, lab),
             n_both_determinate = sum(both), frac_identical_2class = mean(lab[both] == base_lab[both]),
             kappa_2class = cohen_kappa(base_lab[both], lab[both]), n_baseline_reversals = sum(br),
             n_baseline_reversals_kept = sum(br & lab == "reversal"),
             frac_baseline_reversals_kept = sum(br & lab == "reversal") / sum(br),
             n_baseline_reversals_to_concordant = sum(br & lab == "concordant"),
             n_baseline_reversals_to_undetermined = sum(br & lab == "undetermined"),
             n_new_reversals = sum(!br & lab == "reversal"), stringsAsFactors = FALSE)
}))
wcsv(rule_agreement, "rule_agreement.csv")

stab_group <- function(id, name, sel) {
  s <- stability[sel]; f <- rev_flag_stability[sel]
  data.frame(row_id = id, group = name, n_pairs = sum(sel), frac_majority_stable_ge_0.9 = mean(s >= stab_thresh),
             median_majority_stability = median(s), q25_majority_stability = unname(quantile(s, .25)),
             q75_majority_stability = unname(quantile(s, .75)), mean_majority_stability = mean(s),
             frac_reversal_flag_ge_0.9 = mean(f >= stab_thresh), median_reversal_flag_frequency = median(f),
             mean_reversal_flag_frequency = mean(f), stringsAsFactors = FALSE)
}
boot_stability <- rbind(
  stab_group("T01", "baseline reversals (|r_global| >= 0.2)", reversal),
  stab_group("T02", "baseline concordant (|r_global| >= 0.2)", concordant),
  stab_group("T03", "all |r_global| >= 0.2 pairs (incl. ties)", D0),
  stab_group("T04", "all high-variance pairs", rep(TRUE, P_hv)))
wcsv(boot_stability, "bootstrap_stability.csv")
toc("c_reversal_rules", t0)

## ---------------------------------------------------------------------------
## Per-pair summaries and session info
## ---------------------------------------------------------------------------
t0 <- tic()
saveRDS(list(
  params = list(G = G, M = M, eps = eps, t_strong = t_strong, min_k = min_k, min_n_per_cohort = min_n_per_cohort,
                t_switch = t_switch, stab_thresh = stab_thresh, B_pair = B_pair, B_samp = B_samp,
                seed_main = seed_main, seed_boot_base = seed_boot_base, matched_file = matched_file),
  genes = genes, mirs = mirs, cohorts = cohorts, n_by_ct = n_by_ct, r_global = rg, k_valid = k_valid,
  maj = maj, maj_w = maj_w, maj_ci = maj_ci, maj_p = maj_p, maj_q = maj_q, r_RE = r_RE, r_fixed = r_fixed,
  I2 = I2, p_het = p_het, fdr_het = fdr_het, mixed = mixed, max_abs = max_abs, r_max = r_max, r_min = r_min,
  reversal = reversal, elig_simpson = elig_simpson, n_same_maj_boot = n_same_maj, n_rev_boot = n_rev_boot,
  n_elig_boot = n_elig_boot, switch_nominal = switch_nom, switch_bh_within_pair = switch_bh,
  gw_subsample_pair_index = gw_sub_p), file.path(out_dir, "pair_summaries.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
toc("save", t0)

total_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
timing_df <- data.frame(row_id = sprintf("X%02d", seq_along(timing)), step = names(timing),
                        seconds = round(unlist(timing), 1), stringsAsFactors = FALSE)
timing_df <- rbind(timing_df, data.frame(row_id = "X99", step = "total", seconds = round(total_secs, 1)))
wcsv(timing_df, "timing.csv")
