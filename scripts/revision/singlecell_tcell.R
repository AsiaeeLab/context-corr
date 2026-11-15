#!/usr/bin/env Rscript
# =============================================================================
# singlecell_tcell.R  (revision analysis for the BMC Genomics manuscript)
#
# (a) Corrected Fig. S3: Simpson-reversal rate vs |r_global| threshold for
#     coarse PBMC cell types (sc_100) vs T-cell subtypes (sc_110), eps grid.
#     For each panel, r_global and r_c use the same cells: only contexts with
#     at least N_MIN_MAIN cells enter the analysis population.
# (b) Controlled T-cell refinement comparisons:
#     (i)  same retained cells and same pair set; partition {coarse types} vs
#          {coarse types with T split into subtypes}; the common cell set is
#          fixed before either partition is analyzed;
#     (ii) T cells only, fixed pair set; pooled-over-T r vs subtype majority.
# (c) Sparse-aware check of the common-population residualized flip rates (PBMC
#     11.8%, CITE-seq joint WNN 37.1%): Spearman and aggregation versions next
#     to Pearson.
# (d) PBMC 3k donor / metadata check.
#
# Run from the repository root:
#   Rscript scripts/revision/singlecell_tcell.R
#
# Inputs (read only): results/sc/*/analysis_ready.rds, results/sc/*/
# reversal_rate_curve.csv, and the raw Seurat .rds files whose paths are
# recorded in those caches. Everything is written to
# results/revision/singlecell_tcell/ only.
#
# Only SeuratObject is needed to read the .rds files. LogNormalize and CLR are
# re-implemented below. Within-context and residualized Pearson quantities are
# checked against the cached pipeline output (to 1e-10); pooled quantities are
# recomputed on the retained-context population adopted for the revision.
# =============================================================================

t_start <- Sys.time()
suppressPackageStartupMessages({
  library(Matrix)
  stopifnot(requireNamespace("SeuratObject", quietly = TRUE))
})
source("00-paths.R")
source(file.path("R", "contextcorr.R"))
source(file.path("R", "sc_context_pipeline.R"))

## ---------------------------------------------------------------- settings --
SEED <- 20260903
B_BOOT <- 1000
FIGURE_ONLY <- identical(Sys.getenv("SINGLECELL_TCELL_FIGURE_ONLY", "0"), "1")
EPS_GRID <- c(0, 0.05, 0.1)
THR_GRID <- seq(0, 0.5, by = 0.05)
THR_REPORT <- c(0.05, 0.2)
N_MIN_MAIN <- 200          # paper convention: a context needs >= 200 cells
N_MIN_SENS <- 50           # sensitivity: lets the small T subtypes vote
PB_UNIT_SIZES <- c(20, 40) # cells per pseudobulk unit (within context)
PB_N_MIN_UNITS <- 5        # units a context needs to contribute a within-context r
SIMPSON_THR <- 0.2         # for the Simpson columns of sparse_check.csv
SIMPSON_EPS <- 0.05

res_root <- if (dir.exists(file.path(paths$results, "sc"))) paths$results else "results"
out_dir <- file.path("results", "revision", "singlecell_tcell")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
log_msg("results root used for caches: %s", normalizePath(res_root))
log_msg("output directory: %s", normalizePath(out_dir))

## ----------------------------------------------------------------- helpers --
check_close <- function(a, b, tol, label) {
  stopifnot(length(a) == length(b))
  if (!all(is.na(a) == is.na(b))) stop("Validation failed (NA pattern differs): ", label)
  d <- abs(a - b)
  mx <- if (any(!is.na(d))) max(d, na.rm = TRUE) else 0
  if (mx > tol) stop(sprintf("Validation failed: %s (max |diff| = %g)", label, mx))
  log_msg("check ok: %-60s max |diff| = %.1e over %d values", label, mx, sum(!is.na(d)))
  invisible(mx)
}

qc_keep_cells <- function(counts_rna, md, min_features = 200, max_mt = 20) {
  # Mirrors sc_preprocess(): nFeature_RNA >= 200 and percent.mt <= 20 ("^MT-").
  mt <- grepl("^MT-", rownames(counts_rna))
  pct_mt <- 100 * Matrix::colSums(counts_rna[mt, , drop = FALSE]) / Matrix::colSums(counts_rna)
  nfeat <- if ("nFeature_RNA" %in% colnames(md)) md$nFeature_RNA else Matrix::colSums(counts_rna > 0)
  list(keep = nfeat >= min_features & pct_mt <= max_mt, pct_mt = pct_mt, nfeat = nfeat)
}

log_normalize_sparse <- function(counts, scale_factor = 1e4) {
  # Seurat LogNormalize: log1p(count / total_counts_of_cell * 1e4)
  stopifnot(inherits(counts, "dgCMatrix"))
  cs <- Matrix::colSums(counts)
  out <- counts
  out@x <- log1p(out@x / rep.int(cs, diff(out@p)) * scale_factor)
  out
}

log_normalize_dense_subset <- function(counts_all, features, scale_factor = 1e4) {
  # Same normalization for a dense subset; totals are taken over ALL features.
  cs <- Matrix::colSums(counts_all)
  m <- as.matrix(counts_all[features, , drop = FALSE])
  log1p(sweep(m, 2, cs, "/") * scale_factor)
}

clr_by_cell <- function(m) {
  # Seurat CLR with margin = 2: per cell, log1p(x / exp(mean(log1p(x)))).
  m <- as.matrix(m)
  out <- apply(m, 2, function(x) log1p(x / exp(mean(log1p(x)))))
  out <- matrix(out, nrow = nrow(m), dimnames = dimnames(m))
  out
}

rank_rows <- function(m) {
  r <- t(apply(m, 1, rank, ties.method = "average"))
  r <- matrix(r, nrow = nrow(m), dimnames = dimnames(m))
  r
}

pair_cor_from_z <- function(zx, zy, idx_x, idx_y, chunk = 4000L) {
  n <- length(idx_x)
  out <- rep(NA_real_, n)
  if (n == 0) return(out)
  for (a in seq(1L, n, by = chunk)) {
    b <- min(n, a + chunk - 1L)
    ii <- a:b
    out[ii] <- sc_pairwise_cor_from_z(zx, zy, idx_x[ii], idx_y[ii])
  }
  out
}

pair_stats_generic <- function(X, Y = NULL, pairs, ctx, n_min, method = c("pearson", "spearman"),
                               resid_mode = c("residualize_then_rank", "rank_then_residualize")) {
  # Common revision convention: r_global, within-context r, and r_resid all use
  # cells from contexts with >= n_min cells. Spearman = Pearson on ranks (global:
  # ranks over all retained cells;
  # within context: ranks within the context; residualized: ranks of the
  # mean-residualized values), which matches cor(method = "spearman") and
  # mean_residualized_cor(method = "spearman"). resid_mode = "rank_then_residualize"
  # is the alternative rank-based version: rank over the retained cells first,
  # center the ranks within context, then Pearson (only r_resid differs).
  method <- match.arg(method)
  resid_mode <- match.arg(resid_mode)
  symmetric <- is.null(Y)
  if (symmetric) Y <- X
  stopifnot(is.matrix(X), is.matrix(Y), ncol(X) == ncol(Y), length(ctx) == ncol(X))
  ctx <- as.character(ctx)
  tab <- table(ctx)
  levels_use <- names(tab)[tab >= n_min & tab >= 3]
  keep <- which(ctx %in% levels_use)
  f <- if (method == "spearman") rank_rows else identity

  zx <- sc_row_zscore_matrix(f(X[, keep, drop = FALSE]))
  zy <- if (symmetric) zx else sc_row_zscore_matrix(f(Y[, keep, drop = FALSE]))
  r_global <- pair_cor_from_z(zx, zy, pairs$idx_x, pairs$idx_y)

  r_ctx <- matrix(NA_real_, nrow(pairs), length(levels_use), dimnames = list(NULL, levels_use))
  for (lv in levels_use) {
    cells <- which(ctx == lv)
    zxc <- sc_row_zscore_matrix(f(X[, cells, drop = FALSE]))
    zyc <- if (symmetric) zxc else sc_row_zscore_matrix(f(Y[, cells, drop = FALSE]))
    r_ctx[, lv] <- pair_cor_from_z(zxc, zyc, pairs$idx_x, pairs$idx_y)
  }

  r_resid <- rep(NA_real_, nrow(pairs))
  if (length(keep) >= 4) {
    fk <- factor(ctx[keep])
    if (method == "spearman" && resid_mode == "rank_then_residualize") {
      Xr <- sc_residualize_matrix_by_context(rank_rows(X[, keep, drop = FALSE]), fk)
      Yr <- if (symmetric) Xr else sc_residualize_matrix_by_context(rank_rows(Y[, keep, drop = FALSE]), fk)
      zxr <- sc_row_zscore_matrix(Xr)
      zyr <- if (symmetric) zxr else sc_row_zscore_matrix(Yr)
    } else {
      Xr <- sc_residualize_matrix_by_context(X[, keep, drop = FALSE], fk)
      Yr <- if (symmetric) Xr else sc_residualize_matrix_by_context(Y[, keep, drop = FALSE], fk)
      zxr <- sc_row_zscore_matrix(f(Xr))
      zyr <- if (symmetric) zxr else sc_row_zscore_matrix(f(Yr))
    }
    r_resid <- pair_cor_from_z(zxr, zyr, pairs$idx_x, pairs$idx_y)
  }
  method_label <- if (method == "spearman" && resid_mode == "rank_then_residualize") "spearman_rankfirst" else method

  list(r_global = r_global, r_ctx = r_ctx, r_resid = r_resid,
       levels_use = levels_use, n_by_level = as.integer(tab[levels_use]),
       context_table = tab, n_cells = length(keep), n_cells_resid = length(keep),
       method = method_label, n_min = n_min)
}

sign_flags <- function(r_global, r_ctx, eps) {
  # eps rule of the paper: s_c = +1 if r_c >= eps, -1 if r_c <= -eps, else 0.
  pos <- rowSums(r_ctx >= eps, na.rm = TRUE)
  neg <- rowSums(r_ctx <= -eps, na.rm = TRUE)
  maj <- sign(pos - neg)
  gs <- sign(r_global)
  gs[!is.finite(r_global)] <- 0
  list(pos = pos, neg = neg, maj = maj,
       eligible = is.finite(r_global) & gs != 0 & maj != 0,
       reversal = is.finite(r_global) & gs != 0 & maj != 0 & gs != maj,
       mixed = pos > 0 & neg > 0)
}

flip_resid <- function(r_global, r_resid) {
  is.finite(r_global) & is.finite(r_resid) & r_global != 0 & r_resid != 0 & sign(r_global) != sign(r_resid)
}

boot_rate_ci <- function(x, B = B_BOOT, seed = SEED) {
  x <- as.logical(x)
  n <- length(x)
  if (n == 0) return(c(NA_real_, NA_real_))
  set.seed(seed)
  m <- vapply(seq_len(B), function(b) mean(x[sample.int(n, n, replace = TRUE)]), numeric(1))
  as.numeric(stats::quantile(m, c(0.025, 0.975)))
}

boot_paired_ci <- function(a, b, B = B_BOOT, seed = SEED) {
  a <- as.logical(a); b <- as.logical(b)
  n <- length(a)
  if (n == 0) return(list(a = c(NA_real_, NA_real_), b = c(NA_real_, NA_real_), diff = c(NA_real_, NA_real_)))
  set.seed(seed)
  ra <- numeric(B); rb <- numeric(B)
  for (i in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    ra[i] <- mean(a[idx]); rb[i] <- mean(b[idx])
  }
  q <- function(v) as.numeric(stats::quantile(v, c(0.025, 0.975)))
  list(a = q(ra), b = q(rb), diff = q(ra - rb))
}

exact_ci <- function(k, n) {
  if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
  as.numeric(stats::binom.test(k, n)$conf.int)
}

ctx_summary <- function(ctx, n_min) {
  tab <- table(ctx)
  used <- tab[tab >= n_min]
  excl <- tab[tab < n_min]
  fmt <- function(t) if (length(t)) paste(sprintf("%s(%d)", names(t), as.integer(t)), collapse = ";") else ""
  list(used = fmt(used), excluded = fmt(excl), k = length(used))
}

make_pseudobulk <- function(counts_list, ctx, n_min_cells, m, seed) {
  # Within each context with >= n_min_cells cells, shuffle cells (fixed seed)
  # and cut them into floor(n_c / m) units of m or m + 1 cells (all cells of
  # the context are used). Raw counts are summed per unit.
  ctx <- as.character(ctx)
  tab <- table(ctx)
  use <- names(tab)[tab >= n_min_cells]
  set.seed(seed)
  unit <- rep(NA_character_, length(ctx))
  for (lv in use) {
    cells <- which(ctx == lv)
    cells <- cells[sample.int(length(cells))]
    n_units <- max(1L, floor(length(cells) / m))
    grp <- if (n_units == 1L) rep(1L, length(cells)) else cut(seq_along(cells), breaks = n_units, labels = FALSE)
    unit[cells] <- sprintf("%s#%03d", lv, grp)
  }
  ok <- !is.na(unit)
  u <- factor(unit[ok])
  ind <- Matrix::t(Matrix::fac2sparse(u))  # cells(ok) x units
  agg <- lapply(counts_list, function(cm) {
    a <- cm[, ok, drop = FALSE] %*% ind
    a <- as(a, "CsparseMatrix")
    colnames(a) <- levels(u)
    a
  })
  unit_ctx <- sub("#.*$", "", levels(u))
  list(agg = agg, unit_context = unit_ctx, n_units = length(levels(u)),
       units_per_context = table(unit_ctx), n_cells_used = sum(ok))
}

## ------------------------------------------------------------ load caches --
ar100 <- readRDS(file.path(res_root, "sc", "sc_100_pbmc_coarse_celltypes", "analysis_ready.rds"))
ar110 <- readRDS(file.path(res_root, "sc", "sc_110_pbmc_tcell_subtypes", "analysis_ready.rds"))
ar200 <- readRDS(file.path(res_root, "sc", "sc_200_citeseq_protein_rna_contexts", "analysis_ready.rds"))
# =============================================================================
# Load PBMC 3k (Seurat-free), replicate QC, normalize, and validate vs caches
# =============================================================================
log_msg("loading PBMC 3k and replicating the pipeline preprocessing")
pbmc_src <- sub("^local_rds(_list)?:", "", ar100$inputs$extra$dataset$source)
obj <- readRDS(pbmc_src)
cnt <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
md <- obj@meta.data
stopifnot(identical(colnames(cnt), rownames(md)))
qc <- qc_keep_cells(cnt, md)
stopifnot(sum(qc$keep) == length(ar100$context))
cnt_qc <- cnt[, qc$keep, drop = FALSE]
ln <- log_normalize_sparse(cnt_qc)

X100 <- as.matrix(ln[ar100$inputs$feature_x, , drop = FALSE])
ctx_coarse <- as.character(ar100$context)
is_T <- ctx_coarse == "T"
stopifnot(sum(is_T) == length(ar110$context))
ctx_tsub <- as.character(ar110$context)
ctx_split <- ctx_coarse
ctx_split[is_T] <- paste0("T_", ctx_tsub)
X110 <- as.matrix(ln[ar110$inputs$feature_x, is_T, drop = FALSE])
pairs100 <- ar100$pairs
pairs110 <- ar110$pairs

st_coarse <- pair_stats_generic(X100, NULL, pairs100, ctx_coarse, N_MIN_MAIN, "pearson")
for (lv in st_coarse$levels_use) check_close(st_coarse$r_ctx[, lv], ar100$r_by_context[, lv], 1e-10, sprintf("PBMC r_c[%s] vs sc_100 cache", lv))
check_close(st_coarse$r_resid, ar100$summary_pairs$r_resid, 1e-10, "PBMC r_resid vs sc_100 cache")

st_t110 <- pair_stats_generic(X110, NULL, pairs110, ctx_tsub, N_MIN_MAIN, "pearson")
for (lv in st_t110$levels_use) check_close(st_t110$r_ctx[, lv], ar110$r_by_context[, lv], 1e-10, sprintf("T-cell r_c[%s] vs sc_110 cache", lv))
check_close(st_t110$r_resid, ar110$summary_pairs$r_resid, 1e-10, "T-cell r_resid vs sc_110 cache")

# =============================================================================
# (a) Fig. S3 (corrected) + tcell_rates.csv
# =============================================================================
log_msg("(a) reversal-rate curves on common retained-context populations")
keep100 <- ctx_coarse %in% st_coarse$levels_use
keep110 <- ctx_tsub %in% st_t110$levels_use
st_coarse_uniform <- pair_stats_generic(X100[, keep100, drop = FALSE], NULL, pairs100,
                                        ctx_coarse[keep100], N_MIN_MAIN, "pearson")
st_t110_uniform <- pair_stats_generic(X110[, keep110, drop = FALSE], NULL, pairs110,
                                      ctx_tsub[keep110], N_MIN_MAIN, "pearson")
for (lv in st_coarse$levels_use) {
  check_close(st_coarse_uniform$r_ctx[, lv], st_coarse$r_ctx[, lv], 1e-10,
              sprintf("uniform PBMC r_c[%s] vs cached-context value", lv))
}
for (lv in st_t110$levels_use) {
  check_close(st_t110_uniform$r_ctx[, lv], st_t110$r_ctx[, lv], 1e-10,
              sprintf("uniform T-cell r_c[%s] vs cached-context value", lv))
}

c100 <- sc_reversal_rate_curve(st_coarse_uniform$r_global, st_coarse_uniform$r_ctx,
                               THR_GRID, EPS_GRID)
c110 <- sc_reversal_rate_curve(st_t110_uniform$r_global, st_t110_uniform$r_ctx,
                               THR_GRID, EPS_GRID)
cs100 <- ctx_summary(ctx_coarse, N_MIN_MAIN)
cs110 <- ctx_summary(ctx_tsub, N_MIN_MAIN)
rates <- rbind(
  data.frame(analysis = "coarse_cell_types", experiment_id = "sc_100_pbmc_coarse_celltypes",
             n_cells = sum(keep100), n_pairs_total = nrow(pairs100),
             contexts_used = cs100$used, contexts_excluded_below_n_min = cs100$excluded,
             analysis_population = "cells in retained contexts", c100, stringsAsFactors = FALSE),
  data.frame(analysis = "t_cell_subtypes", experiment_id = "sc_110_pbmc_tcell_subtypes",
             n_cells = sum(keep110), n_pairs_total = nrow(pairs110),
             contexts_used = cs110$used, contexts_excluded_below_n_min = cs110$excluded,
             analysis_population = "cells in retained contexts", c110, stringsAsFactors = FALSE)
)
rates$n_reversal <- as.integer(round(rates$reversal_rate * rates$eligible_pairs))
rates$low_count_flag <- rates$eligible_pairs < 50
utils::write.csv(rates, file.path(out_dir, "tcell_rates.csv"), row.names = FALSE)

fig_file <- file.path(out_dir, "Figure_S3_TcellRefinement.pdf")
eps_cols <- c("0" = "#1B4F72", "0.05" = "#E67E22", "0.1" = "#7F8C8D")
panel_titles <- c(coarse_cell_types = "Coarse Cell Types", t_cell_subtypes = "T-Cell Subtypes")
ymax <- max(0.30, ceiling(max(rates$reversal_rate, na.rm = TRUE) * 20) / 20)
grDevices::cairo_pdf(fig_file, width = 8.5, height = 4.0)
par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3.0, 1.0), font.main = 2,
    cex.main = 1.15, cex.lab = 1.05, las = 1)
for (i in seq_along(panel_titles)) {
  d <- rates[rates$analysis == names(panel_titles)[i], ]
  plot(NA, xlim = c(0, 0.5), ylim = c(0, ymax),
       xlab = expression(paste("|", r[global], "| threshold")),
       ylab = "Simpson-reversal rate", main = panel_titles[[i]])
  grid(col = "grey92", lty = 1)
  for (e in EPS_GRID) {
    s <- d[d$eps == e, ]
    s <- s[order(s$abs_r_global_threshold), ]
    col_e <- eps_cols[[as.character(e)]]
    lines(s$abs_r_global_threshold, s$reversal_rate, lwd = 2, col = col_e)
    points(s$abs_r_global_threshold, s$reversal_rate,
           pch = ifelse(s$low_count_flag, 1, 16), col = col_e, cex = 1.0, lwd = 1.6)
  }
  usr <- par("usr")
  mtext(paste0(LETTERS[i], "."), side = 3, line = 1.1, adj = 0,
        at = usr[1] - 0.19 * diff(usr[1:2]), font = 2, cex = 1.35)
  legend("topright", legend = as.expression(lapply(EPS_GRID, function(e) bquote(epsilon == .(e)))),
         col = eps_cols, lwd = 2, pch = 16, bty = "n", cex = 0.92)
  if (any(d$low_count_flag)) {
    legend("right", legend = "fewer than 50 eligible pairs", pch = 1,
           pt.lwd = 1.6, pt.cex = 1.0, bty = "n", cex = 0.85)
  }
}
invisible(dev.off())
log_msg("wrote %s", fig_file)

if (FIGURE_ONLY) {
  log_msg("figure-only mode complete")
  quit(save = "no", status = 0L)
}

# Spearman implementation check against the canonical helpers on 200 pairs.
st_coarse_sp <- pair_stats_generic(X100, NULL, pairs100, ctx_coarse, N_MIN_MAIN, "spearman")
set.seed(SEED)
vi <- sample.int(nrow(pairs100), 200)
ref_g <- vapply(vi, function(i) suppressWarnings(stats::cor(X100[pairs100$idx_x[i], keep100], X100[pairs100$idx_y[i], keep100], method = "spearman")), numeric(1))
check_close(st_coarse_sp$r_global[vi], ref_g, 1e-10, "retained-population Spearman r_global vs cor(method = 'spearman')")
ref_T <- vapply(vi, function(i) suppressWarnings(stats::cor(X100[pairs100$idx_x[i], is_T], X100[pairs100$idx_y[i], is_T], method = "spearman")), numeric(1))
check_close(st_coarse_sp$r_ctx[vi, "T"], ref_T, 1e-10, "Spearman r_c[T] vs cor(method = 'spearman')")
ctx_ok <- ctx_coarse
ctx_ok[!(ctx_ok %in% st_coarse_sp$levels_use)] <- NA
ref_r <- vapply(vi, function(i) mean_residualized_cor(X100[pairs100$idx_x[i], ], X100[pairs100$idx_y[i], ], ctx_ok, method = "spearman"), numeric(1))
check_close(st_coarse_sp$r_resid[vi], ref_r, 1e-10, "Spearman r_resid vs mean_residualized_cor(method = 'spearman')")

# =============================================================================
# (b) Controlled T-cell comparisons
# =============================================================================
log_msg("(b) controlled comparisons")
CTRL_COLS <- c("design", "pair_universe", "cells", "n_cells", "partition", "n_min", "k_contexts_used",
               "contexts_used", "contexts_excluded_below_n_min", "n_cells_resid", "eps", "threshold",
               "n_pairs_strong", "n_eligible_partition", "n_reversal_partition_all_eligible",
               "reversal_rate_partition_all_eligible", "n_eligible_common",
               "n_reversal", "reversal_rate", "reversal_boot_lo", "reversal_boot_hi", "reversal_exact_lo", "reversal_exact_hi",
               "n_mixed", "mixed_rate", "mixed_boot_lo", "mixed_boot_hi",
               "n_flip_resid", "flip_resid_rate",
               "n_rev_coarse_only", "n_rev_split_only", "n_rev_both", "note")
ctrl_row <- function(...) {
  v <- list(...)
  base <- setNames(as.list(rep(NA, length(CTRL_COLS))), CTRL_COLS)
  base[names(v)] <- v
  as.data.frame(base, stringsAsFactors = FALSE)
}
ctrl <- list()

# (i) same retained cells, same pairs; partition coarse vs coarse-with-T-split.
# Fix the analysis population first to cells in split-partition contexts meeting
# n_min, then apply both partitions to that identical population. This prevents
# the split partition's smaller T-cell contexts from changing the cell set.
controlled_stats_by_n <- list()
for (n_min in c(N_MIN_MAIN, N_MIN_SENS)) {
  split_tab_full <- table(ctx_split)
  split_levels_keep <- names(split_tab_full)[split_tab_full >= n_min & split_tab_full >= 3]
  keep_controlled <- ctx_split %in% split_levels_keep
  X_controlled <- X100[, keep_controlled, drop = FALSE]
  ctx_coarse_controlled <- ctx_coarse[keep_controlled]
  ctx_split_controlled <- ctx_split[keep_controlled]
  sc_ <- pair_stats_generic(X_controlled, NULL, pairs100, ctx_coarse_controlled, n_min, "pearson")
  ss_ <- pair_stats_generic(X_controlled, NULL, pairs100, ctx_split_controlled, n_min, "pearson")
  stopifnot(isTRUE(all.equal(sc_$r_global, ss_$r_global)))
  stopifnot(sc_$n_cells_resid == ncol(X_controlled), ss_$n_cells_resid == ncol(X_controlled))
  controlled_stats_by_n[[as.character(n_min)]] <- list(
    coarse = sc_, split = ss_, keep = keep_controlled,
    ctx_coarse = ctx_coarse_controlled, ctx_split = ctx_split_controlled)
  r_g <- sc_$r_global
  csc <- ctx_summary(ctx_coarse_controlled, n_min)
  css <- ctx_summary(ctx_split_controlled, n_min)
  for (eps in EPS_GRID) {
    fc <- sign_flags(r_g, sc_$r_ctx, eps)
    fs <- sign_flags(r_g, ss_$r_ctx, eps)
    frc <- flip_resid(r_g, sc_$r_resid)
    frs <- flip_resid(r_g, ss_$r_resid)
    for (thr in THR_REPORT) {
      strong <- is.finite(r_g) & abs(r_g) >= thr & r_g != 0
      common <- strong & fc$eligible & fs$eligible
      idx <- which(common)
      n <- length(idx)
      bp <- boot_paired_ci(fc$reversal[idx], fs$reversal[idx])
      bm <- boot_paired_ci(fc$mixed[idx], fs$mixed[idx])
      kc <- sum(fc$reversal[idx]); ks <- sum(fs$reversal[idx])
      common_args <- list(design = "i_same_cells_same_pairs", pair_universe = "sc_100_pairs", cells = "common_retained_pbmc",
                          n_cells = ncol(X_controlled), n_min = n_min, eps = eps, threshold = thr,
                          n_pairs_strong = sum(strong), n_eligible_common = n)
      ctrl[[length(ctrl) + 1]] <- do.call(ctrl_row, c(common_args, list(
        partition = "coarse_6types", k_contexts_used = csc$k, contexts_used = csc$used,
        contexts_excluded_below_n_min = csc$excluded, n_cells_resid = sc_$n_cells_resid,
        n_eligible_partition = sum(strong & fc$eligible),
        n_reversal_partition_all_eligible = sum(fc$reversal[strong & fc$eligible]),
        reversal_rate_partition_all_eligible = mean(fc$reversal[strong & fc$eligible]),
        n_reversal = kc, reversal_rate = mean(fc$reversal[idx]), reversal_boot_lo = bp$a[1], reversal_boot_hi = bp$a[2],
        reversal_exact_lo = exact_ci(kc, n)[1], reversal_exact_hi = exact_ci(kc, n)[2],
        n_mixed = sum(fc$mixed[idx]), mixed_rate = mean(fc$mixed[idx]), mixed_boot_lo = bm$a[1], mixed_boot_hi = bm$a[2],
        n_flip_resid = sum(frc[idx]), flip_resid_rate = mean(frc[idx]),
        note = "rates on pairs eligible under BOTH partitions")))
      ctrl[[length(ctrl) + 1]] <- do.call(ctrl_row, c(common_args, list(
        partition = "coarse_with_T_split", k_contexts_used = css$k, contexts_used = css$used,
        contexts_excluded_below_n_min = css$excluded, n_cells_resid = ss_$n_cells_resid,
        n_eligible_partition = sum(strong & fs$eligible),
        n_reversal_partition_all_eligible = sum(fs$reversal[strong & fs$eligible]),
        reversal_rate_partition_all_eligible = mean(fs$reversal[strong & fs$eligible]),
        n_reversal = ks, reversal_rate = mean(fs$reversal[idx]), reversal_boot_lo = bp$b[1], reversal_boot_hi = bp$b[2],
        reversal_exact_lo = exact_ci(ks, n)[1], reversal_exact_hi = exact_ci(ks, n)[2],
        n_mixed = sum(fs$mixed[idx]), mixed_rate = mean(fs$mixed[idx]), mixed_boot_lo = bm$b[1], mixed_boot_hi = bm$b[2],
        n_flip_resid = sum(frs[idx]), flip_resid_rate = mean(frs[idx]),
        note = "rates on pairs eligible under BOTH partitions")))
      ctrl[[length(ctrl) + 1]] <- do.call(ctrl_row, c(common_args, list(
        partition = "difference_coarse_minus_split",
        n_reversal = kc - ks, reversal_rate = mean(fc$reversal[idx]) - mean(fs$reversal[idx]),
        reversal_boot_lo = bp$diff[1], reversal_boot_hi = bp$diff[2],
        n_mixed = sum(fc$mixed[idx]) - sum(fs$mixed[idx]), mixed_rate = mean(fc$mixed[idx]) - mean(fs$mixed[idx]),
        mixed_boot_lo = bm$diff[1], mixed_boot_hi = bm$diff[2],
        n_rev_coarse_only = sum(fc$reversal[idx] & !fs$reversal[idx]),
        n_rev_split_only = sum(!fc$reversal[idx] & fs$reversal[idx]),
        n_rev_both = sum(fc$reversal[idx] & fs$reversal[idx]),
        note = "paired bootstrap (same resampled pairs for both partitions)")))
    }
  }
}

# Pair-level table for design (i) at the headline setting (thr 0.2, n_min 200, eps 0.05).
{
  primary_controlled <- controlled_stats_by_n[[as.character(N_MIN_MAIN)]]
  st_controlled_coarse <- primary_controlled$coarse
  st_controlled_split <- primary_controlled$split
  r_g <- st_controlled_coarse$r_global
  fc <- sign_flags(r_g, st_controlled_coarse$r_ctx, 0.05)
  fs <- sign_flags(r_g, st_controlled_split$r_ctx, 0.05)
  strong <- which(is.finite(r_g) & abs(r_g) >= 0.2 & r_g != 0)
  pl <- data.frame(pair_id = pairs100$pair_id[strong], feature_x = pairs100$feature_x[strong], feature_y = pairs100$feature_y[strong],
                   r_global = r_g[strong], stringsAsFactors = FALSE)
  for (lv in st_controlled_coarse$levels_use) pl[[paste0("r_coarse_", lv)]] <- st_controlled_coarse$r_ctx[strong, lv]
  for (lv in setdiff(st_controlled_split$levels_use, st_controlled_coarse$levels_use)) pl[[paste0("r_split_", lv)]] <- st_controlled_split$r_ctx[strong, lv]
  pl$majority_coarse <- fc$maj[strong]; pl$majority_split <- fs$maj[strong]
  pl$eligible_coarse <- fc$eligible[strong]; pl$eligible_split <- fs$eligible[strong]
  pl$reversal_coarse <- fc$reversal[strong]; pl$reversal_split <- fs$reversal[strong]
  pl$mixed_coarse <- fc$mixed[strong]; pl$mixed_split <- fs$mixed[strong]
  pl$r_resid_coarse <- st_controlled_coarse$r_resid[strong]; pl$r_resid_split <- st_controlled_split$r_resid[strong]
  utils::write.csv(pl, file.path(out_dir, "tcell_controlled_pairs_thr0.2.csv"), row.names = FALSE)
}

# (ii) T cells only, fixed pair set; pooled-over-T r vs within-subtype majority
st_t110_50 <- pair_stats_generic(X110, NULL, pairs110, ctx_tsub, N_MIN_SENS, "pearson")
X100T <- X100[, is_T, drop = FALSE]
st_t100 <- pair_stats_generic(X100T, NULL, pairs100, ctx_tsub, N_MIN_MAIN, "pearson")
st_t100_50 <- pair_stats_generic(X100T, NULL, pairs100, ctx_tsub, N_MIN_SENS, "pearson")

tonly <- list(
  list(universe = "sc_110_pairs", st = list(`200` = st_t110, `50` = st_t110_50), n_feat = length(ar110$inputs$feature_x)),
  list(universe = "sc_100_pairs", st = list(`200` = st_t100, `50` = st_t100_50), n_feat = length(ar100$inputs$feature_x))
)
for (u in tonly) {
  for (n_min in c(N_MIN_MAIN, N_MIN_SENS)) {
    st <- u$st[[as.character(n_min)]]
    r_g <- st$r_global
    cs <- ctx_summary(ctx_tsub, n_min)
    for (eps in EPS_GRID) {
      f <- sign_flags(r_g, st$r_ctx, eps)
      fr <- flip_resid(r_g, st$r_resid)
      for (thr in THR_REPORT) {
        strong <- is.finite(r_g) & abs(r_g) >= thr & r_g != 0
        idx <- which(strong & f$eligible)
        n <- length(idx)
        k <- sum(f$reversal[idx])
        ctrl[[length(ctrl) + 1]] <- ctrl_row(
          design = "ii_tcells_only", pair_universe = u$universe, cells = "t_cells", n_cells = ncol(X100T),
          partition = "none_vs_t_subtypes", n_min = n_min, k_contexts_used = cs$k, contexts_used = cs$used,
          contexts_excluded_below_n_min = cs$excluded, n_cells_resid = st$n_cells_resid, eps = eps, threshold = thr,
          n_pairs_strong = sum(strong), n_eligible_partition = n, n_eligible_common = n,
          n_reversal_partition_all_eligible = k, reversal_rate_partition_all_eligible = if (n > 0) k / n else NA_real_,
          n_reversal = k, reversal_rate = if (n > 0) k / n else NA_real_,
          reversal_boot_lo = boot_rate_ci(f$reversal[idx])[1], reversal_boot_hi = boot_rate_ci(f$reversal[idx])[2],
          reversal_exact_lo = exact_ci(k, n)[1], reversal_exact_hi = exact_ci(k, n)[2],
          n_mixed = sum(f$mixed[idx]), mixed_rate = if (n > 0) mean(f$mixed[idx]) else NA_real_,
          mixed_boot_lo = boot_rate_ci(f$mixed[idx])[1], mixed_boot_hi = boot_rate_ci(f$mixed[idx])[2],
          n_flip_resid = sum(fr[idx]), flip_resid_rate = if (n > 0) mean(fr[idx]) else NA_real_,
          note = sprintf("r_global pooled over T cells; %d features in universe", u$n_feat))
      }
    }
  }
}
controlled <- do.call(rbind, ctrl)
utils::write.csv(controlled, file.path(out_dir, "tcell_controlled.csv"), row.names = FALSE)

# Pair-level table for design (ii), sc_110 universe, |r_T| >= 0.2 (the "41 pairs").
{
  r_g <- st_t110$r_global
  f200 <- sign_flags(r_g, st_t110$r_ctx, 0.05)
  f50 <- sign_flags(r_g, st_t110_50$r_ctx, 0.05)
  strong <- which(is.finite(r_g) & abs(r_g) >= 0.2 & r_g != 0)
  pt <- data.frame(pair_id = pairs110$pair_id[strong], feature_x = pairs110$feature_x[strong], feature_y = pairs110$feature_y[strong],
                   r_T_pooled = r_g[strong], stringsAsFactors = FALSE)
  for (lv in st_t110_50$levels_use) pt[[paste0("r_", lv)]] <- st_t110_50$r_ctx[strong, lv]
  pt$majority_nmin200 <- f200$maj[strong]; pt$eligible_nmin200 <- f200$eligible[strong]
  pt$reversal_nmin200 <- f200$reversal[strong]; pt$mixed_nmin200 <- f200$mixed[strong]
  pt$majority_nmin50 <- f50$maj[strong]; pt$eligible_nmin50 <- f50$eligible[strong]
  pt$reversal_nmin50 <- f50$reversal[strong]; pt$mixed_nmin50 <- f50$mixed[strong]
  pt$r_resid_nmin200 <- st_t110$r_resid[strong]
  utils::write.csv(pt, file.path(out_dir, "tcell_only_strong_pairs_thr0.2.csv"), row.names = FALSE)
}

# =============================================================================
# (c) Sparse-aware check of residualized flip rates
# =============================================================================
log_msg("(c) sparse-aware check")
SPARSE_COLS <- c("dataset", "context_definition", "level", "unit_size_m", "n_units", "units_per_context",
                 "n_contexts_used", "contexts_used", "method", "n_cells", "n_pairs",
                 "n_flip_resid", "flip_resid_rate", "flip_boot_lo", "flip_boot_hi",
                 "n_pairs_absr_ge_0.05", "n_flip_absr_ge_0.05", "flip_rate_absr_ge_0.05",
                 "n_pairs_absr_ge_0.2", "n_flip_absr_ge_0.2", "flip_rate_absr_ge_0.2",
                 "simpson_thr", "simpson_eps", "simpson_n_eligible", "simpson_n_reversal", "simpson_reversal_rate", "note")
sparse_row <- function(...) {
  v <- list(...)
  base <- setNames(as.list(rep(NA, length(SPARSE_COLS))), SPARSE_COLS)
  base[names(v)] <- v
  as.data.frame(base, stringsAsFactors = FALSE)
}
sparse_from_stats <- function(st, dataset, ctxdef, level, m = NA, n_units = NA, units_per_context = "", n_cells, note = "") {
  # Match the cross-setting convention: analyze only pairs with finite pooled
  # and residualized correlations and at least two valid context correlations.
  r_ctx_valid <- st$r_ctx
  r_ctx_valid[!is.finite(r_ctx_valid) | abs(r_ctx_valid) >= 1] <- NA_real_
  k_valid <- rowSums(is.finite(r_ctx_valid))
  ok <- is.finite(st$r_global) & is.finite(st$r_resid) & k_valid >= 2
  fr <- flip_resid(st$r_global, st$r_resid)
  f <- sign_flags(st$r_global, r_ctx_valid, SIMPSON_EPS)
  elig <- ok & abs(st$r_global) >= SIMPSON_THR & f$eligible
  ci <- boot_rate_ci(fr[ok])
  # Flip rates restricted to pairs whose (method-specific) |r_global| clears a threshold,
  # because the all-pairs rate is dominated by near-zero correlations whose sign is noise.
  ok05 <- ok & abs(st$r_global) >= 0.05
  ok20 <- ok & abs(st$r_global) >= 0.2
  sparse_row(dataset = dataset,
             n_pairs_absr_ge_0.05 = sum(ok05), n_flip_absr_ge_0.05 = sum(fr[ok05]),
             flip_rate_absr_ge_0.05 = if (sum(ok05) > 0) mean(fr[ok05]) else NA_real_,
             n_pairs_absr_ge_0.2 = sum(ok20), n_flip_absr_ge_0.2 = sum(fr[ok20]),
             flip_rate_absr_ge_0.2 = if (sum(ok20) > 0) mean(fr[ok20]) else NA_real_, context_definition = ctxdef, level = level, unit_size_m = m, n_units = n_units,
             units_per_context = units_per_context, n_contexts_used = length(st$levels_use),
             contexts_used = paste(sprintf("%s(%d)", st$levels_use, st$n_by_level), collapse = ";"),
             method = st$method, n_cells = n_cells, n_pairs = sum(ok),
             n_flip_resid = sum(fr[ok]), flip_resid_rate = mean(fr[ok]), flip_boot_lo = ci[1], flip_boot_hi = ci[2],
             simpson_thr = SIMPSON_THR, simpson_eps = SIMPSON_EPS, simpson_n_eligible = sum(elig),
             simpson_n_reversal = sum(f$reversal[elig]), simpson_reversal_rate = if (sum(elig) > 0) mean(f$reversal[elig]) else NA_real_,
             note = note)
}
sp <- list()
SPARSE_METHODS <- c("pearson", "spearman", "spearman_rankfirst")
run_stats <- function(X, Y, pairs, ctx, n_min, meth) {
  if (meth == "spearman_rankfirst") pair_stats_generic(X, Y, pairs, ctx, n_min, "spearman", "rank_then_residualize")
  else pair_stats_generic(X, Y, pairs, ctx, n_min, meth)
}
SP_NOTE <- c(pearson = "",
             spearman = "Spearman: ranks over all retained cells (global), within context (r_c); r_resid = Spearman of context-mean-centered values (pipeline definition)",
             spearman_rankfirst = "rank-first variant: r_resid = Pearson of context-mean-centered ranks (ranks taken over the retained cells)")

## PBMC scRNA, coarse cell types -------------------------------------------
sp[[length(sp) + 1]] <- sparse_from_stats(st_coarse, "PBMC_scRNA", "coarse_cell_types", "cell", n_cells = st_coarse$n_cells,
                                          note = "recomputed on the common retained-context population")
sp[[length(sp) + 1]] <- sparse_from_stats(st_coarse_sp, "PBMC_scRNA", "coarse_cell_types", "cell", n_cells = st_coarse_sp$n_cells,
                                          note = SP_NOTE[["spearman"]])
st_coarse_rf <- run_stats(X100, NULL, pairs100, ctx_coarse, N_MIN_MAIN, "spearman_rankfirst")
check_close(st_coarse_rf$r_global, st_coarse_sp$r_global, 1e-12, "rank-first Spearman r_global equals Spearman r_global")
check_close(st_coarse_rf$r_global[vi], ref_g, 1e-10, "retained-population rank-first r_global vs slow Spearman")
one_ctx <- ctx_coarse == st_coarse_rf$levels_use[1]
st_one <- pair_stats_generic(X100[, one_ctx, drop = FALSE], NULL, pairs100[vi, , drop = FALSE],
                             ctx_coarse[one_ctx], N_MIN_MAIN, "spearman", "rank_then_residualize")
check_close(st_one$r_global, st_one$r_resid, 1e-10, "one-context rank-first r_global equals r_resid")
check_close(st_one$r_global, st_one$r_ctx[, 1], 1e-10, "one-context rank-first r_global equals r_context")
sp[[length(sp) + 1]] <- sparse_from_stats(st_coarse_rf, "PBMC_scRNA", "coarse_cell_types", "cell", n_cells = st_coarse_rf$n_cells,
                                          note = SP_NOTE[["spearman_rankfirst"]])
for (m in PB_UNIT_SIZES) {
  pb <- make_pseudobulk(list(RNA = cnt_qc), ctx_coarse, N_MIN_MAIN, m, seed = SEED + m)
  Xpb <- log_normalize_dense_subset(pb$agg$RNA, ar100$inputs$feature_x)
  upc <- paste(sprintf("%s(%d)", names(pb$units_per_context), as.integer(pb$units_per_context)), collapse = ";")
  for (meth in SPARSE_METHODS) {
    stpb <- run_stats(Xpb, NULL, pairs100, pb$unit_context, PB_N_MIN_UNITS, meth)
    sp[[length(sp) + 1]] <- sparse_from_stats(stpb, "PBMC_scRNA", "coarse_cell_types", "pseudobulk_within_context", m = m,
                                              n_units = pb$n_units, units_per_context = upc, n_cells = pb$n_cells_used,
                                              note = paste0(sprintf("raw counts summed over %d-%d cells per unit within each cell type (n >= %d cells), then LogNormalize; contexts need >= %d units",
                                                                    m, m + 1, N_MIN_MAIN, PB_N_MIN_UNITS),
                                                            if (nzchar(SP_NOTE[[meth]])) paste0("; ", SP_NOTE[[meth]]) else ""))
  }
}
sp[[length(sp) + 1]] <- sparse_row(dataset = "PBMC_scRNA", context_definition = "coarse_cell_types", level = "pseudobulk_by_context",
                                   n_units = length(st_coarse$levels_use), n_contexts_used = length(st_coarse$levels_use),
                                   contexts_used = paste(sprintf("%s(%d)", st_coarse$levels_use, st_coarse$n_by_level), collapse = ";"),
                                   method = "pearson", n_cells = st_coarse$n_cells_resid, n_pairs = nrow(pairs100),
                                   note = "infeasible: one pseudobulk unit per cell type gives 4 units (7 if all types kept); residualizing one unit per context leaves zero variance, so r_resid is undefined and a 4-point r_global is not a per-pair statistic")

## CITE-seq, joint WNN clusters ---------------------------------------------
log_msg("loading CITE-seq and replicating the pipeline preprocessing")
cite_src <- sub("^local_rds(_list)?:", "", ar200$inputs$dataset$source)
obj2 <- readRDS(cite_src)
cnt2 <- SeuratObject::LayerData(obj2[["RNA"]], layer = "counts")
adt2 <- SeuratObject::LayerData(obj2[["ADT"]], layer = "counts")
md2 <- obj2@meta.data
stopifnot(identical(colnames(cnt2), rownames(md2)), identical(colnames(adt2), colnames(cnt2)))
qc2 <- qc_keep_cells(cnt2, md2)
ctx_joint <- as.character(ar200$context_defs$Joint_RNA_ADT)
stopifnot(sum(qc2$keep) == length(ctx_joint))
cnt2_qc <- cnt2[, qc2$keep, drop = FALSE]
adt2_qc <- as(adt2[, qc2$keep, drop = FALSE], "CsparseMatrix")
ln2 <- log_normalize_sparse(cnt2_qc)
rna_feats <- ar200$inputs$feature_sets$rna_features
adt_feats <- ar200$inputs$feature_sets$adt_features
R200 <- as.matrix(ln2[rna_feats, , drop = FALSE])
A200 <- clr_by_cell(adt2_qc)[adt_feats, , drop = FALSE]
pairs200 <- ar200$pairs
ps200 <- ar200$per_definition_stats$Joint_RNA_ADT

st_cite <- pair_stats_generic(A200, R200, pairs200, ctx_joint, N_MIN_MAIN, "pearson")
for (lv in st_cite$levels_use) check_close(st_cite$r_ctx[, lv], ps200$r_by_context[, lv], 1e-10, sprintf("CITE-seq r_c[cluster %s] vs sc_200 cache", lv))
check_close(st_cite$r_resid, ps200$summary_pairs$r_resid, 1e-10, "CITE-seq r_resid vs sc_200 cache")
st_cite_sp <- pair_stats_generic(A200, R200, pairs200, ctx_joint, N_MIN_MAIN, "spearman")
set.seed(SEED)
vi2 <- sample.int(nrow(pairs200), 200)
keep200 <- ctx_joint %in% st_cite_sp$levels_use
ref_g2 <- vapply(vi2, function(i) suppressWarnings(stats::cor(A200[pairs200$idx_x[i], keep200], R200[pairs200$idx_y[i], keep200], method = "spearman")), numeric(1))
check_close(st_cite_sp$r_global[vi2], ref_g2, 1e-10, "CITE-seq retained-population Spearman r_global vs cor(method = 'spearman')")
ctx_ok2 <- ctx_joint
ctx_ok2[!(ctx_ok2 %in% st_cite_sp$levels_use)] <- NA
ref_r2 <- vapply(vi2, function(i) mean_residualized_cor(A200[pairs200$idx_x[i], ], R200[pairs200$idx_y[i], ], ctx_ok2, method = "spearman"), numeric(1))
check_close(st_cite_sp$r_resid[vi2], ref_r2, 1e-10, "CITE-seq Spearman r_resid vs mean_residualized_cor(method = 'spearman')")

sp[[length(sp) + 1]] <- sparse_from_stats(st_cite, "CITEseq_protein_RNA", "Joint_RNA_ADT_WNN", "cell", n_cells = st_cite$n_cells,
                                          note = "recomputed on the common retained-context population")
sp[[length(sp) + 1]] <- sparse_from_stats(st_cite_sp, "CITEseq_protein_RNA", "Joint_RNA_ADT_WNN", "cell", n_cells = st_cite_sp$n_cells,
                                          note = SP_NOTE[["spearman"]])
st_cite_rf <- run_stats(A200, R200, pairs200, ctx_joint, N_MIN_MAIN, "spearman_rankfirst")
check_close(st_cite_rf$r_global, st_cite_sp$r_global, 1e-12, "CITE-seq rank-first Spearman r_global equals Spearman r_global")
check_close(st_cite_rf$r_global[vi2], ref_g2, 1e-10, "CITE-seq retained-population rank-first r_global vs slow Spearman")
sp[[length(sp) + 1]] <- sparse_from_stats(st_cite_rf, "CITEseq_protein_RNA", "Joint_RNA_ADT_WNN", "cell", n_cells = st_cite_rf$n_cells,
                                          note = SP_NOTE[["spearman_rankfirst"]])
for (m in PB_UNIT_SIZES) {
  pb <- make_pseudobulk(list(RNA = cnt2_qc, ADT = adt2_qc), ctx_joint, N_MIN_MAIN, m, seed = SEED + m)
  Rpb <- log_normalize_dense_subset(pb$agg$RNA, rna_feats)
  Apb <- clr_by_cell(pb$agg$ADT)[adt_feats, , drop = FALSE]
  upc <- paste(sprintf("%s(%d)", names(pb$units_per_context), as.integer(pb$units_per_context)), collapse = ";")
  for (meth in SPARSE_METHODS) {
    stpb <- run_stats(Apb, Rpb, pairs200, pb$unit_context, PB_N_MIN_UNITS, meth)
    sp[[length(sp) + 1]] <- sparse_from_stats(stpb, "CITEseq_protein_RNA", "Joint_RNA_ADT_WNN", "pseudobulk_within_context", m = m,
                                              n_units = pb$n_units, units_per_context = upc, n_cells = pb$n_cells_used,
                                              note = paste0(sprintf("RNA and ADT raw counts summed over %d-%d cells per unit within each WNN cluster (n >= %d cells); RNA LogNormalize and ADT CLR on the sums; contexts need >= %d units",
                                                                    m, m + 1, N_MIN_MAIN, PB_N_MIN_UNITS),
                                                            if (nzchar(SP_NOTE[[meth]])) paste0("; ", SP_NOTE[[meth]]) else ""))
  }
}
sp[[length(sp) + 1]] <- sparse_row(dataset = "CITEseq_protein_RNA", context_definition = "Joint_RNA_ADT_WNN", level = "pseudobulk_by_context",
                                   n_units = length(st_cite$levels_use), n_contexts_used = length(st_cite$levels_use),
                                   contexts_used = paste(sprintf("%s(%d)", st_cite$levels_use, st_cite$n_by_level), collapse = ";"),
                                   method = "pearson", n_cells = st_cite$n_cells_resid, n_pairs = nrow(pairs200),
                                   note = "infeasible: one pseudobulk unit per WNN cluster gives 11 units (20 if all clusters kept); residualizing one unit per context leaves zero variance, so r_resid is undefined")
sparse <- do.call(rbind, sp)
utils::write.csv(sparse, file.path(out_dir, "sparse_check.csv"), row.names = FALSE)

# =============================================================================
# (d) PBMC 3k donor / metadata check
# =============================================================================
log_msg("(d) donor / metadata check")
suffixes <- table(sub("^.*-", "", colnames(cnt)))
donor <- data.frame(
  item = c("source_file_used_by_loader", "loader_source_tag_in_cache", "n_cells_raw", "n_cells_after_qc",
           "metadata_columns", "orig_ident_levels", "project_name", "barcode_suffixes",
           "donor_or_sample_column_present", "raw_input_tarball", "dataset_provenance"),
  value = c(pbmc_src, ar100$inputs$extra$dataset$source, ncol(cnt), sum(qc$keep),
            paste(colnames(md), collapse = ";"), paste(unique(as.character(md$orig.ident)), collapse = ";"),
            obj@project.name, paste(sprintf("-%s(%d)", names(suffixes), as.integer(suffixes)), collapse = ";"),
            "no", file.path(dirname(pbmc_src), "pbmc3k.tar.gz"),
            "10x Genomics pbmc3k (filtered_gene_bc_matrices/hg19); one Cell Ranger sample, one GEM well (all barcodes -1); no donor field in metadata"),
  stringsAsFactors = FALSE
)
utils::write.csv(donor, file.path(out_dir, "pbmc_donor_check.csv"), row.names = FALSE)

# =============================================================================
# wrap up
# =============================================================================
writeLines(utils::capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
writeLines(sprintf("runtime_minutes,%.2f", elapsed), con = file.path(out_dir, "runtime.txt"))

cat("\n==== key numbers ====\n")
key <- controlled[controlled$eps == 0.05 & controlled$n_min == N_MIN_MAIN,
                  c("design", "pair_universe", "partition", "threshold", "n_pairs_strong", "n_eligible_common",
                    "n_reversal", "reversal_rate", "reversal_boot_lo", "reversal_boot_hi", "reversal_exact_hi", "n_mixed", "mixed_rate")]
print(key, row.names = FALSE, digits = 4)
print(sparse[, c("dataset", "level", "unit_size_m", "n_units", "method", "n_pairs", "n_flip_resid", "flip_resid_rate",
                 "n_pairs_absr_ge_0.05", "flip_rate_absr_ge_0.05", "n_pairs_absr_ge_0.2", "flip_rate_absr_ge_0.2",
                 "simpson_n_eligible", "simpson_n_reversal", "simpson_reversal_rate")], row.names = FALSE, digits = 4)
print(donor, row.names = FALSE)
log_msg("done in %.1f minutes", elapsed)
