#!/usr/bin/env Rscript

# Multi-donor analysis of single-cell correlation reversals, with alternative
# context definitions and donor-by-cell-type pseudobulk summaries.
#
# Data: Kang et al. 2018 (Nat Biotechnol 36:89-94), GEO GSE96583, batch 2:
#   GSM2560248 (control PBMCs, 8 lupus donors, multiplexed) and
#   GSM2560249 (the same 8 donors after 6 h IFN-beta stimulation).
#   Donor assignments (demuxlet), cell-type labels, and singlet/doublet calls
#   come from the GEO supplementary file GSE96583_batch2.total.tsne.df.tsv.gz.
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583
#
# Run from the repo root:  Rscript scripts/revision/multidonor.R
# Quick smoke test (few pairs, few bootstrap replicates, scratch output dir):
#   MULTIDONOR_QUICK=1 MULTIDONOR_OUT=/some/dir Rscript scripts/revision/multidonor.R
#
# Conventions copied from the paper's PBMC 3k pipeline
# (scripts/sc_100_pbmc_coarse_celltypes.R, R/sc_context_pipeline.R):
#   - Seurat v5: min 200 detected features, max 20% mitochondrial reads,
#     LogNormalize, 2,000 HVGs (vst)
#   - features: marker panel + HVGs detected in >= 1% of cells, capped at 1,800
#   - 200,000 random feature pairs (sc_make_feature_pairs, seed 100)
#   - a context enters the within-context statistics only when it holds at
#     least n_min cells (n_min = 200 in the paper; 50 is run as sensitivity)
#   - sign tolerance eps = 0.05; mixed-sign = both signs among contexts;
#     Simpson reversal = sign(r_global) != sign(context majority), ties and
#     r_global = 0 excluded, reported among |r_global| >= 0.2;
#     residualized flip = sign(r_resid) != sign(r_global), where r_resid is the
#     correlation after subtracting per-context means from both genes
#   - Q with weights n_c - 3 on Fisher z; I^2 = max(0, (Q - (k - 1)) / Q)
#
# Outputs (results/revision/multidonor/):
#   dataset_summary.csv          what was used (accession, counts, parameters)
#   cell_counts.csv              cells per donor x condition x cell type (post QC)
#   features.csv                 the genes entering the pair sample
#   cell_level_rates.csv         one row per (context definition, n_min), with CIs
#   cell_level_reversal_curve.csv Simpson reversal vs |r_global| threshold per definition
#   per_donor_rates.csv          single-donor analyses (donor = pooled population)
#   donor_replication.csv        do the all-cells flips replicate donor by donor
#   pseudobulk_units.csv         the pseudobulk units and their cell counts
#   pseudobulk_rates.csv         pseudobulk rates per (unit set, aggregation, context),
#                                with a within-context permutation null
#   pseudobulk_concordance.csv   cell-level vs pseudobulk agreement on the same pairs
#   pbmc3k_reference.csv         the paper's PBMC 3k rates recomputed with this code
#   Figure_MultiDonor.pdf
#   runtime.txt, sessionInfo.txt
#   cache/                       downloaded GEO files, processed Seurat object,
#                                pair statistics per run

t_start <- Sys.time()
section_times <- list()
tic <- function() Sys.time()
toc <- function(label, t0) {
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  section_times[[label]] <<- dt
  message(sprintf("[%s] %.1f s", label, dt))
}

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

source("00-paths.R")
source(file.path("R", "contextcorr.R"))
source(file.path("R", "sc_datasets.R"))
source(file.path("R", "sc_context_pipeline.R"))

QUICK <- identical(Sys.getenv("MULTIDONOR_QUICK"), "1")
FORCE <- identical(Sys.getenv("MULTIDONOR_FORCE"), "1")
out_dir <- Sys.getenv("MULTIDONOR_OUT", unset = file.path("results", "revision", "multidonor"))
cache_dir <- file.path(out_dir, "cache")
geo_dir <- file.path(cache_dir, "geo")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

n_cores <- max(1L, min(24L, parallel::detectCores() - 4L))

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
accession <- "GSE96583"
geo_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583"
geo_suppl_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/"

min_features <- 200
max_mt <- 20
n_hvg <- 2000
min_detect_frac <- 0.01
max_features <- 1800
n_pairs <- if (QUICK) 5000L else 200000L
seed_pairs <- 100                 # same seed as scripts/sc_100_pbmc_coarse_celltypes.R
eps <- 0.05
t_main <- 0.2
n_min_paper <- 200
n_min_grid <- c(200, 50)
B_boot <- if (QUICK) 50L else 500L
seed_boot <- 20260905
min_cells_per_unit <- 20          # pseudobulk units with fewer cells are dropped
n_min_units <- 4                  # within-context pseudobulk r needs >= 4 units
n_perm <- if (QUICK) 4L else 20L  # within-context permutations for the pseudobulk null
seed_perm <- 20260907
thresholds <- seq(0, 0.5, by = 0.05)
eps_grid <- c(0, 0.05, 0.1)
n_plot_points <- 50000L
seed_plot <- 20260906

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Two-sided p-value of a Pearson r on df residual degrees of freedom.
p_from_r <- function(r, df) {
  p <- rep(NA_real_, length(r))
  ok <- is.finite(r) & abs(r) < 1 & df >= 1
  tt <- r[ok] * sqrt(df / (1 - r[ok]^2))
  p[ok] <- 2 * stats::pt(-abs(tt), df = df)
  p[is.finite(r) & abs(r) >= 1 & df >= 1] <- 0
  p
}

# Rates for one run from the summary_pairs data frame produced by
# sc_compute_pair_statistics_fast_pearson, with pair-level bootstrap CIs.
# n_obs: retained observations (cells or units) behind r_global; n_obs_resid
# and k_retained: the same observations and contexts behind r_resid, giving its residual
# df = n_obs_resid - k_retained - 1. A "confident flip" is a sign flip between
# an r_global and an r_resid that are each individually significant at
# alpha = 0.05; it matters for pseudobulk, where r_resid rests on few units.
rates_from_summary <- function(sp, n_obs, n_obs_resid, k_retained, thr = t_main, B = B_boot, seed = seed_boot, alpha = 0.05) {
  rg <- sp$r_global
  rr <- sp$r_resid
  flip_elig <- is.finite(rg) & is.finite(rr) & rg != 0 & rr != 0
  flip_ind <- flip_elig & (sign(rg) != sign(rr))
  flip02_elig <- flip_elig & abs(rg) >= thr
  flip02_ind <- flip02_elig & flip_ind
  simp_elig <- is.finite(rg) & rg != 0 & abs(rg) >= thr & sp$majority_sign != 0
  simp_ind <- simp_elig & sp$flag_simpson
  mixed_elig <- sp$k_contexts_used >= 2
  mixed_ind <- mixed_elig & sp$n_pos > 0 & sp$n_neg > 0
  I2 <- sp$I2
  I2_ok <- is.finite(I2)
  df_g <- n_obs - 2
  df_r <- n_obs_resid - k_retained - 1
  p_g <- p_from_r(rg, df_g)
  p_r <- p_from_r(rr, df_r)
  conf_elig <- flip_elig & is.finite(p_g) & is.finite(p_r) & p_g < alpha & p_r < alpha
  conf_ind <- conf_elig & flip_ind

  rate <- function(ind, elig) if (sum(elig) > 0) sum(ind) / sum(elig) else NA_real_
  point <- c(
    flip_rate = rate(flip_ind, flip_elig),
    flip_rate_abs02 = rate(flip02_ind, flip02_elig),
    simpson_rate_abs02 = rate(simp_ind, simp_elig),
    mixed_sign_rate = rate(mixed_ind, mixed_elig),
    median_I2 = if (any(I2_ok)) stats::median(I2[I2_ok]) else NA_real_,
    confident_flip_rate_cond = rate(conf_ind, conf_elig),
    confident_flip_rate_all = rate(conf_ind, flip_elig)
  )

  ind_list <- list(flip_ind, flip02_ind, simp_ind, mixed_ind, conf_ind, conf_ind)
  elig_list <- list(flip_elig, flip02_elig, simp_elig, mixed_elig, conf_elig, flip_elig)
  n <- nrow(sp)
  ci <- matrix(NA_real_, nrow = 2, ncol = 7)
  if (B > 0 && n > 1) {
    set.seed(seed)
    bm <- matrix(NA_real_, nrow = B, ncol = 7)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      for (j in seq_along(ind_list)) {
        e <- sum(elig_list[[j]][idx])
        bm[b, j] <- if (e > 0) sum(ind_list[[j]][idx]) / e else NA_real_
      }
      v <- I2[idx]
      v <- v[is.finite(v)]
      bm[b, 7] <- if (length(v) > 0) stats::median(v) else NA_real_
    }
    ci <- apply(bm, 2, function(col) {
      if (all(is.na(col))) c(NA_real_, NA_real_) else stats::quantile(col, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    })
  }

  data.frame(
    n_pairs = n,
    n_pairs_flip_eval = sum(flip_elig),
    flip_rate = point[["flip_rate"]], flip_lo = ci[1, 1], flip_hi = ci[2, 1],
    n_pairs_flip_eval_abs02 = sum(flip02_elig),
    flip_rate_abs02 = point[["flip_rate_abs02"]], flip_abs02_lo = ci[1, 2], flip_abs02_hi = ci[2, 2],
    n_pairs_simpson_elig_abs02 = sum(simp_elig),
    simpson_rate_abs02 = point[["simpson_rate_abs02"]], simpson_lo = ci[1, 3], simpson_hi = ci[2, 3],
    n_pairs_mixed_eval = sum(mixed_elig),
    mixed_sign_rate = point[["mixed_sign_rate"]], mixed_lo = ci[1, 4], mixed_hi = ci[2, 4],
    n_pairs_I2 = sum(I2_ok),
    median_I2 = point[["median_I2"]], I2_lo = ci[1, 7], I2_hi = ci[2, 7],
    df_global = df_g, df_resid = df_r,
    n_pairs_both_significant = sum(conf_elig),
    n_confident_flips = sum(conf_ind),
    confident_flip_rate_cond = point[["confident_flip_rate_cond"]], confident_cond_lo = ci[1, 5], confident_cond_hi = ci[2, 5],
    confident_flip_rate_all = point[["confident_flip_rate_all"]], confident_all_lo = ci[1, 6], confident_all_hi = ci[2, 6],
    eps = eps, abs_r_threshold = thr, alpha = alpha, B_boot = B, seed_boot = seed,
    stringsAsFactors = FALSE
  )
}

# Context-size bookkeeping for one run.
context_summary <- function(context, n_min) {
  tab <- table(context)
  keep <- tab >= n_min
  data.frame(
    n_cells = length(context),
    k_contexts_total = length(tab),
    k_contexts_retained = sum(keep),
    n_cells_in_retained_contexts = sum(tab[keep]),
    min_context_size = as.integer(min(tab)),
    median_context_size = as.numeric(stats::median(as.integer(tab))),
    max_context_size = as.integer(max(tab)),
    contexts_retained = paste(names(tab)[keep], collapse = ";"),
    stringsAsFactors = FALSE
  )
}

# Q and I2 recomputed from the per-context matrices with exactly the validity
# rule of heterogeneity_Q_I2 in R/contextcorr.R (finite r_c, n_c >= 4, |r_c| < 1).
# The fast path clamps |r_c| = 1 instead of excluding it; that never matters
# with hundreds of cells per context but can with 8 pseudobulk units.
recompute_heterogeneity <- function(st) {
  R <- st$r_by_context
  N <- st$n_by_context
  valid <- is.finite(R) & is.finite(N) & N >= 4 & abs(R) < 1
  W <- N - 3
  W[!valid] <- 0
  Z <- contextcorr_fisher_z(R)
  Z[!valid] <- 0
  sum_w <- rowSums(W)
  sum_wz <- rowSums(W * Z)
  sum_wz2 <- rowSums(W * Z^2)
  k <- rowSums(valid)
  df <- pmax(k - 1L, 0L)
  Q <- rep(NA_real_, nrow(R))
  I2 <- Q
  p <- Q
  ok <- k >= 2 & sum_w > 0
  Q[ok] <- sum_wz2[ok] - sum_wz[ok]^2 / sum_w[ok]
  p[ok] <- stats::pchisq(Q[ok], df = df[ok], lower.tail = FALSE)
  I2[ok] <- ifelse(is.finite(Q[ok]) & Q[ok] > 0, pmax(0, (Q[ok] - df[ok]) / Q[ok]), NA_real_)
  st$summary_pairs$Q <- Q
  st$summary_pairs$p_heterogeneity <- p
  st$summary_pairs$I2 <- I2
  st
}

# Check the vectorized results against the paper's per-pair helpers
# (R/contextcorr.R) on n_check random pairs. This replaces the library's
# built-in check, which compares Q before the |r_c| < 1 exclusion and therefore
# fails whenever a sparse gene pair is perfectly correlated inside a context.
validate_run <- function(st, expr, pairs, context, n_min, n_check = 25, tol = 1e-6) {
  n <- nrow(pairs)
  if (n == 0 || n_check == 0) return(invisible(TRUE))
  idx <- seq(1L, n, length.out = min(n_check, n))
  idx <- unique(as.integer(round(idx)))
  context <- as.factor(context)
  lev <- st$context_levels
  same <- function(a, b, tol) (!is.finite(a) & !is.finite(b)) | (is.finite(a) & is.finite(b) & abs(a - b) <= tol)
  for (i in idx) {
    x <- as.numeric(expr[pairs$idx_x[i], ])
    y <- as.numeric(expr[pairs$idx_y[i], ])
    cc <- context_correlations(x, y, context = context, method = "pearson")
    r_use <- cc$r_by_context
    r_use[cc$n_by_context < n_min] <- NA_real_
    sp <- st$summary_pairs
    retained <- is.finite(x) & is.finite(y) & !is.na(context) &
      as.character(context) %in% names(cc$n_by_context)[cc$n_by_context >= n_min]
    r_global_use <- if (sum(retained) >= 2L) {
      suppressWarnings(stats::cor(x[retained], y[retained], method = "pearson"))
    } else {
      NA_real_
    }
    if (!same(r_global_use, sp$r_global[i], tol)) stop("validation: r_global mismatch at pair ", i)
    if (!all(same(as.numeric(r_use[lev]), as.numeric(st$r_by_context[i, ]), tol))) stop("validation: r_by_context mismatch at pair ", i)
    sf <- simpson_flag(r_global_use, r_use, eps = eps, majority = "count")
    if (!identical(isTRUE(sf$flag), isTRUE(sp$flag_simpson[i]))) stop("validation: simpson flag mismatch at pair ", i)
    ht <- heterogeneity_Q_I2(r_use, cc$n_by_context)
    if (!same(ht$Q, sp$Q[i], 1e-5)) stop("validation: Q mismatch at pair ", i)
    if (!same(ht$I2, sp$I2[i], 1e-5)) stop("validation: I2 mismatch at pair ", i)
    ctx_ok <- as.character(context)
    ctx_ok[!(ctx_ok %in% names(cc$n_by_context)[cc$n_by_context >= n_min])] <- NA_character_
    r_res <- if (any(!is.na(ctx_ok))) suppressWarnings(mean_residualized_cor(x, y, context = ctx_ok, method = "pearson")) else NA_real_
    if (!same(r_res, sp$r_resid[i], 1e-5)) stop("validation: r_resid mismatch at pair ", i)
  }
  invisible(TRUE)
}

run_pair_stats <- function(expr, pairs, context, n_min, cache_file, force = FALSE, validate_n = 25) {
  st <- sc_compute_pair_statistics_fast_pearson(
    expr_x = expr, expr_y = expr, pairs = pairs, context = context,
    eps = eps, majority = "count", n_min = n_min,
    cache_file = cache_file, force = force, progress_every = 1e9, validate_n = 0
  )
  st <- recompute_heterogeneity(st)
  validate_run(st, expr, pairs, context, n_min, n_check = validate_n)
  st
}

fmt_pct <- function(x, d = 1) ifelse(is.finite(x), sprintf(paste0("%.", d, "f%%"), 100 * x), "NA")
fmt_num <- function(x, d = 3) ifelse(is.finite(x), formatC(x, format = "f", digits = d, big.mark = ","), "NA")
fmt_int <- function(x) format(x, big.mark = ",", trim = TRUE)

# ---------------------------------------------------------------------------
# 1. Download (only if missing) and load GSE96583 batch 2
# ---------------------------------------------------------------------------
t0 <- tic()
suppl_files <- c("GSE96583_RAW.tar", "GSE96583_batch2.genes.tsv.gz", "GSE96583_batch2.total.tsne.df.tsv.gz")
for (f in suppl_files) {
  dest <- file.path(geo_dir, f)
  if (!file.exists(dest)) {
    message("Downloading ", f)
    ok <- tryCatch(utils::download.file(paste0(geo_suppl_url, f), dest, mode = "wb", quiet = TRUE), error = function(e) 1L)
    if (!identical(as.integer(ok), 0L) || !file.exists(dest)) stop("Download failed for ", f, " from ", geo_suppl_url)
  }
}
raw_members <- c("GSM2560248_2.1.mtx.gz", "GSM2560248_barcodes.tsv.gz",
                 "GSM2560249_2.2.mtx.gz", "GSM2560249_barcodes.tsv.gz")
if (!all(file.exists(file.path(geo_dir, raw_members)))) {
  utils::untar(file.path(geo_dir, "GSE96583_RAW.tar"), files = raw_members, exdir = geo_dir)
}
stopifnot(all(file.exists(file.path(geo_dir, raw_members))))

seurat_cache <- file.path(cache_dir, "seurat_kang2018_batch2_processed.rds")
qc_cache <- file.path(cache_dir, "qc_counts.rds")

if (file.exists(seurat_cache) && file.exists(qc_cache)) {
  seu <- readRDS(seurat_cache)
  qc <- readRDS(qc_cache)
} else {
  genes_file <- file.path(geo_dir, "GSE96583_batch2.genes.tsv.gz")
  m_ctrl <- Seurat::ReadMtx(
    mtx = file.path(geo_dir, "GSM2560248_2.1.mtx.gz"),
    cells = file.path(geo_dir, "GSM2560248_barcodes.tsv.gz"),
    features = genes_file, feature.column = 2, cell.column = 1,
    unique.features = TRUE, strip.suffix = FALSE
  )
  m_stim <- Seurat::ReadMtx(
    mtx = file.path(geo_dir, "GSM2560249_2.2.mtx.gz"),
    cells = file.path(geo_dir, "GSM2560249_barcodes.tsv.gz"),
    features = genes_file, feature.column = 2, cell.column = 1,
    unique.features = TRUE, strip.suffix = FALSE
  )
  stopifnot(identical(rownames(m_ctrl), rownames(m_stim)))
  bc_ctrl <- colnames(m_ctrl)
  bc_stim <- colnames(m_stim)

  # Per-cell metadata (donor = demuxlet "ind", condition = "stim", cell type =
  # "cell", singlet/doublet call = "multiplets"). Rows are the control barcodes
  # in file order followed by the stimulated barcodes in file order; the 313
  # stimulated barcodes that collide with a control barcode carry an extra "1"
  # ("-11" instead of "-1") in this file. Both facts are verified here.
  md <- utils::read.delim(gzfile(file.path(geo_dir, "GSE96583_batch2.total.tsne.df.tsv.gz")),
                          header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(nrow(md) == length(bc_ctrl) + length(bc_stim))
  rn <- rownames(md)
  stopifnot(identical(rn[seq_along(bc_ctrl)], bc_ctrl))
  rs <- rn[length(bc_ctrl) + seq_along(bc_stim)]
  stopifnot(identical(sub("-11$", "-1", rs), bc_stim))
  stopifnot(all(md$stim[seq_along(bc_ctrl)] == "ctrl"), all(md$stim[length(bc_ctrl) + seq_along(bc_stim)] == "stim"))

  colnames(m_ctrl) <- paste0("ctrl_", bc_ctrl)
  colnames(m_stim) <- paste0("stim_", bc_stim)
  counts <- cbind(m_ctrl, m_stim)
  rm(m_ctrl, m_stim)

  meta <- data.frame(
    donor = paste0("D", md$ind),
    condition = md$stim,
    cell_type = md$cell,
    multiplets = md$multiplets,
    kang_cluster = md$cluster,
    gsm = c(rep("GSM2560248", length(bc_ctrl)), rep("GSM2560249", length(bc_stim))),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  qc <- list(
    n_cells_raw = ncol(counts),
    n_cells_raw_ctrl = length(bc_ctrl),
    n_cells_raw_stim = length(bc_stim),
    n_genes_raw = nrow(counts),
    n_singlets = sum(meta$multiplets == "singlet"),
    n_doublets = sum(meta$multiplets == "doublet"),
    n_ambiguous = sum(meta$multiplets == "ambs"),
    n_singlets_missing_celltype = sum(meta$multiplets == "singlet" & is.na(meta$cell_type))
  )

  keep <- meta$multiplets == "singlet" & !is.na(meta$cell_type)
  counts <- counts[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  qc$n_cells_singlet_labeled <- ncol(counts)

  seu <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta, project = "kang2018_batch2", min.cells = 0, min.features = 0)
  rm(counts)
  invisible(gc())

  seu <- sc_preprocess(
    seu, assay = "RNA", min_features = min_features, min_counts = NULL, max_mt = max_mt,
    n_hvg = n_hvg, normalize_mode = "log", residualize = FALSE, vars_to_regress = NULL,
    seed = 100, verbose = FALSE
  )
  qc$n_cells_after_qc <- ncol(seu)
  qc$n_removed_by_qc <- qc$n_cells_singlet_labeled - ncol(seu)
  saveRDS(seu, seurat_cache)
  saveRDS(qc, qc_cache)
}
md <- seu@meta.data
md$donor <- as.character(md$donor)
md$condition <- as.character(md$condition)
md$cell_type <- as.character(md$cell_type)

# Paper-style coarse labels (T, B, NK, Mono, DC, Platelet) from the provided
# ones. Explicit lookup: sc_map_coarse_label's regexes were written for the
# marker-based labels and miss the plural "B cells".
coarse_lookup <- c(
  "B cells" = "B", "CD14+ Monocytes" = "Mono", "FCGR3A+ Monocytes" = "Mono",
  "CD4 T cells" = "T", "CD8 T cells" = "T", "NK cells" = "NK",
  "Dendritic cells" = "DC", "Megakaryocytes" = "Platelet"
)
stopifnot(all(md$cell_type %in% names(coarse_lookup)))
md$cell_type_coarse <- unname(coarse_lookup[md$cell_type])
seu$cell_type_coarse <- md$cell_type_coarse
donors <- sort(unique(md$donor))
conditions <- sort(unique(md$condition))
cell_types <- sort(unique(md$cell_type))
coarse_map <- unique(md[, c("cell_type", "cell_type_coarse")])
coarse_map <- coarse_map[order(coarse_map$cell_type), ]
stopifnot(all(md$multiplets == "singlet"))
toc("load_data", t0)

message(sprintf("Cells after QC: %d; donors: %d; conditions: %s; cell types: %d",
                ncol(seu), length(donors), paste(conditions, collapse = "/"), length(cell_types)))

# ---------------------------------------------------------------------------
# 2. Features and pairs (same selection as the PBMC 3k analysis)
# ---------------------------------------------------------------------------
t0 <- tic()
feat_info <- sc_select_rna_features(
  seu, assay = "RNA", n_hvg = n_hvg, marker_genes = sc_default_marker_panel(),
  min_detect_frac = min_detect_frac, max_features = max_features
)
seu <- feat_info$seurat
features <- feat_info$features
expr <- sc_get_expression_matrix(seu, features = features, assay = "RNA", slot = "data")
features <- rownames(expr)
stopifnot(ncol(expr) == nrow(md), identical(colnames(expr), rownames(md)))

pair_file <- file.path(cache_dir, sprintf("pairs_%d.rds", n_pairs))
if (file.exists(pair_file)) {
  pairs <- readRDS(pair_file)
  stopifnot(identical(pairs$feature_x, features[pairs$idx_x]), identical(pairs$feature_y, features[pairs$idx_y]))
} else {
  pairs <- sc_make_feature_pairs(features, features, n_pairs = n_pairs, seed = seed_pairs)
  saveRDS(pairs, pair_file)
}
utils::write.csv(
  data.frame(feature = features,
             is_marker = features %in% feat_info$marker_genes_used,
             is_hvg = features %in% feat_info$hvg,
             detect_frac = as.numeric(feat_info$detect_frac[features]),
             stringsAsFactors = FALSE),
  file.path(out_dir, "features.csv"), row.names = FALSE
)
toc("features_pairs", t0)
message(sprintf("Features: %d (%d markers, %d HVGs); pairs: %d", length(features),
                length(feat_info$marker_genes_used), sum(features %in% feat_info$hvg), nrow(pairs)))

# ---------------------------------------------------------------------------
# 3. Cell-level analysis under each context definition
# ---------------------------------------------------------------------------
t0 <- tic()
ctx_defs <- list(
  cell_type = md$cell_type,
  cell_type_coarse = md$cell_type_coarse,
  donor = md$donor,
  condition = md$condition,
  donor_x_cell_type = paste(md$donor, md$cell_type, sep = "|"),
  condition_x_cell_type = paste(md$condition, md$cell_type, sep = "|"),
  donor_x_condition = paste(md$donor, md$condition, sep = "|"),
  donor_x_condition_x_cell_type = paste(md$donor, md$condition, md$cell_type, sep = "|")
)
cell_grid <- expand.grid(def = names(ctx_defs), n_min = n_min_grid, stringsAsFactors = FALSE)

cell_runs <- parallel::mclapply(seq_len(nrow(cell_grid)), function(i) {
  def <- cell_grid$def[i]
  nm <- cell_grid$n_min[i]
  ctx <- ctx_defs[[def]]
  st <- run_pair_stats(expr, pairs, ctx, n_min = nm,
                       cache_file = file.path(cache_dir, sprintf("cell_%s_nmin%d_pairs%d.rds", def, nm, n_pairs)),
                       force = FORCE)
  cs <- context_summary(ctx, nm)
  rates <- rates_from_summary(st$summary_pairs, n_obs = cs$n_cells_in_retained_contexts, n_obs_resid = cs$n_cells_in_retained_contexts,
                              k_retained = cs$k_contexts_retained)
  curve <- sc_reversal_rate_curve(st$summary_pairs$r_global, st$r_by_context, thresholds = thresholds, eps_grid = eps_grid)
  curve <- cbind(data.frame(context_definition = def, n_min = nm, stringsAsFactors = FALSE), curve)
  list(
    row = cbind(data.frame(context_definition = def, n_min = nm, stringsAsFactors = FALSE), cs, rates),
    curve = curve,
    summary = st$summary_pairs
  )
}, mc.cores = min(n_cores, nrow(cell_grid)), mc.preschedule = FALSE)

for (r in cell_runs) if (inherits(r, "try-error") || is.null(r$row)) stop("A cell-level run failed: ", paste(r, collapse = " "))
cell_level <- do.call(rbind, lapply(cell_runs, `[[`, "row"))
cell_level$row_id <- seq_len(nrow(cell_level))
cell_level <- cell_level[, c("row_id", setdiff(colnames(cell_level), "row_id"))]
utils::write.csv(cell_level, file.path(out_dir, "cell_level_rates.csv"), row.names = FALSE)
cell_curve <- do.call(rbind, lapply(cell_runs, `[[`, "curve"))
utils::write.csv(cell_curve, file.path(out_dir, "cell_level_reversal_curve.csv"), row.names = FALSE)
cell_summaries <- setNames(lapply(cell_runs, `[[`, "summary"), paste(cell_grid$def, cell_grid$n_min, sep = "_nmin"))
rm(cell_runs)
invisible(gc())
toc("cell_level", t0)

# ---------------------------------------------------------------------------
# 4. Single-donor analyses (donor = pooled population, cell type = context)
# ---------------------------------------------------------------------------
t0 <- tic()
donor_grid <- expand.grid(
  donor = donors, condition_subset = c("ctrl", "stim", "both"),
  celltype_definition = c("cell_type", "cell_type_coarse"), n_min = n_min_grid,
  stringsAsFactors = FALSE
)
donor_runs <- parallel::mclapply(seq_len(nrow(donor_grid)), function(i) {
  d <- donor_grid$donor[i]
  cs <- donor_grid$condition_subset[i]
  cd <- donor_grid$celltype_definition[i]
  nm <- donor_grid$n_min[i]
  cells <- which(md$donor == d & (cs == "both" | md$condition == cs))
  ctx <- md[[cd]][cells]
  st <- run_pair_stats(expr[, cells, drop = FALSE], pairs, ctx, n_min = nm,
                       cache_file = file.path(cache_dir, sprintf("donor_%s_%s_%s_nmin%d_pairs%d.rds", d, cs, cd, nm, n_pairs)),
                       force = FORCE)
  csum <- context_summary(ctx, nm)
  cbind(data.frame(donor = d, condition_subset = cs, celltype_definition = cd, n_min = nm, stringsAsFactors = FALSE),
        csum, rates_from_summary(st$summary_pairs, n_obs = csum$n_cells_in_retained_contexts, n_obs_resid = csum$n_cells_in_retained_contexts,
                                 k_retained = csum$k_contexts_retained))
}, mc.cores = min(n_cores, nrow(donor_grid)), mc.preschedule = FALSE)
for (r in donor_runs) if (inherits(r, "try-error") || !is.data.frame(r)) stop("A per-donor run failed: ", paste(r, collapse = " "))
per_donor <- do.call(rbind, donor_runs)
per_donor$row_id <- seq_len(nrow(per_donor))
per_donor <- per_donor[, c("row_id", setdiff(colnames(per_donor), "row_id"))]
utils::write.csv(per_donor, file.path(out_dir, "per_donor_rates.csv"), row.names = FALSE)
rm(donor_runs)
toc("per_donor", t0)

# ---------------------------------------------------------------------------
# 4b. Do the cell-level flips replicate across donors?
# For every pair, the all-cells analysis (cell type as context) gives a pooled
# and a residualized sign. Each donor analyzed alone gives its own residualized
# sign (within cell type, within donor). If flips were artifacts of one donor's
# cell-level dependence, the residualized sign would not agree across donors.
# ---------------------------------------------------------------------------
t0 <- tic()
sp_full <- cell_summaries[[sprintf("cell_type_nmin%d", n_min_paper)]]
row_full <- cell_level[cell_level$context_definition == "cell_type" & cell_level$n_min == n_min_paper, ]
p_g_full <- p_from_r(sp_full$r_global, row_full$df_global)
p_r_full <- p_from_r(sp_full$r_resid, row_full$df_resid)
ok_full <- is.finite(sp_full$r_global) & is.finite(sp_full$r_resid) & sp_full$r_global != 0 & sp_full$r_resid != 0
flip_full <- ok_full & sign(sp_full$r_global) != sign(sp_full$r_resid)
both_sig_full <- ok_full & is.finite(p_g_full) & is.finite(p_r_full) & p_g_full < 0.05 & p_r_full < 0.05
pair_groups <- list(
  all_flips = flip_full,
  confident_flips = flip_full & both_sig_full,
  confident_nonflips = ok_full & !flip_full & both_sig_full,
  all_evaluable = ok_full
)
rep_rows <- list()
for (cs in c("both", "ctrl", "stim")) {
  for (nm in n_min_grid) {
    agree <- matrix(NA, nrow = nrow(sp_full), ncol = length(donors), dimnames = list(NULL, donors))
    flip_in_donor <- agree
    for (d in donors) {
      f <- file.path(cache_dir, sprintf("donor_%s_%s_%s_nmin%d_pairs%d.rds", d, cs, "cell_type", nm, n_pairs))
      spd <- readRDS(f)$summary_pairs
      stopifnot(identical(spd$pair_id, sp_full$pair_id))
      okd <- is.finite(spd$r_resid) & spd$r_resid != 0 & is.finite(spd$r_global) & spd$r_global != 0
      agree[, d] <- ifelse(okd & ok_full, sign(spd$r_resid) == sign(sp_full$r_resid), NA)
      flip_in_donor[, d] <- ifelse(okd, sign(spd$r_global) != sign(spd$r_resid), NA)
    }
    n_fin <- rowSums(!is.na(agree))
    frac_agree <- rowSums(agree, na.rm = TRUE) / n_fin
    frac_flip <- rowSums(flip_in_donor, na.rm = TRUE) / rowSums(!is.na(flip_in_donor))
    for (g in names(pair_groups)) {
      sel <- pair_groups[[g]] & n_fin >= 6
      k <- n_fin[sel]
      rep_rows[[length(rep_rows) + 1]] <- data.frame(
        condition_subset = cs, n_min = nm, pair_group = g,
        n_pairs = sum(sel),
        mean_n_donors_evaluable = mean(k),
        mean_frac_donors_resid_sign_agrees = mean(frac_agree[sel]),
        frac_pairs_ge_75pct_donors_agree = mean(frac_agree[sel] >= 0.75),
        frac_pairs_all_donors_agree = mean(frac_agree[sel] == 1),
        frac_pairs_majority_donors_agree = mean(frac_agree[sel] > 0.5),
        mean_frac_donors_flip_within_donor = mean(frac_flip[sel], na.rm = TRUE),
        expected_frac_ge_75pct_if_random_signs = mean(stats::pbinom(ceiling(0.75 * k) - 1, k, 0.5, lower.tail = FALSE)),
        expected_frac_all_if_random_signs = mean(0.5^k),
        expected_frac_majority_if_random_signs = mean(stats::pbinom(floor(k / 2), k, 0.5, lower.tail = FALSE)),
        stringsAsFactors = FALSE
      )
    }
  }
}
donor_replication <- do.call(rbind, rep_rows)
donor_replication$row_id <- seq_len(nrow(donor_replication))
donor_replication <- donor_replication[, c("row_id", setdiff(colnames(donor_replication), "row_id"))]
utils::write.csv(donor_replication, file.path(out_dir, "donor_replication.csv"), row.names = FALSE)
toc("donor_replication", t0)

# ---------------------------------------------------------------------------
# 5. Pseudobulk analyses
# ---------------------------------------------------------------------------
t0 <- tic()
data_layer <- tryCatch(Seurat::GetAssayData(seu, assay = "RNA", layer = "data"),
                       error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "data"))
counts_layer <- tryCatch(Seurat::GetAssayData(seu, assay = "RNA", layer = "counts"),
                         error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "counts"))
data_feat <- data_layer[features, , drop = FALSE]
counts_feat <- counts_layer[features, , drop = FALSE]
lib_size <- Matrix::colSums(counts_layer)   # all genes, per cell
rm(data_layer, counts_layer)
invisible(gc())

unit_sets <- list(
  donor_x_cell_type = c("donor", "cell_type"),
  donor_x_condition_x_cell_type = c("donor", "condition", "cell_type")
)
unit_context_defs <- list(
  donor_x_cell_type = c("cell_type", "donor"),
  donor_x_condition_x_cell_type = c("cell_type", "donor", "condition", "condition_x_cell_type", "donor_x_condition", "donor_x_cell_type")
)

make_units <- function(keys) {
  lab <- do.call(paste, c(md[keys], sep = "|"))
  tab <- table(lab)
  keep_units <- names(tab)[tab >= min_cells_per_unit]
  cell_keep <- lab %in% keep_units
  f <- factor(lab[cell_keep], levels = keep_units)
  ind <- Matrix::sparse.model.matrix(~ 0 + f)          # cells x units
  colnames(ind) <- levels(f)
  n_per_unit <- Matrix::colSums(ind)
  mean_log <- as.matrix((data_feat[, cell_keep, drop = FALSE] %*% ind) %*% Matrix::Diagonal(x = 1 / n_per_unit))
  sum_counts <- as.matrix(counts_feat[, cell_keep, drop = FALSE] %*% ind)
  unit_lib <- as.numeric(lib_size[cell_keep] %*% ind)
  log_cpm <- log2(sweep(sum_counts, 2, unit_lib, "/") * 1e6 + 1)
  dimnames(mean_log) <- dimnames(log_cpm) <- list(features, levels(f))
  um <- unique(md[cell_keep, keys, drop = FALSE])
  um$unit <- do.call(paste, c(um[keys], sep = "|"))
  um <- um[match(levels(f), um$unit), , drop = FALSE]
  if (!("condition" %in% colnames(um))) um$condition <- NA_character_
  um$n_cells <- as.integer(n_per_unit[um$unit])
  um$library_size <- unit_lib[match(um$unit, levels(f))]
  rownames(um) <- NULL
  list(mean_log = mean_log, log_cpm = log_cpm, units = um, n_units_dropped = sum(tab < min_cells_per_unit),
       n_cells_dropped = sum(tab[tab < min_cells_per_unit]))
}

unit_context <- function(um, def) {
  switch(def,
    cell_type = um$cell_type,
    donor = um$donor,
    condition = um$condition,
    condition_x_cell_type = paste(um$condition, um$cell_type, sep = "|"),
    donor_x_condition = paste(um$donor, um$condition, sep = "|"),
    donor_x_cell_type = paste(um$donor, um$cell_type, sep = "|"),
    stop("unknown unit context ", def)
  )
}

# Permutation null for a pseudobulk run: within every context, and for every
# gene independently, shuffle the unit values among that context's units. This
# destroys within-context covariation between genes (the quantity behind r_c
# and r_resid) while preserving every context's mean and spread, so the pooled
# correlation keeps its between-context structure. The flip and reversal rates
# under this null show what noise alone produces with so few units per context.
permute_within_context <- function(mat, ctx, seed) {
  set.seed(seed)
  out <- mat
  for (lv in unique(ctx)) {
    idx <- which(ctx == lv)
    if (length(idx) < 2) next
    for (g in seq_len(nrow(mat))) out[g, idx] <- mat[g, idx[sample.int(length(idx))]]
  }
  out
}

null_rates_pseudobulk <- function(mat, ctx, n_obs, n_obs_resid, k_retained, R = n_perm, seed = seed_perm) {
  res <- parallel::mclapply(seq_len(R), function(r) {
    pm <- permute_within_context(mat, ctx, seed = seed + r)
    st <- sc_compute_pair_statistics_fast_pearson(
      expr_x = pm, expr_y = pm, pairs = pairs, context = ctx, eps = eps, majority = "count",
      n_min = n_min_units, cache_file = NULL, progress_every = 1e9, validate_n = 0
    )
    st <- recompute_heterogeneity(st)
    rates_from_summary(st$summary_pairs, n_obs = n_obs, n_obs_resid = n_obs_resid, k_retained = k_retained, B = 0L)
  }, mc.cores = min(n_cores, R), mc.preschedule = FALSE)
  for (x in res) if (!is.data.frame(x)) stop("A permutation run failed: ", paste(x, collapse = " "))
  res <- do.call(rbind, res)
  summ <- function(v) c(mean = mean(v, na.rm = TRUE), lo = stats::quantile(v, 0.025, na.rm = TRUE, names = FALSE), hi = stats::quantile(v, 0.975, na.rm = TRUE, names = FALSE))
  out <- data.frame(n_permutations = R, seed_perm = seed)
  for (col in c("flip_rate", "flip_rate_abs02", "simpson_rate_abs02", "mixed_sign_rate", "median_I2", "confident_flip_rate_cond", "confident_flip_rate_all", "n_pairs_both_significant", "n_confident_flips")) {
    s <- summ(res[[col]])
    out[[paste0("null_", col, "_mean")]] <- s[["mean"]]
    out[[paste0("null_", col, "_lo")]] <- s[["lo"]]
    out[[paste0("null_", col, "_hi")]] <- s[["hi"]]
  }
  out
}

pb_units_tables <- list()
pb_rows <- list()
pb_summaries <- list()
for (us in names(unit_sets)) {
  pb <- make_units(unit_sets[[us]])
  ut <- pb$units
  ut$unit_set <- us
  pb_units_tables[[us]] <- ut
  for (agg in c("mean_log", "log_cpm")) {
    mat <- pb[[agg]]
    for (def in unit_context_defs[[us]]) {
      ctx <- unit_context(pb$units, def)
      tab <- table(ctx)
      st <- run_pair_stats(mat, pairs, ctx, n_min = n_min_units,
                           cache_file = file.path(cache_dir, sprintf("pb_%s_%s_%s_pairs%d.rds", us, agg, def, n_pairs)),
                           force = FORCE)
      k_ret <- sum(tab >= n_min_units)
      n_ret <- sum(tab[tab >= n_min_units])
      estimable <- k_ret >= 1
      rates <- rates_from_summary(st$summary_pairs, n_obs = n_ret, n_obs_resid = n_ret, k_retained = k_ret,
                                  B = if (estimable) B_boot else 0L)
      null_df <- if (estimable) {
        null_rates_pseudobulk(mat, ctx, n_obs = n_ret, n_obs_resid = n_ret, k_retained = k_ret)
      } else {
        data.frame(n_permutations = 0L, seed_perm = seed_perm)
      }
      pb_rows[[length(pb_rows) + 1]] <- cbind(
        data.frame(unit_set = us, aggregation = agg, context_definition = def,
                   n_units = ncol(mat), n_units_dropped_lt_min_cells = pb$n_units_dropped,
                   n_cells_in_dropped_units = pb$n_cells_dropped, n_min_units = n_min_units,
                   k_contexts_total = length(tab), k_contexts_retained = k_ret, n_units_in_retained_contexts = n_ret,
                   min_units_per_context = as.integer(min(tab)), median_units_per_context = as.numeric(stats::median(as.integer(tab))),
                   max_units_per_context = as.integer(max(tab)),
                   within_context_estimable = estimable,
                   n_pairs_finite_pooled = sum(is.finite(st$summary_pairs$r_global)),
                   stringsAsFactors = FALSE),
        rates, null_df)
      pb_summaries[[paste(us, agg, def, sep = "__")]] <- st$summary_pairs
    }
  }
}
pb_rows <- lapply(pb_rows, function(d) { miss <- setdiff(unique(unlist(lapply(pb_rows, colnames))), colnames(d)); for (m in miss) d[[m]] <- NA; d })
pb_rows <- lapply(pb_rows, function(d) d[, colnames(pb_rows[[1]])])
pb_units <- do.call(rbind, lapply(pb_units_tables, function(u) u[, c("unit_set", "unit", "donor", "condition", "cell_type", "n_cells", "library_size")]))
utils::write.csv(pb_units, file.path(out_dir, "pseudobulk_units.csv"), row.names = FALSE)
pseudobulk <- do.call(rbind, pb_rows)
pseudobulk$row_id <- seq_len(nrow(pseudobulk))
pseudobulk <- pseudobulk[, c("row_id", setdiff(colnames(pseudobulk), "row_id"))]
utils::write.csv(pseudobulk, file.path(out_dir, "pseudobulk_rates.csv"), row.names = FALSE)

# Concordance between the cell-level (cell type as context, n_min = 200) and the
# pseudobulk (donor x cell type units, cell type as context) results on the same pairs.
sp_cell <- cell_summaries[[sprintf("cell_type_nmin%d", n_min_paper)]]
saveRDS(sp_cell, file.path(out_dir, "cell_pair_summary.rds"))
saveRDS(pb_summaries[["donor_x_cell_type__mean_log__cell_type"]],
        file.path(out_dir, "pseudobulk_pair_summary.rds"))
conc_rows <- list()
for (agg in c("mean_log", "log_cpm")) {
  sp_pb <- pb_summaries[[paste("donor_x_cell_type", agg, "cell_type", sep = "__")]]
  stopifnot(identical(sp_cell$pair_id, sp_pb$pair_id))
  ok <- is.finite(sp_cell$r_global) & is.finite(sp_cell$r_resid) & is.finite(sp_pb$r_global) & is.finite(sp_pb$r_resid) &
    sp_cell$r_global != 0 & sp_cell$r_resid != 0 & sp_pb$r_global != 0 & sp_pb$r_resid != 0
  flip_cell <- sign(sp_cell$r_global) != sign(sp_cell$r_resid)
  flip_pb <- sign(sp_pb$r_global) != sign(sp_pb$r_resid)
  add <- function(metric, value, n, subset) {
    conc_rows[[length(conc_rows) + 1]] <<- data.frame(aggregation = agg, subset = subset, metric = metric, value = value, n_pairs = n, stringsAsFactors = FALSE)
  }
  for (sub in c("all", "abs_r_global_cell_ge_0.2")) {
    s <- ok & (sub == "all" | abs(sp_cell$r_global) >= t_main)
    add("cor_r_global_cell_vs_pseudobulk", stats::cor(sp_cell$r_global[s], sp_pb$r_global[s]), sum(s), sub)
    add("cor_r_resid_cell_vs_pseudobulk", stats::cor(sp_cell$r_resid[s], sp_pb$r_resid[s]), sum(s), sub)
    add("sign_agreement_r_global", mean(sign(sp_cell$r_global[s]) == sign(sp_pb$r_global[s])), sum(s), sub)
    add("sign_agreement_r_resid", mean(sign(sp_cell$r_resid[s]) == sign(sp_pb$r_resid[s])), sum(s), sub)
    add("flip_rate_cell", mean(flip_cell[s]), sum(s), sub)
    add("flip_rate_pseudobulk", mean(flip_pb[s]), sum(s), sub)
    add("n_flip_cell", sum(flip_cell[s]), sum(s), sub)
    add("n_flip_both", sum(flip_cell[s] & flip_pb[s]), sum(s), sub)
    add("frac_cell_flips_also_flip_in_pseudobulk", if (sum(flip_cell[s]) > 0) sum(flip_cell[s] & flip_pb[s]) / sum(flip_cell[s]) else NA_real_, sum(flip_cell[s]), sub)
    add("frac_pseudobulk_flips_also_flip_at_cell_level", if (sum(flip_pb[s]) > 0) sum(flip_cell[s] & flip_pb[s]) / sum(flip_pb[s]) else NA_real_, sum(flip_pb[s]), sub)
    add("flip_rate_pseudobulk_among_cell_nonflips", if (sum(!flip_cell[s]) > 0) mean(flip_pb[s & !flip_cell]) else NA_real_, sum(s & !flip_cell), sub)
    # Simpson concordance among pairs eligible at both levels
    e <- s & abs(sp_cell$r_global) >= t_main & sp_cell$majority_sign != 0 & abs(sp_pb$r_global) >= t_main & sp_pb$majority_sign != 0
    add("simpson_rate_cell_among_jointly_eligible", mean(sp_cell$flag_simpson[e]), sum(e), sub)
    add("simpson_rate_pseudobulk_among_jointly_eligible", mean(sp_pb$flag_simpson[e]), sum(e), sub)
    add("frac_cell_simpson_also_simpson_in_pseudobulk", if (sum(sp_cell$flag_simpson[e]) > 0) sum(sp_cell$flag_simpson[e] & sp_pb$flag_simpson[e]) / sum(sp_cell$flag_simpson[e]) else NA_real_, sum(sp_cell$flag_simpson[e]), sub)
  }
}
concordance <- do.call(rbind, conc_rows)
concordance$row_id <- seq_len(nrow(concordance))
concordance <- concordance[, c("row_id", "aggregation", "subset", "metric", "value", "n_pairs")]
utils::write.csv(concordance, file.path(out_dir, "pseudobulk_concordance.csv"), row.names = FALSE)
toc("pseudobulk", t0)

# ---------------------------------------------------------------------------
# 6. PBMC 3k reference under the paper's retained-context convention
# ---------------------------------------------------------------------------
t0 <- tic()
pbmc_file <- file.path("results", "sc", "sc_100_pbmc_coarse_celltypes", "summary_pairs.csv")
uniform_metrics_file <- file.path("results", "revision", "uniform_table", "table1_metrics.csv")
pbmc_ref <- NULL
if (file.exists(uniform_metrics_file)) {
  uniform_metrics <- utils::read.csv(uniform_metrics_file, stringsAsFactors = FALSE)
  uniform_pbmc <- uniform_metrics[uniform_metrics$setting == "PBMC", , drop = FALSE]
  if (nrow(uniform_pbmc) == 1L) {
    pbmc_ref <- data.frame(
      dataset = "PBMC 3k (single donor)", context_definition = "cell_type_coarse",
      n_min = 200L, n_cells = uniform_pbmc$n_samples_analyzed,
      k_contexts_total = 7L, k_contexts_retained = uniform_pbmc$n_contexts_retained,
      n_cells_in_retained_contexts = uniform_pbmc$n_samples_analyzed,
      contexts_retained = "B;Mono;NK;T",
      n_pairs_flip_eval = uniform_pbmc$n_pairs_analyzed,
      flip_rate = uniform_pbmc$m5_pct_resid_flip_all / 100,
      flip_lo = uniform_pbmc$resid_flip_all_ci95_lo / 100,
      flip_hi = uniform_pbmc$resid_flip_all_ci95_hi / 100,
      n_pairs_flip_eval_abs02 = uniform_pbmc$n_pairs_abs_r_global_ge_0.2,
      flip_rate_abs02 = uniform_pbmc$m5_pct_resid_flip_ge_0.2 / 100,
      flip_abs02_lo = uniform_pbmc$resid_flip_ge_0.2_ci95_lo / 100,
      flip_abs02_hi = uniform_pbmc$resid_flip_ge_0.2_ci95_hi / 100,
      n_pairs_simpson_elig_abs02 = uniform_pbmc$n_pairs_eligible_simpson_ge_0.2,
      simpson_rate_abs02 = uniform_pbmc$m5_pct_simpson_reversal_ge_0.2 / 100,
      simpson_lo = uniform_pbmc$simpson_reversal_ge_0.2_ci95_lo / 100,
      simpson_hi = uniform_pbmc$simpson_reversal_ge_0.2_ci95_hi / 100,
      mixed_sign_rate = uniform_pbmc$m2_pct_mixed_sign_eps0.05 / 100,
      median_I2 = uniform_pbmc$m3_median_I2,
      confident_flip_rate_all = NA_real_, stringsAsFactors = FALSE
    )
  }
} else if (file.exists(pbmc_file)) {
  sp3k <- utils::read.csv(pbmc_file, stringsAsFactors = FALSE)
  inp3k <- readRDS(file.path("results", "sc", "sc_100_pbmc_coarse_celltypes", "inputs.rds"))
  tab3k <- inp3k$context_table
  keep3k <- tab3k >= inp3k$params$n_min
  pbmc_ref <- cbind(
    data.frame(dataset = "PBMC 3k (single donor)", context_definition = "cell_type_coarse", n_min = inp3k$params$n_min,
               n_cells = inp3k$n_cells, k_contexts_total = length(tab3k), k_contexts_retained = sum(keep3k),
               n_cells_in_retained_contexts = sum(tab3k[keep3k]),
               contexts_retained = paste(names(tab3k)[keep3k], collapse = ";"),
               stringsAsFactors = FALSE),
    rates_from_summary(sp3k, n_obs = inp3k$n_cells, n_obs_resid = sum(tab3k[keep3k]), k_retained = sum(keep3k))
  )
}
if (!is.null(pbmc_ref)) utils::write.csv(pbmc_ref, file.path(out_dir, "pbmc3k_reference.csv"), row.names = FALSE)
toc("pbmc3k_reference", t0)

# The reference line uses the same retained-context population as the revision.
pbmc_plot_ref <- if (!is.null(pbmc_ref)) pbmc_ref$flip_rate else NA_real_

# ---------------------------------------------------------------------------
# 7. Dataset summary and cell counts
# ---------------------------------------------------------------------------
cc <- as.data.frame(table(donor = md$donor, condition = md$condition, cell_type = md$cell_type), stringsAsFactors = FALSE)
colnames(cc)[4] <- "n_cells"
cc <- cc[order(cc$donor, cc$condition, cc$cell_type), ]
utils::write.csv(cc, file.path(out_dir, "cell_counts.csv"), row.names = FALSE)

ct_counts <- table(md$cell_type)
donor_counts <- table(md$donor)
cond_counts <- table(md$condition)
ds <- rbind(
  data.frame(item = "accession", value = accession),
  data.frame(item = "url", value = geo_url),
  data.frame(item = "supplementary_url", value = geo_suppl_url),
  data.frame(item = "samples", value = "GSM2560248 (control, batch 2.1); GSM2560249 (IFN-beta 6 h, batch 2.2)"),
  data.frame(item = "files", value = paste(c(suppl_files, raw_members), collapse = "; ")),
  data.frame(item = "reference", value = "Kang HM et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nat Biotechnol 2018;36:89-94"),
  data.frame(item = "n_cells_raw", value = qc$n_cells_raw),
  data.frame(item = "n_cells_raw_ctrl", value = qc$n_cells_raw_ctrl),
  data.frame(item = "n_cells_raw_stim", value = qc$n_cells_raw_stim),
  data.frame(item = "n_genes_raw", value = qc$n_genes_raw),
  data.frame(item = "n_demuxlet_singlets", value = qc$n_singlets),
  data.frame(item = "n_demuxlet_doublets_removed", value = qc$n_doublets),
  data.frame(item = "n_demuxlet_ambiguous_removed", value = qc$n_ambiguous),
  data.frame(item = "n_singlets_missing_celltype_removed", value = qc$n_singlets_missing_celltype),
  data.frame(item = "n_cells_singlet_labeled", value = qc$n_cells_singlet_labeled),
  data.frame(item = "n_cells_removed_by_qc", value = qc$n_removed_by_qc),
  data.frame(item = "n_cells_after_qc", value = qc$n_cells_after_qc),
  data.frame(item = "qc_min_features", value = min_features),
  data.frame(item = "qc_max_percent_mt", value = max_mt),
  data.frame(item = "normalization", value = "Seurat LogNormalize (scale factor 10,000)"),
  data.frame(item = "n_hvg", value = n_hvg),
  data.frame(item = "n_donors", value = length(donors)),
  data.frame(item = "donors", value = paste(sprintf("%s=%d", names(donor_counts), as.integer(donor_counts)), collapse = "; ")),
  data.frame(item = "n_conditions", value = length(conditions)),
  data.frame(item = "conditions", value = paste(sprintf("%s=%d", names(cond_counts), as.integer(cond_counts)), collapse = "; ")),
  data.frame(item = "n_cell_types_provided", value = length(cell_types)),
  data.frame(item = "cell_types_provided", value = paste(sprintf("%s=%d", names(ct_counts), as.integer(ct_counts)), collapse = "; ")),
  data.frame(item = "cell_type_coarse_map", value = paste(sprintf("%s->%s", coarse_map$cell_type, coarse_map$cell_type_coarse), collapse = "; ")),
  data.frame(item = "n_features_used", value = length(features)),
  data.frame(item = "n_marker_genes_used", value = length(feat_info$marker_genes_used)),
  data.frame(item = "n_hvg_used", value = sum(features %in% feat_info$hvg)),
  data.frame(item = "min_detect_frac", value = min_detect_frac),
  data.frame(item = "max_features", value = max_features),
  data.frame(item = "n_pairs", value = nrow(pairs)),
  data.frame(item = "seed_pairs", value = seed_pairs),
  data.frame(item = "eps", value = eps),
  data.frame(item = "abs_r_threshold_simpson", value = t_main),
  data.frame(item = "n_min_grid", value = paste(n_min_grid, collapse = ";")),
  data.frame(item = "bootstrap_B", value = B_boot),
  data.frame(item = "bootstrap_seed", value = seed_boot),
  data.frame(item = "pseudobulk_min_cells_per_unit", value = min_cells_per_unit),
  data.frame(item = "pseudobulk_n_min_units", value = n_min_units),
  data.frame(item = "quick_mode", value = QUICK),
  stringsAsFactors = FALSE
)
utils::write.csv(ds, file.path(out_dir, "dataset_summary.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 8. Figure
# ---------------------------------------------------------------------------
t0 <- tic()
col_flip <- "#E8601C"
col_keep <- "grey70"
col_ctrl <- "#1b9e77"
col_stim <- "#7570b3"

scatter_panel <- function(sp, main, xlab, ylab) {
  ok <- is.finite(sp$r_global) & is.finite(sp$r_resid) & sp$r_global != 0 & sp$r_resid != 0
  d <- sp[ok, c("r_global", "r_resid")]
  flip <- sign(d$r_global) != sign(d$r_resid)
  frac <- mean(flip)
  set.seed(seed_plot)
  show <- if (nrow(d) > n_plot_points) sort(sample.int(nrow(d), n_plot_points)) else seq_len(nrow(d))
  d <- d[show, ]
  flip <- flip[show]
  lim <- range(c(d$r_global, d$r_resid, -0.5, 0.5))
  plot(NA, xlim = lim, ylim = lim, xlab = xlab, ylab = ylab, main = main, asp = 1)
  abline(h = 0, v = 0, col = "grey85", lwd = 0.6)
  abline(0, 1, lty = 2, col = "grey45", lwd = 0.8)
  points(d$r_global[!flip], d$r_resid[!flip], pch = 16, cex = 0.25, col = grDevices::adjustcolor(col_keep, 0.45))
  points(d$r_global[flip], d$r_resid[flip], pch = 16, cex = 0.25, col = grDevices::adjustcolor(col_flip, 0.6))
  legend("topleft", bty = "n", pch = 16, pt.cex = 1, col = c(col_flip, col_keep), cex = 0.85,
         legend = c(sprintf("sign flip (%.1f%% of %s pairs)", 100 * frac, fmt_int(sum(ok))), "no flip"))
  invisible(frac)
}

panel_label <- function(lab) {
  usr <- par("usr")
  mtext(lab, side = 3, line = 1.3, at = usr[1] - 0.16 * diff(usr[1:2]), adj = 0, font = 2, cex = 1.35)
}

pdf_file <- file.path(out_dir, "Figure_MultiDonor.pdf")
grDevices::cairo_pdf(pdf_file, width = 10.5, height = 9.2)
layout(matrix(1:4, 2, 2, byrow = TRUE))
par(mar = c(4.4, 4.6, 3.2, 1.2), font.main = 2, cex.main = 1.2, cex.lab = 1.1, cex.axis = 0.95, mgp = c(2.6, 0.7, 0))

# A: cell level, cell type as context
scatter_panel(cell_summaries[[sprintf("cell_type_nmin%d", n_min_paper)]],
              main = "Cell Level, Cell Type as Context",
              xlab = expression(r[global] ~ "(retained-context cells)"), ylab = expression(r[resid] ~ "(cell-type means removed)"))
panel_label("A.")

# B: pseudobulk, donor x cell type units, cell type as context
scatter_panel(pb_summaries[["donor_x_cell_type__mean_log__cell_type"]],
              main = "Pseudobulk, Cell Type as Context",
              xlab = expression(r[global] ~ "(donor" %*% "cell-type units)"), ylab = expression(r[resid] ~ "(cell-type means removed)"))
panel_label("B.")

# C: per-donor residualized flip rates (each donor analyzed alone, coarse cell types, paper's n_min)
pdC <- per_donor[per_donor$celltype_definition == "cell_type_coarse" & per_donor$n_min == n_min_paper &
                   per_donor$condition_subset %in% c("ctrl", "stim"), ]
pdC <- pdC[order(pdC$donor, pdC$condition_subset), ]
pdC_estimable <- pdC$k_contexts_retained >= 2L
ymax <- max(c(pdC$flip_hi[pdC_estimable], pdC$flip_rate[pdC_estimable], if (is.finite(pbmc_plot_ref)) pbmc_plot_ref else 0), na.rm = TRUE) * 1.15
bar_mat <- matrix(ifelse(pdC_estimable, pdC$flip_rate, NA_real_), nrow = 2, dimnames = list(c("ctrl", "stim"), donors))
bp <- barplot(bar_mat, beside = TRUE, col = c(col_ctrl, col_stim), border = NA, ylim = c(0, ymax),
              names.arg = sub("^D", "", donors), xlab = "Donor (analyzed alone)", ylab = "Residualized flip rate",
              main = "Single-Donor Flip Rates", las = 1, yaxt = "n")
axis(2, at = pretty(c(0, ymax)), labels = paste0(100 * pretty(c(0, ymax)), "%"), las = 1)
lo <- matrix(ifelse(pdC_estimable, pdC$flip_lo, NA_real_), nrow = 2)
hi <- matrix(ifelse(pdC_estimable, pdC$flip_hi, NA_real_), nrow = 2)
arrows(bp, lo, bp, hi, angle = 90, code = 3, length = 0.02, col = "grey30", lwd = 0.8)
if (is.finite(pbmc_plot_ref)) abline(h = pbmc_plot_ref, lty = 2, lwd = 1.2, col = "black")
k_lab <- matrix(pdC$k_contexts_retained, nrow = 2)
lab_y <- ifelse(is.finite(pmax(hi, bar_mat)), pmax(hi, bar_mat) + 0.02 * ymax, 0.02 * ymax)
text(bp, lab_y, labels = ifelse(k_lab >= 2L, k_lab, "N/A"), cex = 0.6, col = "grey30")
if (is.finite(pbmc_plot_ref)) {
  legend("topright", bty = "n", cex = 0.85,
         legend = c("control cells", expression("IFN-" * beta * " cells"), sprintf("PBMC 3k, single donor (%.1f%%)", 100 * pbmc_plot_ref)),
         fill = c(col_ctrl, col_stim, NA), border = NA, lty = c(NA, NA, 2), lwd = c(NA, NA, 1.2), col = c(NA, NA, "black"))
} else {
  legend("topright", bty = "n", cex = 0.85, legend = c("control cells", expression("IFN-" * beta * " cells")),
         fill = c(col_ctrl, col_stim), border = NA)
}
panel_label("C.")

# D: rates by context definition at the cell level (paper's n_min)
pdD <- cell_level[cell_level$n_min == n_min_paper, ]
def_labels <- list(
  cell_type = "cell type (8)", cell_type_coarse = "cell type, coarse (6)", donor = "donor (8)", condition = "condition (2)",
  donor_x_cell_type = expression("donor" %*% "cell type (64)"), condition_x_cell_type = expression("condition" %*% "cell type (16)"),
  donor_x_condition = expression("donor" %*% "condition (16)"), donor_x_condition_x_cell_type = expression("donor" %*% "condition" %*% "cell type (128)")
)
ordD <- match(names(def_labels), pdD$context_definition)
pdD <- pdD[ordD, ]
nD <- nrow(pdD)
ypos <- rev(seq_len(nD))
par(mar = c(4.4, 11.5, 3.2, 1.2))
xmax <- max(c(pdD$flip_hi, pdD$simpson_hi, pdD$mixed_hi), na.rm = TRUE) * 1.08
plot(NA, xlim = c(0, xmax), ylim = c(0.5, nD + 0.5), yaxt = "n", ylab = "", xaxt = "n",
     xlab = "Rate among evaluable pairs", main = "Rates by Context Definition")
axis(1, at = pretty(c(0, xmax)), labels = paste0(100 * pretty(c(0, xmax)), "%"))
abline(h = ypos, col = "grey92")
for (i in seq_len(nD)) {
  axis(2, at = ypos[i], labels = def_labels[[i]], las = 1, cex.axis = 0.85, tick = FALSE)
}
off <- 0.22
segments(pdD$flip_lo, ypos + off, pdD$flip_hi, ypos + off, col = col_flip, lwd = 1.5)
points(pdD$flip_rate, ypos + off, pch = 16, col = col_flip, cex = 1.1)
segments(pdD$simpson_lo, ypos, pdD$simpson_hi, ypos, col = "#1f78b4", lwd = 1.5)
points(pdD$simpson_rate_abs02, ypos, pch = 17, col = "#1f78b4", cex = 1.1)
segments(pdD$mixed_lo, ypos - off, pdD$mixed_hi, ypos - off, col = "grey40", lwd = 1.5)
points(pdD$mixed_sign_rate, ypos - off, pch = 15, col = "grey40", cex = 1.0)
if (is.finite(pbmc_plot_ref)) abline(v = pbmc_plot_ref, lty = 2, col = "black", lwd = 1)
legend("topright", bty = "n", cex = 0.8, pch = c(16, 17, 15, NA), lty = c(NA, NA, NA, 2), col = c(col_flip, "#1f78b4", "grey40", "black"),
       legend = c("residualized flip", expression("Simpson reversal, |" * r[global] * "|" >= 0.2), expression("mixed sign (" * epsilon == 0.05 * ")"), "PBMC 3k flip rate"))
panel_label("D.")
invisible(dev.off())
toc("figure", t0)

# ---------------------------------------------------------------------------
# 9. Runtime and session information
# ---------------------------------------------------------------------------
total_secs <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
rt <- c(sprintf("total_seconds\t%.1f", total_secs),
        sprintf("%s\t%.1f", names(section_times), unlist(section_times)),
        sprintf("n_cores\t%d", n_cores), sprintf("quick_mode\t%s", QUICK), sprintf("finished\t%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
writeLines(rt, file.path(out_dir, "runtime.txt"))
writeLines(c(capture.output(sessionInfo())), file.path(out_dir, "sessionInfo.txt"))
