#!/usr/bin/env Rscript

# Effect-size–conditioned robustness checks for context-dependent correlation:
# - mixed-sign fraction vs max(|r_cohort|)
# - Simpson reversal fraction vs |r_global|
#
# Outputs:
# - results/robustness_signflip_summary.csv
# - doc/bioinformatics/figures/Figure_RobustnessSignFlip.pdf
#
# This script is designed to run from the repo root and uses the local path
# configuration in $HOME/Paths/context.json (via 00-paths.R).

source("00-paths.R")
source(file.path("R", "contextcorr.R"))

out_csv <- file.path("results", "robustness_signflip_summary.csv")
out_pdf <- file.path("doc", "genomebiology", "figures", "Figure_RobustnessSignFlip.pdf")
cache_rds <- file.path("results", "cache", "robustness_signflip.rds")
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
eps_grid <- c(0, 0.05, 0.10)
t_grid <- c(0.0, 0.05, 0.10, 0.20, 0.30)

# Size choices: large enough to be convincing, small enough to be tractable.
G <- 5000
M <- 200

cache_params <- list(
  G = G,
  M = M,
  min_n_per_cohort = min_n_per_cohort,
  min_k = min_k,
  eps_grid = eps_grid,
  t_grid = t_grid,
  matched_file = if (file.exists(f1)) normalizePath(f1) else normalizePath(f2)
)

df_all <- NULL
if (file.exists(cache_rds)) {
  cache <- readRDS(cache_rds)
  if (is.list(cache) && identical(cache$params, cache_params) && is.data.frame(cache$df_all)) {
    df_all <- cache$df_all
  }
}

if (is.null(df_all)) {
  gene_var <- apply(mRNAdata, 1, var)
  mir_var <- apply(miRdata, 1, var)

  genes <- names(sort(gene_var, decreasing = TRUE))[1:G]
  mirs <- names(sort(mir_var, decreasing = TRUE))[1:M]

  mRNA_sub <- mRNAdata[genes, , drop = FALSE]
  miR_sub <- miRdata[mirs, , drop = FALSE]

  rm(gene_var, mir_var)
  gc()

  # Global correlations (pan-cancer)
  scaled_mir_global <- scale(t(miR_sub)) / (ncol(miR_sub) - 1)
  std_rna_global <- t(scale(t(mRNA_sub)))
  r_global <- std_rna_global %*% scaled_mir_global

  rm(scaled_mir_global, std_rna_global)
  gc()

  cohorts <- levels(cancerType)
  idx_by_ct <- lapply(cohorts, function(ct) which(cancerType == ct))
  names(idx_by_ct) <- cohorts
  n_by_ct <- vapply(idx_by_ct, length, integer(1))

  k_valid <- matrix(0L, nrow = G, ncol = M, dimnames = list(genes, mirs))
  max_abs <- matrix(0, nrow = G, ncol = M, dimnames = list(genes, mirs))

  pos_counts <- lapply(eps_grid, function(eps) matrix(0L, nrow = G, ncol = M, dimnames = list(genes, mirs)))
  neg_counts <- lapply(eps_grid, function(eps) matrix(0L, nrow = G, ncol = M, dimnames = list(genes, mirs)))

  for (ct in cohorts) {
    idx <- idx_by_ct[[ct]]
    if (length(idx) < min_n_per_cohort) next

    X <- mRNA_sub[, idx, drop = FALSE]
    Y <- miR_sub[, idx, drop = FALSE]

    scaled_mir <- scale(t(Y)) / (ncol(Y) - 1)
    std_rna <- t(scale(t(X)))
    r_ct <- std_rna %*% scaled_mir

    valid <- is.finite(r_ct) & abs(r_ct) < 1
    k_valid[valid] <- k_valid[valid] + 1L
    max_abs[valid] <- pmax(max_abs[valid], abs(r_ct[valid]))

    for (i in seq_along(eps_grid)) {
      eps <- eps_grid[i]
      pos <- valid & (r_ct >= eps)
      neg <- valid & (r_ct <= -eps)
      pos_counts[[i]][pos] <- pos_counts[[i]][pos] + 1L
      neg_counts[[i]][neg] <- neg_counts[[i]][neg] + 1L
    }

    rm(X, Y, scaled_mir, std_rna, r_ct, valid, pos, neg)
    gc()
  }

flatten <- function(mat) as.vector(mat)

  df_all <- data.frame()

  for (i in seq_along(eps_grid)) {
    eps <- eps_grid[i]
    pos <- pos_counts[[i]]
    neg <- neg_counts[[i]]
    used <- (pos + neg) > 0L
    maj <- matrix(0L, nrow = G, ncol = M)
    maj[used] <- sign(pos[used] - neg[used])

    mixed <- (pos > 0L) & (neg > 0L)

    for (t in t_grid) {
      elig_mixed <- (k_valid >= min_k) & (max_abs >= t)
      frac_mixed <- if (sum(elig_mixed) > 0) sum(mixed & elig_mixed) / sum(elig_mixed) else NA_real_

      elig_simpson <- (k_valid >= min_k) & (abs(r_global) >= t) & (maj != 0L) & (sign(r_global) != 0)
      simpson <- elig_simpson & (sign(r_global) != maj)
      frac_simpson <- if (sum(elig_simpson) > 0) sum(simpson) / sum(elig_simpson) else NA_real_

      df_all <- rbind(
        df_all,
        data.frame(
          eps = eps,
          threshold = t,
          min_k = min_k,
          min_n_per_cohort = min_n_per_cohort,
          domain_pairs = G * M,
          eligible_mixed = sum(elig_mixed),
          frac_mixed = frac_mixed,
          eligible_simpson = sum(elig_simpson),
          frac_simpson = frac_simpson,
          stringsAsFactors = FALSE
        )
      )
  }
}

  saveRDS(list(params = cache_params, df_all = df_all), cache_rds)
}

utils::write.csv(df_all, out_csv, row.names = FALSE)

pdf(out_pdf, width = 8.5, height = 4.0)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

cols <- c("#1b9e77", "#d95f02", "#7570b3")

# Panel A: mixed sign vs max abs within-cohort
plot(NA,
  xlim = range(t_grid),
  ylim = c(0, 1),
  xlab = expression(paste("Effect-size threshold on ", max[c](abs(r[c])))),
  ylab = "Fraction mixed-sign",
  main = "A  Mixed-sign vs effect size"
)
for (i in seq_along(eps_grid)) {
  eps <- eps_grid[i]
  dat <- df_all[df_all$eps == eps, ]
  lines(dat$threshold, dat$frac_mixed, lwd = 2, col = cols[i])
  points(dat$threshold, dat$frac_mixed, pch = 16, col = cols[i], cex = 0.7)
}
legend("topright",
  legend = paste0("eps = ", eps_grid),
  col = cols, lwd = 2, bty = "n", cex = 0.9
)

# Panel B: Simpson reversal vs |r_global|
plot(NA,
  xlim = range(t_grid),
  ylim = c(0, 1),
  xlab = expression(paste("Effect-size threshold on ", abs(r[global]))),
  ylab = "Fraction Simpson reversal",
  main = "B  Simpson reversals vs effect size"
)
for (i in seq_along(eps_grid)) {
  eps <- eps_grid[i]
  dat <- df_all[df_all$eps == eps, ]
  lines(dat$threshold, dat$frac_simpson, lwd = 2, col = cols[i])
  points(dat$threshold, dat$frac_simpson, pch = 16, col = cols[i], cex = 0.7)
}
legend("topright",
  legend = paste0("eps = ", eps_grid),
  col = cols, lwd = 2, bty = "n", cex = 0.9
)

dev.off()

message("Wrote:")
message("- ", out_csv)
message("- ", out_pdf)
