#!/usr/bin/env Rscript

# Pearson vs Spearman sensitivity check on a manageable high-variance subset.
#
# Outputs:
# - results/spearman_sensitivity.csv
# - doc/bioinformatics/figures/Figure_PearsonVsSpearman.pdf
#
# This script is designed to run from the repo root and uses the local path
# configuration in $HOME/Paths/context.json (via 00-paths.R).

source("00-paths.R")
source(file.path("R", "contextcorr.R"))

out_csv <- file.path("results", "spearman_sensitivity.csv")
out_pdf <- file.path("doc", "bioinformatics", "figures", "Figure_PearsonVsSpearman.pdf")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

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

G <- 1000
M <- 50
min_n_per_cohort <- 20
min_k <- 10
eps <- 0.05
thr <- 0.2

gene_var <- apply(mRNAdata, 1, var)
mir_var <- apply(miRdata, 1, var)
genes <- names(sort(gene_var, decreasing = TRUE))[1:G]
mirs <- names(sort(mir_var, decreasing = TRUE))[1:M]

mRNA_sub <- mRNAdata[genes, , drop = FALSE]
miR_sub <- miRdata[mirs, , drop = FALSE]

rm(gene_var, mir_var)
gc()

X <- t(mRNA_sub) # samples x genes
Y <- t(miR_sub)  # samples x miRs

r_global_pearson <- suppressWarnings(stats::cor(X, Y, method = "pearson"))
r_global_spearman <- suppressWarnings(stats::cor(X, Y, method = "spearman"))

pearson_vs_spearman <- suppressWarnings(stats::cor(as.vector(r_global_pearson), as.vector(r_global_spearman)))

cohorts <- levels(cancerType)
idx_by_ct <- lapply(cohorts, function(ct) which(cancerType == ct))
names(idx_by_ct) <- cohorts

pos_counts <- matrix(0L, nrow = G, ncol = M)
neg_counts <- matrix(0L, nrow = G, ncol = M)
k_valid <- matrix(0L, nrow = G, ncol = M)

for (ct in cohorts) {
  idx <- idx_by_ct[[ct]]
  if (length(idx) < min_n_per_cohort) next
  Xct <- t(mRNA_sub[, idx, drop = FALSE])
  Yct <- t(miR_sub[, idx, drop = FALSE])
  r_ct <- suppressWarnings(stats::cor(Xct, Yct, method = "spearman"))
  valid <- is.finite(r_ct) & abs(r_ct) < 1
  k_valid[valid] <- k_valid[valid] + 1L
  pos <- valid & (r_ct >= eps)
  neg <- valid & (r_ct <= -eps)
  pos_counts[pos] <- pos_counts[pos] + 1L
  neg_counts[neg] <- neg_counts[neg] + 1L
  rm(Xct, Yct, r_ct, valid, pos, neg)
  gc()
}

used <- (pos_counts + neg_counts) > 0L
maj <- matrix(0L, nrow = G, ncol = M)
maj[used] <- sign(pos_counts[used] - neg_counts[used])

eligible <- (k_valid >= min_k) & (abs(r_global_spearman) >= thr) & (maj != 0L) & (sign(r_global_spearman) != 0)
reversal <- eligible & (sign(r_global_spearman) != maj)
reversal_rate <- if (sum(eligible) > 0) sum(reversal) / sum(eligible) else NA_real_

utils::write.csv(
  data.frame(
    G = G,
    M = M,
    domain_pairs = G * M,
    min_n_per_cohort = min_n_per_cohort,
    min_k = min_k,
    eps = eps,
    thr_abs_global = thr,
    cor_pearson_vs_spearman = pearson_vs_spearman,
    simpson_reversal_rate_spearman = reversal_rate,
    eligible_pairs = sum(eligible),
    stringsAsFactors = FALSE
  ),
  out_csv,
  row.names = FALSE
)

pdf(out_pdf, width = 4.5, height = 4.5)
set.seed(1)
ix <- which(is.finite(r_global_pearson) & is.finite(r_global_spearman))
if (length(ix) > 20000) ix <- sample(ix, 20000)
plot(
  as.vector(r_global_pearson)[ix],
  as.vector(r_global_spearman)[ix],
  pch = 16,
  col = rgb(0, 0, 0, 0.07),
  xlab = "Pearson r (global)",
  ylab = "Spearman rho (global)",
  main = sprintf("Pearson vs Spearman (corr=%.3f)", pearson_vs_spearman)
)
abline(a = 0, b = 1, col = "#d95f02", lwd = 2)
dev.off()

message("Wrote:")
message("- ", out_csv)
message("- ", out_pdf)
