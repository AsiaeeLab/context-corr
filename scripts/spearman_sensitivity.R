#!/usr/bin/env Rscript

# Pearson vs Spearman sensitivity check on a manageable high-variance subset.
#
# Outputs:
# - results/spearman_sensitivity.csv
# - results/spearman_sensitivity_by_cohort.csv
# - figures/Figure_PearsonVsSpearman.pdf
#
# This script is designed to run from the repo root and uses the local path
# configuration in $HOME/Paths/context.json (via 00-paths.R).

source("00-paths.R")
source(file.path("R", "contextcorr.R"))

out_csv <- file.path("results", "spearman_sensitivity.csv")
out_cohort_csv <- file.path("results", "spearman_sensitivity_by_cohort.csv")
out_pdf <- file.path("figures", "Figure_PearsonVsSpearman.pdf")
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

pos_counts_pearson <- matrix(0L, nrow = G, ncol = M)
neg_counts_pearson <- matrix(0L, nrow = G, ncol = M)
k_valid_pearson <- matrix(0L, nrow = G, ncol = M)
pos_counts_spearman <- matrix(0L, nrow = G, ncol = M)
neg_counts_spearman <- matrix(0L, nrow = G, ncol = M)
k_valid_spearman <- matrix(0L, nrow = G, ncol = M)
cohort_concordance <- data.frame(
  cohort = character(),
  n = integer(),
  finite_pairs = integer(),
  cor_pearson_vs_spearman = numeric(),
  stringsAsFactors = FALSE
)

for (ct in cohorts) {
  idx <- idx_by_ct[[ct]]
  if (length(idx) < min_n_per_cohort) next
  Xct <- t(mRNA_sub[, idx, drop = FALSE])
  Yct <- t(miR_sub[, idx, drop = FALSE])
  r_ct_pearson <- suppressWarnings(stats::cor(Xct, Yct, method = "pearson"))
  r_ct_spearman <- suppressWarnings(stats::cor(Xct, Yct, method = "spearman"))

  paired <- is.finite(r_ct_pearson) & is.finite(r_ct_spearman)
  cohort_concordance <- rbind(
    cohort_concordance,
    data.frame(
      cohort = ct,
      n = length(idx),
      finite_pairs = sum(paired),
      cor_pearson_vs_spearman = suppressWarnings(stats::cor(
        r_ct_pearson[paired], r_ct_spearman[paired]
      )),
      stringsAsFactors = FALSE
    )
  )

  valid_pearson <- is.finite(r_ct_pearson) & abs(r_ct_pearson) < 1
  k_valid_pearson[valid_pearson] <- k_valid_pearson[valid_pearson] + 1L
  pos_pearson <- valid_pearson & (r_ct_pearson >= eps)
  neg_pearson <- valid_pearson & (r_ct_pearson <= -eps)
  pos_counts_pearson[pos_pearson] <- pos_counts_pearson[pos_pearson] + 1L
  neg_counts_pearson[neg_pearson] <- neg_counts_pearson[neg_pearson] + 1L

  valid_spearman <- is.finite(r_ct_spearman) & abs(r_ct_spearman) < 1
  k_valid_spearman[valid_spearman] <- k_valid_spearman[valid_spearman] + 1L
  pos_spearman <- valid_spearman & (r_ct_spearman >= eps)
  neg_spearman <- valid_spearman & (r_ct_spearman <= -eps)
  pos_counts_spearman[pos_spearman] <- pos_counts_spearman[pos_spearman] + 1L
  neg_counts_spearman[neg_spearman] <- neg_counts_spearman[neg_spearman] + 1L

  rm(
    Xct, Yct, r_ct_pearson, r_ct_spearman, paired,
    valid_pearson, pos_pearson, neg_pearson,
    valid_spearman, pos_spearman, neg_spearman
  )
  gc()
}

majority_sign <- function(pos_counts, neg_counts) {
  used <- (pos_counts + neg_counts) > 0L
  majority <- matrix(0L, nrow = nrow(pos_counts), ncol = ncol(pos_counts))
  majority[used] <- sign(pos_counts[used] - neg_counts[used])
  majority
}

maj_pearson <- majority_sign(pos_counts_pearson, neg_counts_pearson)
maj_spearman <- majority_sign(pos_counts_spearman, neg_counts_spearman)

eligible_pearson <-
  (k_valid_pearson >= min_k) &
  (abs(r_global_pearson) >= thr) &
  (maj_pearson != 0L) &
  (sign(r_global_pearson) != 0)
reversal_pearson <- eligible_pearson & (sign(r_global_pearson) != maj_pearson)
reversal_rate_pearson <- if (sum(eligible_pearson) > 0) {
  sum(reversal_pearson) / sum(eligible_pearson)
} else {
  NA_real_
}

eligible_spearman <-
  (k_valid_spearman >= min_k) &
  (abs(r_global_spearman) >= thr) &
  (maj_spearman != 0L) &
  (sign(r_global_spearman) != 0)
reversal_spearman <- eligible_spearman & (sign(r_global_spearman) != maj_spearman)
reversal_rate_spearman <- if (sum(eligible_spearman) > 0) {
  sum(reversal_spearman) / sum(eligible_spearman)
} else {
  NA_real_
}

utils::write.csv(cohort_concordance, out_cohort_csv, row.names = FALSE)

utils::write.csv(
  data.frame(
    G = G,
    M = M,
    domain_pairs = G * M,
    min_n_per_cohort = min_n_per_cohort,
    min_k = min_k,
    eps = eps,
    thr_abs_global = thr,
    cor_pooled_pearson_vs_spearman = pearson_vs_spearman,
    mean_cor_within_cohort_pearson_vs_spearman = mean(
      cohort_concordance$cor_pearson_vs_spearman,
      na.rm = TRUE
    ),
    simpson_reversal_rate_pearson = reversal_rate_pearson,
    reversal_pairs_pearson = sum(reversal_pearson),
    eligible_pairs_pearson = sum(eligible_pearson),
    simpson_reversal_rate_spearman = reversal_rate_spearman,
    reversal_pairs_spearman = sum(reversal_spearman),
    eligible_pairs_spearman = sum(eligible_spearman),
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
message("- ", out_cohort_csv)
message("- ", out_pdf)
