#!/usr/bin/env Rscript

# Main Figure 1: global correlation distribution + tissue-driven extreme examples.
#
# Outputs:
# - doc/bioinformatics/figures/Figure_GlobalDistribution.pdf
#
# This script is designed to run from the repo root and uses the local path
# configuration in $HOME/Paths/context.json (via 00-paths.R).

source("00-paths.R")

out_pdf <- file.path("doc", "bioinformatics", "figures", "Figure_GlobalDistribution.pdf")
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

# Panel A: distribution of global correlations on a large random subset (avoid
# computing all ~23M pairs).
G <- 5000
M <- 200

genes <- sample(rownames(mRNAdata), G)
mirs <- sample(rownames(miRdata), M)

mRNA_sub <- mRNAdata[genes, , drop = FALSE]
miR_sub <- miRdata[mirs, , drop = FALSE]

scaled_mir_global <- scale(t(miR_sub)) / (ncol(miR_sub) - 1)
std_rna_global <- t(scale(t(mRNA_sub)))
r_global <- std_rna_global %*% scaled_mir_global

rm(mRNA_sub, miR_sub, scaled_mir_global, std_rna_global)
gc()

# Panels B/C: tissue-driven extremes (as in Kevin's original narrative).
mir_ext <- "hsa-miR-122-3p"
g_pos <- "CFHR2"
g_neg <- "SLC25A36"

if (!(mir_ext %in% rownames(miRdata))) stop("Missing miR: ", mir_ext)
if (!(g_pos %in% rownames(mRNAdata))) stop("Missing gene: ", g_pos)
if (!(g_neg %in% rownames(mRNAdata))) stop("Missing gene: ", g_neg)

cor_by_cohort <- function(gene, mir) {
  sapply(levels(cancerType), function(ct) {
    idx <- which(cancerType == ct)
    x <- mRNAdata[gene, idx]
    y <- miRdata[mir, idx]
    if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) return(NA_real_)
    if (!is.finite(stats::sd(y)) || stats::sd(y) == 0) return(NA_real_)
    suppressWarnings(stats::cor(x, y))
  })
}

r_global_pair <- function(gene, mir) {
  x <- mRNAdata[gene, ]
  y <- miRdata[mir, ]
  if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) return(NA_real_)
  if (!is.finite(stats::sd(y)) || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y))
}

r_pos <- cor_by_cohort(g_pos, mir_ext)
r_neg <- cor_by_cohort(g_neg, mir_ext)

r_pos_global <- r_global_pair(g_pos, mir_ext)
r_neg_global <- r_global_pair(g_neg, mir_ext)

na_pos <- is.na(r_pos)
na_neg <- is.na(r_neg)

pdf(out_pdf, width = 7.0, height = 7.0)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(1.3, 1, 1))

hist(
  as.vector(r_global),
  breaks = 120,
  xlab = "Global Pearson correlation (sampled gene × miR pairs)",
  main = "A  Pan-cancer correlation distribution (sampled)",
  col = "grey80",
  border = "white"
)

opar <- par(mar = c(8, 4, 2, 1))
pts <- barplot(
  r_pos,
  col = cohortColors[names(r_pos)],
  xaxt = "n",
  ylab = "Correlation",
  main = sprintf(
    "B  %s vs %s (within-cohort; global r = %.3f; missing = %d/%d)",
    mir_ext, g_pos, r_pos_global, sum(na_pos), length(na_pos)
  )
)
abline(h = 0, col = "grey50", lwd = 2)
mtext(names(r_pos), side = 1, at = pts, line = 1, las = 2, cex = 0.6)

pts <- barplot(
  r_neg,
  col = cohortColors[names(r_neg)],
  xaxt = "n",
  ylab = "Correlation",
  main = sprintf(
    "C  %s vs %s (within-cohort; global r = %.3f; missing = %d/%d)",
    mir_ext, g_neg, r_neg_global, sum(na_neg), length(na_neg)
  )
)
abline(h = 0, col = "grey50", lwd = 2)
mtext(names(r_neg), side = 1, at = pts, line = 1, las = 2, cex = 0.6)
par(opar)

dev.off()

message("Wrote:")
message("- ", out_pdf)
