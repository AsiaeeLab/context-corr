#!/usr/bin/env Rscript

# Genome Biology Brief Report figure splits:
#   - doc/genomebiology/figures/gb_fig1A_globalhist.pdf
#       Standalone pan-cancer global-correlation histogram (new Fig 1 panel A).
#   - doc/genomebiology/figures/gb_figS2_mir122_bars.pdf
#       miR-122 within-cohort bar panels (Supplementary Fig S2), moved out of
#       the main text.
#
# These are the cleanly separated pieces of the old combined
# Figure_GlobalDistribution.pdf, which baked panels A/B/C into one image.
#
# Run from repo root: Rscript scripts/gb_fig1_histogram_and_mir122.R

source("00-paths.R")

out_dir <- file.path("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
hist_pdf <- file.path(out_dir, "gb_fig1A_globalhist.pdf")
bars_pdf <- file.path(out_dir, "gb_figS2_mir122_bars.pdf")

f1 <- file.path(paths$mirtcga, "matchedNoGBM.Rda")
f2 <- file.path(paths$mirtcga, "matchedNoGBM.rda")
if (file.exists(f1)) {
  load(f1)
} else if (file.exists(f2)) {
  load(f2)
} else {
  stop("Cannot find matchedNoGBM.Rda or matchedNoGBM.rda in ", paths$mirtcga)
}
# cohortColors lives in a sibling file
cc <- file.path(paths$mirtcga, "cohortColors.Rda")
if (file.exists(cc)) load(cc)

set.seed(1)

# ---------------------------------------------------------------------------
# Panel A (Fig 1): distribution of global correlations on a random subset.
# ---------------------------------------------------------------------------
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

pdf(hist_pdf, width = 4.6, height = 3.2)
par(mar = c(4, 4, 1, 1))
hist(
  as.vector(r_global),
  breaks = 120,
  xlab = "Global Pearson correlation",
  main = "",
  col = "steelblue",
  border = "grey20",
  lwd = 0.3
)
dev.off()
message("Wrote: ", hist_pdf)

# ---------------------------------------------------------------------------
# Supplementary Fig S2: miR-122 within-cohort bars (panels B/C of the old fig).
# ---------------------------------------------------------------------------
mir_ext <- "hsa-miR-122-3p"
g_pos <- "CFHR2"
g_neg <- "SLC25A36"

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
  x <- mRNAdata[gene, ]; y <- miRdata[mir, ]
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

bar_cols <- if (exists("cohortColors")) cohortColors else NULL

# Compact side-by-side layout (no baked-in A/B sub-letters): this panel now
# sits as panel (B) of main Figure 1, so it must not carry its own letters.
pdf(bars_pdf, width = 9.6, height = 3.4)
layout(matrix(c(1, 2), nrow = 1))
opar <- par(mar = c(5.5, 4, 2, 0.6))
pts <- barplot(
  r_pos, col = if (!is.null(bar_cols)) bar_cols[names(r_pos)] else "grey60",
  border = NA, xaxt = "n", ylab = "Within-cohort r",
  main = sprintf("vs %s (global r = %.2f)", g_pos, r_pos_global),
  cex.main = 0.95)
abline(h = 0, col = "grey50", lwd = 1.5)
mtext(names(r_pos), side = 1, at = pts, line = 0.5, las = 2, cex = 0.7)

pts <- barplot(
  r_neg, col = if (!is.null(bar_cols)) bar_cols[names(r_neg)] else "grey60",
  border = NA, xaxt = "n", ylab = "Within-cohort r",
  main = sprintf("vs %s (global r = %.2f)", g_neg, r_neg_global),
  cex.main = 0.95)
abline(h = 0, col = "grey50", lwd = 1.5)
mtext(names(r_neg), side = 1, at = pts, line = 0.5, las = 2, cex = 0.7)
par(opar)
dev.off()
message("Wrote: ", bars_pdf)
