#!/usr/bin/env Rscript

# Genome Biology Brief Report -- Figure 1 HTRA3 Simpson example, split into
# FOUR standalone panels (the old fig_simpson.png baked all four into one image).
#
# Pair: hsa-miR-200c-3p vs HTRA3 (fixed to match the main text).
# Produces, in doc/genomebiology/figures/:
#   gb_fig1_htra3_global.pdf  (C) pooled scatter, global regression line (r=0.40)
#   gb_fig1_htra3_bars.pdf    (D) per-cohort within-cohort correlation bar chart
#   gb_fig1_htra3_means.pdf   (E) cohort-means scatter (r=0.52)
#   gb_fig1_htra3_resid.pdf   (F) pooled scatter with per-cohort regression lines
#                                 diverging from the global line (r_resid=-0.11)
#
# Run from repo root: Rscript scripts/gb_fig1_htra3_panels.R

source("00-paths.R")

f1 <- file.path(paths$mirtcga, "matchedNoGBM.Rda")
f2 <- file.path(paths$mirtcga, "matchedNoGBM.rda")
if (file.exists(f1)) load(f1) else if (file.exists(f2)) load(f2) else
  stop("Cannot find matchedNoGBM.{Rda,rda} in ", paths$mirtcga)
cc <- file.path(paths$mirtcga, "cohortColors.Rda")
if (file.exists(cc)) load(cc)
if (!exists("cohortColors")) stop("cohortColors not found")

out_dir <- file.path("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bestMir <- "hsa-miR-200c-3p"
bestGene <- "HTRA3"

mAll <- miRdata[bestMir, ]
gAll <- mRNAdata[bestGene, ]
globalBest <- cor(gAll, mAll)

localBest <- sapply(levels(cancerType), function(ct) {
  idx <- which(cancerType == ct)
  cor(mRNAdata[bestGene, idx], mAll[idx])
})
localBest <- localBest[is.finite(localBest)]
n_neg <- sum(localBest < 0); n_tot <- length(localBest)

# cohort means
coh <- as.character(cancerType)
mx <- tapply(gAll, coh, mean)
my <- tapply(mAll, coh, mean)
r_means <- cor(mx, my)

# mean-residualized correlation
gr <- gAll - mx[coh]
mr <- mAll - my[coh]
r_resid <- cor(gr, mr)

# cohorts with strongest |r| for the divergent-lines panel
K <- 5
ord <- order(abs(localBest), decreasing = TRUE)
keep <- names(localBest)[ord[seq_len(K)]]

message(sprintf("HTRA3 / miR-200c-3p: global r=%.3f, means r=%.3f, resid r=%.3f; negative in %d/%d cohorts",
                globalBest, r_means, r_resid, n_neg, n_tot))

# ---------------------------------------------------------------------------
# Panel C: pooled global scatter
# ---------------------------------------------------------------------------
pdf(file.path(out_dir, "gb_fig1_htra3_global.pdf"), width = 3.0, height = 2.15)
par(mar = c(3.2, 3.4, 1.6, 0.8), mgp = c(1.9, 0.6, 0))
plot(gAll, mAll, col = "grey75", pch = 16, cex = 0.4,
     xlab = "HTRA3", ylab = "miR-200c-3p",
     main = sprintf("Pooled: r = %+.2f", globalBest), cex.main = 1)
abline(lm(mAll ~ gAll), lwd = 2)
dev.off()

# ---------------------------------------------------------------------------
# Panel D: per-cohort correlation bar chart
# ---------------------------------------------------------------------------
lt <- sort(localBest)
pdf(file.path(out_dir, "gb_fig1_htra3_bars.pdf"), width = 3.3, height = 2.15)
par(mar = c(4.6, 4.1, 1.6, 0.8), mgp = c(2.5, 0.6, 0))
pts <- barplot(lt, col = cohortColors[names(lt)], border = NA,
               las = 2, cex.names = 0.45, ylab = "Within-cohort r",
               main = "Per-cohort correlation", cex.main = 1)
abline(h = 0, col = "grey50", lwd = 1)
abline(h = globalBest, lwd = 2, lty = 2, col = "black")
dev.off()

# ---------------------------------------------------------------------------
# Panel E: cohort-means scatter
# ---------------------------------------------------------------------------
pdf(file.path(out_dir, "gb_fig1_htra3_means.pdf"), width = 3.0, height = 2.15)
par(mar = c(3.2, 3.4, 1.6, 0.8), mgp = c(1.9, 0.6, 0))
plot(mx, my, col = cohortColors[names(mx)], pch = 16, cex = 1.1,
     xlab = "HTRA3 (cohort mean)", ylab = "miR-200c-3p (cohort mean)",
     main = sprintf("Cohort means: r = %+.2f", r_means), cex.main = 1, cex.lab = 0.8)
abline(lm(my ~ mx), lwd = 2)
dev.off()

# ---------------------------------------------------------------------------
# Panel F: pooled scatter with divergent per-cohort regression lines
# ---------------------------------------------------------------------------
pdf(file.path(out_dir, "gb_fig1_htra3_resid.pdf"), width = 3.0, height = 2.15)
par(mar = c(3.2, 3.4, 1.6, 0.8), mgp = c(1.9, 0.6, 0))
plot(gAll, mAll, col = "grey82", pch = 16, cex = 0.4,
     xlab = "HTRA3", ylab = "miR-200c-3p",
     main = sprintf("Within-cohort: r = %+.2f", r_resid), cex.main = 1)
abline(lm(mAll ~ gAll), lwd = 2.5, col = "black")
for (ct in keep) {
  idx <- which(cancerType == ct)
  points(gAll[idx], mAll[idx], col = cohortColors[ct], pch = 16, cex = 0.5)
  abline(lm(mAll[idx] ~ gAll[idx]), col = cohortColors[ct], lwd = 1.8)
}
dev.off()

message("Wrote 4 HTRA3 panels to ", out_dir)
