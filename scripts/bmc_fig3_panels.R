#!/usr/bin/env Rscript

# BMC Genomics Figure 3 as FOUR separate panels (A-D) for the LaTeX subfigure
# layout. Panel letters are supplied by the LaTeX subcaptions, so the titles
# here carry NO "A "/"B " prefix. Titles are bold to match the other figures.
# Reads cached results (cheap, no heavy compute). Panel E (let-7a heatmap) is
# the plain pheatmap output (no baked-in letter).
# Writes the four panels into figures/.
#
# Run from repo root: Rscript scripts/bmc_fig3_panels.R

suppressPackageStartupMessages({
  library(utils)
  library(grDevices)
})

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols <- c("#1b9e77", "#d95f02", "#7570b3")
W <- 4.8; H <- 3.4
set_par <- function() par(mar = c(4.2, 4.2, 2.0, 1), font.main = 2,
                          cex.main = 1.25, cex.lab = 1.1, cex.axis = 1.0)

## ---- Panels A (mixed-sign) and B (Simpson reversal): cached CSV ----
df_all <- utils::read.csv(file.path("results", "robustness_signflip_summary.csv"),
                          stringsAsFactors = FALSE)
eps_grid <- sort(unique(df_all$eps))
t_grid <- sort(unique(df_all$threshold))
eps_legend <- as.expression(lapply(eps_grid, function(e) bquote(epsilon == .(e))))

cairo_pdf(file.path(out_dir, "fig3A_mixsign.pdf"), width = W, height = H)
set_par()
y_lo <- 0.97
plot(NA, xlim = range(t_grid), ylim = c(y_lo, 1.00),
     xlab = expression(paste("Effect-size threshold on ", max[c](abs(r[c])))),
     ylab = "Fraction mixed-sign",
     main = "Mixed-sign vs effect size", yaxs = "i")
abline(h = seq(y_lo, 1.0, by = 0.005), col = "grey92", lty = 1)
box()
for (i in seq_along(eps_grid)) {
  dat <- df_all[df_all$eps == eps_grid[i], ]
  lines(dat$threshold, dat$frac_mixed, lwd = 2, col = cols[i])
  points(dat$threshold, dat$frac_mixed, pch = 16, col = cols[i], cex = 0.9)
}
legend("bottomleft", legend = eps_legend, col = cols, lwd = 2, pch = 16, bty = "n", cex = 0.95)
invisible(dev.off())

cairo_pdf(file.path(out_dir, "fig3B_reversal.pdf"), width = W, height = H)
set_par()
plot(NA, xlim = range(t_grid), ylim = c(0, 1),
     xlab = expression(paste("Effect-size threshold on ", abs(r[global]))),
     ylab = "Fraction Simpson reversal",
     main = "Simpson reversals vs effect size")
for (i in seq_along(eps_grid)) {
  dat <- df_all[df_all$eps == eps_grid[i], ]
  lines(dat$threshold, dat$frac_simpson, lwd = 2, col = cols[i])
  points(dat$threshold, dat$frac_simpson, pch = 16, col = cols[i], cex = 0.9)
}
legend("topright", legend = eps_legend, col = cols, lwd = 2, pch = 16, bty = "n", cex = 0.95)
invisible(dev.off())

## ---- Panels C (I^2 distribution) and D (I^2 vs effect size): cached rds ----
h <- readRDS(file.path("results", "cache", "heterogeneity_I2_FDR.rds"))
ok <- h$ok; k_valid <- h$k_valid; I2 <- h$I2; r_fixed <- h$r_fixed; min_k <- h$params$min_k
idx_base <- as.logical(ok) & (k_valid >= min_k)
idx_plot <- which(idx_base & is.finite(I2))
if (length(idx_plot) > 60000) { set.seed(1); idx_plot <- sample(idx_plot, 60000) }

cairo_pdf(file.path(out_dir, "fig3C_i2dist.pdf"), width = W, height = H)
set_par()
hist(I2[idx_plot], breaks = 50, col = "grey85", border = "grey60",
     xlab = expression(I^2),
     main = expression(bold(paste("Heterogeneity (", I^2, ") distribution"))))
invisible(dev.off())

cairo_pdf(file.path(out_dir, "fig3D_i2vseffect.pdf"), width = W, height = H)
set_par()
plot(abs(r_fixed[idx_plot]), I2[idx_plot],
     pch = 16, col = adjustcolor("black", alpha.f = 0.18), cex = 0.4,
     xlab = expression(paste(abs(r["fixed"]), " (within-cohort pooled)")),
     ylab = expression(I^2),
     main = expression(bold(paste(I^2, " vs effect size"))))
abline(h = 0.75, col = "#d95f02", lwd = 2, lty = 2)
invisible(dev.off())

cat("Wrote fig3 panels A-D to", out_dir, "\n")
