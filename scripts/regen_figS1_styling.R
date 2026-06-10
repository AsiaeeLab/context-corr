#!/usr/bin/env Rscript

# Regenerate Fig S1 sub-panels with matched styling (all bold titles) and a
# zoomed y-axis on panel A so the eps curves are distinguishable instead of
# a single flat line near 1.0. Reads from cached results so no heavy compute.

suppressPackageStartupMessages({
  library(utils)
})

fig_dir <- file.path("figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Top row: A (mixed-sign) and B (Simpson reversal) -----

df_all <- utils::read.csv(file.path("results", "robustness_signflip_summary.csv"),
                          stringsAsFactors = FALSE)
eps_grid <- sort(unique(df_all$eps))
t_grid <- sort(unique(df_all$threshold))
cols <- c("#1b9e77", "#d95f02", "#7570b3")

out_top <- file.path(fig_dir, "Figure_RobustnessSignFlip.pdf")
cairo_pdf(out_top, width = 8.5, height = 4.0)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1), font.main = 2, cex.main = 1.05)

# Panel A: zoomed y-axis (the fraction is always close to 1; the interest is
# the separation between eps curves)
y_lo <- 0.97
plot(NA,
  xlim = range(t_grid),
  ylim = c(y_lo, 1.00),
  xlab = expression(paste("Effect-size threshold on ", max[c](abs(r[c])))),
  ylab = "Fraction mixed-sign",
  main = "A  Mixed-sign vs effect size",
  yaxs = "i"
)
abline(h = seq(y_lo, 1.0, by = 0.005), col = "grey92", lty = 1)
box()
for (i in seq_along(eps_grid)) {
  eps <- eps_grid[i]
  dat <- df_all[df_all$eps == eps, ]
  lines(dat$threshold, dat$frac_mixed, lwd = 2, col = cols[i])
  points(dat$threshold, dat$frac_mixed, pch = 16, col = cols[i], cex = 0.9)
}
legend("bottomleft",
  legend = as.expression(lapply(eps_grid, function(e) bquote(epsilon == .(e)))),
  col = cols, lwd = 2, pch = 16, bty = "n", cex = 0.95
)

# Panel B: unchanged
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
  points(dat$threshold, dat$frac_simpson, pch = 16, col = cols[i], cex = 0.9)
}
legend("topright",
  legend = as.expression(lapply(eps_grid, function(e) bquote(epsilon == .(e)))),
  col = cols, lwd = 2, pch = 16, bty = "n", cex = 0.95
)

dev.off()
message("Wrote: ", out_top)

# ----- Bottom row: C (I^2 distribution) and D (I^2 vs effect size) -----

h <- readRDS(file.path("results", "cache", "heterogeneity_I2_FDR.rds"))
ok <- h$ok
k_valid <- h$k_valid
I2 <- h$I2
r_fixed <- h$r_fixed
min_k <- h$params$min_k

idx_base <- as.logical(ok) & (k_valid >= min_k)
idx_plot <- which(idx_base & is.finite(I2))
if (length(idx_plot) > 60000) {
  set.seed(1)
  idx_plot <- sample(idx_plot, 60000)
}

out_bot <- file.path(fig_dir, "Figure_Heterogeneity_I2_CD.pdf")
cairo_pdf(out_bot, width = 8.5, height = 4.0)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1), font.main = 2, cex.main = 1.05)

# Panel C: histogram of I^2 (bold title via par(font.main=2))
hist(I2[idx_plot],
  breaks = 50,
  col = "grey85", border = "grey60",
  xlab = expression(I^2),
  main = expression(paste("C  Heterogeneity (", I^2, ") distribution"))
)

# Panel D: scatter |r_fixed| vs I^2
plot(abs(r_fixed[idx_plot]), I2[idx_plot],
  pch = 16, col = adjustcolor("black", alpha.f = 0.18), cex = 0.4,
  xlab = expression(paste(abs(r["fixed"]), " (within-cohort pooled)")),
  ylab = expression(I^2),
  main = expression(paste("D  ", I^2, " vs effect size"))
)
abline(h = 0.75, col = "#d95f02", lwd = 2, lty = 2)

dev.off()
message("Wrote: ", out_bot)
