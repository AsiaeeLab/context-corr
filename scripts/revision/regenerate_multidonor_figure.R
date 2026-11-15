#!/usr/bin/env Rscript

# Regenerate Figure S5 from the cached full multi-donor analysis without
# rerunning the ten-minute analysis. The PBMC 3k reference line is taken from
# the uniform-convention table, for which r_global and r_resid use the same
# 2,655 cells in retained contexts.

out_dir <- file.path("results", "revision", "multidonor")
out_file <- file.path(out_dir, "Figure_MultiDonor.pdf")

cell_summary <- readRDS(file.path(out_dir, "cell_pair_summary.rds"))
pb_summary <- readRDS(file.path(out_dir, "pseudobulk_pair_summary.rds"))
per_donor <- read.csv(file.path(out_dir, "per_donor_rates.csv"), stringsAsFactors = FALSE)
cell_level <- read.csv(file.path(out_dir, "cell_level_rates.csv"), stringsAsFactors = FALSE)
uniform <- read.csv(file.path("results", "revision", "uniform_table", "table1_metrics.csv"),
                    stringsAsFactors = FALSE)
pbmc_ref <- uniform$m5_pct_resid_flip_all[uniform$setting == "PBMC"] / 100
stopifnot(length(pbmc_ref) == 1L, is.finite(pbmc_ref))

fmt_int <- function(x) format(x, big.mark = ",", trim = TRUE)
col_flip <- "#E67E22"
col_keep <- "#7F8C8D"
col_ctrl <- "#1B4F72"
col_stim <- "#85C1E9"
seed_plot <- 20260906L
n_plot_points <- 50000L

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
  grid(col = "grey92", lty = 1)
  abline(h = 0, v = 0, col = "grey85", lwd = 0.6)
  abline(0, 1, lty = 2, col = "grey45", lwd = 0.8)
  points(d$r_global[!flip], d$r_resid[!flip], pch = 16, cex = 0.25,
         col = grDevices::adjustcolor(col_keep, 0.45))
  points(d$r_global[flip], d$r_resid[flip], pch = 16, cex = 0.25,
         col = grDevices::adjustcolor(col_flip, 0.6))
  legend("topleft", bty = "n", pch = 16, pt.cex = 1,
         col = c(col_flip, col_keep), cex = 0.85,
         legend = c(sprintf("sign flip (%.1f%% of %s pairs)", 100 * frac, fmt_int(sum(ok))),
                    "no flip"))
}

panel_label <- function(lab) {
  usr <- par("usr")
  mtext(lab, side = 3, line = 1.3, at = usr[1] - 0.16 * diff(usr[1:2]),
        adj = 0, font = 2, cex = 1.35)
}

grDevices::cairo_pdf(out_file, width = 10.5, height = 9.2)
layout(matrix(1:4, 2, 2, byrow = TRUE))
par(mar = c(4.4, 4.6, 3.2, 1.2), font.main = 2, cex.main = 1.2,
    cex.lab = 1.1, cex.axis = 0.95, mgp = c(2.6, 0.7, 0))

scatter_panel(cell_summary, "Cell Level, Cell Type as Context",
              expression(r[global] ~ "(retained-context cells)"),
              expression(r[resid] ~ "(cell-type means removed)"))
panel_label("A.")

scatter_panel(pb_summary, "Pseudobulk, Cell Type as Context",
              expression(r[global] ~ "(donor" %*% "cell-type units)"),
              expression(r[resid] ~ "(cell-type means removed)"))
panel_label("B.")

donors <- sort(unique(per_donor$donor))
pdC <- per_donor[per_donor$celltype_definition == "cell_type_coarse" &
                   per_donor$n_min == 200 & per_donor$condition_subset %in% c("ctrl", "stim"), ]
pdC <- pdC[order(pdC$donor, pdC$condition_subset), ]
pdC_estimable <- pdC$k_contexts_retained >= 2L
ymax <- max(c(pdC$flip_hi[pdC_estimable], pdC$flip_rate[pdC_estimable], pbmc_ref), na.rm = TRUE) * 1.15
bar_mat <- matrix(ifelse(pdC_estimable, pdC$flip_rate, NA_real_), nrow = 2, dimnames = list(c("ctrl", "stim"), donors))
bp <- barplot(bar_mat, beside = TRUE, col = c(col_ctrl, col_stim), border = NA,
              ylim = c(0, ymax), names.arg = sub("^D", "", donors),
              xlab = "Donor (analyzed alone)", ylab = "Residualized flip rate",
              main = "Single-Donor Flip Rates", las = 1, yaxt = "n")
axis(2, at = pretty(c(0, ymax)), labels = paste0(100 * pretty(c(0, ymax)), "%"), las = 1)
lo <- matrix(ifelse(pdC_estimable, pdC$flip_lo, NA_real_), nrow = 2)
hi <- matrix(ifelse(pdC_estimable, pdC$flip_hi, NA_real_), nrow = 2)
arrows(bp, lo, bp, hi, angle = 90, code = 3, length = 0.02, col = "grey30", lwd = 0.8)
abline(h = pbmc_ref, lty = 2, lwd = 1.2, col = "black")
k_lab <- matrix(pdC$k_contexts_retained, nrow = 2)
lab_y <- ifelse(is.finite(pmax(hi, bar_mat)), pmax(hi, bar_mat) + 0.02 * ymax, 0.02 * ymax)
text(bp, lab_y, labels = ifelse(k_lab >= 2L, k_lab, "N/A"), cex = 0.6, col = "grey30")
legend("topright", bty = "n", cex = 0.85,
       legend = c("control cells", expression("IFN-" * beta * " cells"),
                  sprintf("PBMC 3k, single donor (%.1f%%)", 100 * pbmc_ref)),
       fill = c(col_ctrl, col_stim, NA), border = NA,
       lty = c(NA, NA, 2), lwd = c(NA, NA, 1.2), col = c(NA, NA, "black"))
panel_label("C.")

pdD <- cell_level[cell_level$n_min == 200, ]
def_labels <- list(
  cell_type = "cell type (8)", cell_type_coarse = "cell type, coarse (6)",
  donor = "donor (8)", condition = "condition (2)",
  donor_x_cell_type = expression("donor" %*% "cell type (64)"),
  condition_x_cell_type = expression("condition" %*% "cell type (16)"),
  donor_x_condition = expression("donor" %*% "condition (16)"),
  donor_x_condition_x_cell_type = expression("donor" %*% "condition" %*% "cell type (128)")
)
pdD <- pdD[match(names(def_labels), pdD$context_definition), ]
nD <- nrow(pdD)
ypos <- rev(seq_len(nD))
par(mar = c(4.4, 11.5, 3.2, 1.2))
xmax <- max(c(pdD$flip_hi, pdD$simpson_hi, pdD$mixed_hi), na.rm = TRUE) * 1.08
plot(NA, xlim = c(0, xmax), ylim = c(0.5, nD + 0.5), yaxt = "n", ylab = "", xaxt = "n",
     xlab = "Rate among evaluable pairs", main = "Rates by Context Definition")
axis(1, at = pretty(c(0, xmax)), labels = paste0(100 * pretty(c(0, xmax)), "%"))
abline(h = ypos, col = "grey92")
for (i in seq_len(nD)) axis(2, at = ypos[i], labels = def_labels[[i]], las = 1,
                            cex.axis = 0.85, tick = FALSE)
off <- 0.22
segments(pdD$flip_lo, ypos + off, pdD$flip_hi, ypos + off, col = col_flip, lwd = 1.5)
points(pdD$flip_rate, ypos + off, pch = 16, col = col_flip, cex = 1.1)
segments(pdD$simpson_lo, ypos, pdD$simpson_hi, ypos, col = "#1B4F72", lwd = 1.5)
points(pdD$simpson_rate_abs02, ypos, pch = 17, col = "#1B4F72", cex = 1.1)
segments(pdD$mixed_lo, ypos - off, pdD$mixed_hi, ypos - off, col = "grey40", lwd = 1.5)
points(pdD$mixed_sign_rate, ypos - off, pch = 15, col = "grey40", cex = 1.0)
abline(v = pbmc_ref, lty = 2, col = "black", lwd = 1)
legend("topright", bty = "n", cex = 0.8, pch = c(16, 17, 15, NA),
       lty = c(NA, NA, NA, 2), col = c(col_flip, "#1B4F72", "grey40", "black"),
       legend = c("residualized flip",
                  expression("Simpson reversal, |" * r[global] * "|" >= 0.2),
                  expression("mixed sign (" * epsilon == 0.05 * ")"),
                  "PBMC 3k flip rate"))
panel_label("D.")
invisible(dev.off())

cat(sprintf("Wrote %s with PBMC 3k reference %.6f (%.1f%%).\n",
            out_file, pbmc_ref, 100 * pbmc_ref))
