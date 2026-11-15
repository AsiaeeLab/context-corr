#!/usr/bin/env Rscript

# Regenerate the downstream-network figure from its result tables, without
# rerunning WGCNA or GENIE3. This keeps the plotted sample-unit label aligned
# with the manuscript (8,890 TCGA samples).

out_dir <- file.path("results", "revision", "downstream")
edge_class <- read.csv(file.path(out_dir, "edge_classification.csv"), stringsAsFactors = FALSE)
benchmark <- read.csv(file.path(out_dir, "benchmark_precision.csv"), stringsAsFactors = FALSE)
edge_overlap <- read.csv(file.path(out_dir, "edge_overlap.csv"), stringsAsFactors = FALSE)
module_agreement <- read.csv(file.path(out_dir, "module_agreement.csv"), stringsAsFactors = FALSE)

fnum <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
k_top <- c(1000L, 5000L, 10000L)
N_total <- 8890L
class_levels <- c("globally_robust", "compositional", "context_switching", "indeterminate")
class_cols <- c(globally_robust = "#1B4F72", compositional = "#E67E22",
                context_switching = "#C0392B", indeterminate = "grey72")
input_cols <- c(pooled = "#E67E22", residualized = "#1B4F72", RE_pooled = "#7F8C8D")
panel_label <- function(lab) mtext(lab, side = 3, line = 1.2, adj = 0,
                                    font = 2, cex = 1.3)

pdf_file <- file.path(out_dir, "Figure_Downstream.pdf")
grDevices::cairo_pdf(pdf_file, width = 10.0, height = 8.4)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.95))
par(mar = c(4.6, 4.8, 3.4, 1.2), oma = c(0, 0.6, 0, 0), mgp = c(2.8, 0.7, 0))

sc <- edge_class[edge_class$in_random_subsample, ]
sc <- sc[order(factor(sc$label,
                      levels = c("indeterminate", "globally_robust", "compositional", "context_switching"))), ]
lim <- range(c(sc$r_pooled, sc$r_resid), na.rm = TRUE)
lim <- c(-1, 1) * max(abs(lim))
plot(sc$r_pooled, sc$r_resid, pch = 16, cex = 0.45,
     col = grDevices::adjustcolor(class_cols[sc$label], alpha.f = 0.55),
     xlim = lim, ylim = lim,
     xlab = bquote(r[pooled] ~ .(sprintf("(all %s samples)", fnum(N_total)))),
     ylab = expression(r[resid] ~ "(cohort means removed)"), main = "", las = 1)
abline(h = 0, v = 0, col = "grey60", lty = 3)
abline(0, 1, col = "grey40", lty = 2)
title(main = "Pooled vs Residualized Edge Weights", font.main = 2, cex.main = 1.1)
panel_label("A.")
cnt <- table(factor(sc$label, levels = class_levels))
legend("topleft", bty = "n", cex = 0.85, pch = 16, pt.cex = 1.1,
       col = class_cols[class_levels],
       legend = sprintf("%s (n = %s)",
                        c("Globally robust", "Compositional", "Context switching", "Indeterminate"),
                        fnum(as.integer(cnt))))

prec_curve <- function(rn) {
  b <- benchmark[benchmark$ranking == rn & benchmark$subset == "top_k" &
                   benchmark$edge_set == "gene_pairs", ]
  b[match(k_top, b$k), "precision_mappable"]
}
curves <- list(
  pooled = prec_curve("pooled_abs_r"), residualized = prec_curve("residualized_abs_r"),
  RE_pooled = prec_curve("RE_pooled_abs_r"), pooled_TOM = prec_curve("pooled_TOM"),
  residualized_TOM = prec_curve("residualized_TOM"), RE_pooled_TOM = prec_curve("RE_pooled_TOM")
)
base_row <- benchmark[benchmark$ranking == "random_pairs_exact" &
                        benchmark$subset == "all_mappable_pairs" &
                        benchmark$edge_set == "gene_pairs", ]
stopifnot(nrow(base_row) == 1L)
baseline_string <- base_row$precision_mappable
ymax <- max(unlist(curves), baseline_string, na.rm = TRUE) * 1.15
plot(NA, xlim = c(800, 12500), ylim = c(0, ymax), log = "x", xaxt = "n", las = 1,
     xlab = expression("Top-" * italic(k) * " edges by weight"),
     ylab = "Precision (fraction with STRING physical link)", main = "")
axis(1, at = k_top, labels = fnum(k_top))
abline(h = baseline_string, col = "grey35", lty = 3, lwd = 1.6)
for (nm in c("pooled", "residualized", "RE_pooled")) {
  lines(k_top, curves[[nm]], col = input_cols[nm], lwd = 2.4)
  points(k_top, curves[[nm]], col = input_cols[nm], pch = 16, cex = 1.2)
  lines(k_top, curves[[paste0(nm, "_TOM")]], col = input_cols[nm], lwd = 2, lty = 2)
  points(k_top, curves[[paste0(nm, "_TOM")]], col = input_cols[nm], pch = 17, cex = 1.1)
}
title(main = "STRING Support of Top-k Edges", font.main = 2, cex.main = 1.1)
panel_label("B.")
legend("topright", bty = "n", cex = 0.82, ncol = 1,
       legend = c("Pooled", "Residualized", expression("RE-pooled " * r[RE]),
                  expression("Ranking by |" * r * "|"), "Ranking by WGCNA TOM", "Random gene pairs"),
       col = c(input_cols, "grey20", "grey20", "grey35"),
       lwd = c(2.4, 2.4, 2.4, 2.4, 2, 1.6), lty = c(1, 1, 1, 1, 2, 3),
       pch = c(16, 16, 16, 16, 17, NA), pt.cex = c(1.2, 1.2, 1.2, 1.2, 1.1, NA))

get_ov <- function(level, b, col) {
  r <- edge_overlap[edge_overlap$level == level & edge_overlap$input_a == "pooled" &
                      edge_overlap$input_b == b, ]
  if (nrow(r) == 0) NA_real_ else r[[col]]
}
get_ari <- function(b) module_agreement$ari_all_genes_unassigned_as_module[
  module_agreement$input_a == "pooled" & module_agreement$input_b == b]
bar_mat <- rbind(
  residualized = c(get_ov("similarity_r", "residualized", "jaccard_top1000"),
                   get_ov("similarity_r", "residualized", "jaccard_top10000"),
                   get_ov("WGCNA_TOM", "residualized", "jaccard_top1000"),
                   get_ov("WGCNA_TOM", "residualized", "jaccard_top10000"),
                   get_ari("residualized"), get_ov("GENIE3_weight", "residualized", "jaccard_top1000"),
                   get_ov("GENIE3_weight", "residualized", "jaccard_top10000")),
  RE_pooled = c(get_ov("similarity_r", "RE_pooled", "jaccard_top1000"),
                get_ov("similarity_r", "RE_pooled", "jaccard_top10000"),
                get_ov("WGCNA_TOM", "RE_pooled", "jaccard_top1000"),
                get_ov("WGCNA_TOM", "RE_pooled", "jaccard_top10000"),
                get_ari("RE_pooled"), NA, NA)
)
colnames(bar_mat) <- c("|r|\ntop 1,000", "|r|\ntop 10,000", "TOM\ntop 1,000",
                       "TOM\ntop 10,000", "Modules\n(ARI)", "GENIE3\ntop 1,000", "GENIE3\ntop 10,000")
par(mar = c(5.2, 4.8, 3.4, 1.2))
bp_x <- barplot(bar_mat, beside = TRUE,
                col = input_cols[c("residualized", "RE_pooled")], border = NA,
                ylim = c(0, 1), las = 1, names.arg = rep("", ncol(bar_mat)),
                ylab = "Agreement with the pooled output", main = "")
axis(1, at = colMeans(bp_x), labels = FALSE, tick = FALSE)
mtext(colnames(bar_mat), side = 1, at = colMeans(bp_x), line = 2.4, cex = 0.8)
title(main = "Agreement of Pooled and Context-Aware Outputs", font.main = 2, cex.main = 1.1)
panel_label("C.")
ceil <- c(get_ov("GENIE3_weight_subset", "pooled_second_seed", "jaccard_top1000"),
          get_ov("GENIE3_weight_subset", "pooled_second_seed", "jaccard_top10000"))
if (all(is.finite(ceil))) {
  for (m in 1:2) {
    xg <- bp_x[, ncol(bar_mat) - 2 + m]
    segments(min(xg) - 0.45, ceil[m], max(xg) + 0.45, ceil[m],
             lwd = 2.2, lty = 1, col = "grey20")
  }
}
legend("topleft", bty = "n", cex = 0.85, ncol = 1, inset = c(0.01, 0),
       fill = c(input_cols[c("residualized", "RE_pooled")], NA), border = NA,
       legend = c("Pooled vs residualized", expression("Pooled vs RE-pooled " * r[RE]),
                  "GENIE3 same-input seed replicate"),
       lty = c(NA, NA, 1), lwd = c(NA, NA, 2.2), col = c(NA, NA, "grey20"), merge = FALSE)
invisible(dev.off())

cat(sprintf("Wrote %s with sample-unit label '%s'.\n", pdf_file,
            sprintf("all %s samples", fnum(N_total))))
