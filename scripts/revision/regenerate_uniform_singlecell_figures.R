#!/usr/bin/env Rscript

# Regenerate the PBMC/CITE-seq global-versus-residualized figures under the
# revision's uniform population convention: r_global and r_resid are both
# computed on cells belonging to contexts with at least 200 cells.

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(ggplot2)
})
source("00-paths.R")

out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

get_counts <- function(seu, assay) SeuratObject::LayerData(seu[[assay]], layer = "counts")
qc_keep <- function(cnt, min_features = 200L, max_mt = 20) {
  nfeat <- Matrix::colSums(cnt > 0)
  mt <- grepl("^MT-", rownames(cnt))
  pmt <- 100 * Matrix::colSums(cnt[mt, , drop = FALSE]) / Matrix::colSums(cnt)
  nfeat >= min_features & pmt <= max_mt
}
lognorm <- function(cnt) {
  cs <- Matrix::colSums(cnt)
  out <- cnt
  out@x <- log1p(out@x / rep.int(cs, diff(out@p)) * 1e4)
  out
}
clr_cell <- function(x) log1p(x / exp(sum(log1p(x[x > 0]), na.rm = TRUE) / length(x)))
rowcor_pairs <- function(X, Y, ix, iy, chunk = 20000L) {
  zx <- suppressWarnings(t(scale(t(X))))
  zy <- if (identical(X, Y)) zx else suppressWarnings(t(scale(t(Y))))
  denom <- ncol(X) - 1
  out <- numeric(length(ix))
  for (a in seq(1L, length(ix), by = chunk)) {
    b <- min(length(ix), a + chunk - 1L)
    ii <- a:b
    out[ii] <- rowSums(zx[ix[ii], , drop = FALSE] * zy[iy[ii], , drop = FALSE]) / denom
  }
  out[!is.finite(out)] <- NA_real_
  out
}

message("Recomputing PBMC pooled correlations on retained contexts")
pbmc <- readRDS(file.path(paths$scratch, "sc_inputs", "pbmc_scrna.rds"))
cnt <- get_counts(pbmc, "RNA")
keep <- qc_keep(cnt)
ln <- lognorm(cnt[, keep, drop = FALSE])
ar <- readRDS(file.path("results", "sc", "sc_100_pbmc_coarse_celltypes", "analysis_ready.rds"))
ctx <- ar$context
tab <- table(ctx)
retained_pbmc <- names(tab)[tab >= 200]
pop <- ctx %in% retained_pbmc
E <- as.matrix(ln[ar$inputs$feature_x, pop, drop = FALSE])
rg_pbmc <- rowcor_pairs(E, E, ar$pairs$idx_x, ar$pairs$idx_y)
rr_pbmc <- ar$summary_pairs$r_resid
R_pbmc <- ar$r_by_context[, retained_pbmc, drop = FALSE]
rm(pbmc, cnt, ln, E)

message("Recomputing CITE-seq pooled correlations on retained contexts")
cite <- readRDS(file.path(paths$scratch, "sc_inputs", "pbmc_citeseq.rds"))
rna <- get_counts(cite, "RNA")
adt <- get_counts(cite, "ADT")
keep <- qc_keep(rna)
lnr <- lognorm(rna[, keep, drop = FALSE])
adt_clr <- apply(as.matrix(adt[, keep, drop = FALSE]), 2, clr_cell)
dimnames(adt_clr) <- dimnames(adt[, keep, drop = FALSE])
ac <- readRDS(file.path("results", "sc", "sc_200_citeseq_protein_rna_contexts", "analysis_ready.rds"))
ctx <- ac$context_defs$Joint_RNA_ADT
st <- ac$per_definition_stats$Joint_RNA_ADT
tab <- table(ctx)
retained_cite <- st$context_levels[tab[st$context_levels] >= 200]
pop <- ctx %in% retained_cite
af <- ac$inputs$feature_sets$adt_features
rf <- ac$inputs$feature_sets$rna_features
EA <- adt_clr[af, pop, drop = FALSE]
ER <- as.matrix(lnr[rf, pop, drop = FALSE])
rg_cite <- rowcor_pairs(EA, ER, ac$pairs$idx_x, ac$pairs$idx_y)
rr_cite <- st$summary_pairs$r_resid
R_cite <- st$r_by_context[, retained_cite, drop = FALSE]
rm(cite, rna, adt, lnr, adt_clr, EA, ER)

make_df <- function(rg, rr, R) {
  k_valid <- rowSums(is.finite(R) & abs(R) < 1)
  ok <- is.finite(rg) & is.finite(rr) & rg != 0 & rr != 0 & k_valid >= 2
  d <- data.frame(rg = rg[ok], rr = rr[ok])
  d$flip <- sign(d$rg) != sign(d$rr)
  d
}
d_pbmc <- make_df(rg_pbmc, rr_pbmc, R_pbmc)
d_cite <- make_df(rg_cite, rr_cite, R_cite)
rate_pbmc <- mean(d_pbmc$flip)
rate_cite <- mean(d_cite$flip)
message(sprintf("Computed uniform flip rates: PBMC %.12f%%; CITE-seq %.12f%%",
                100 * rate_pbmc, 100 * rate_cite))
stopifnot(abs(100 * rate_pbmc - 11.8366445694557) < 1e-8,
          abs(100 * rate_cite - 37.1115049073064) < 1e-8)

draw_panel <- function(d, title, xlim, ylim, panel = NULL) {
  flip_col <- "#E67E22"
  keep_col <- "#7F8C8D"
  rate <- mean(d$flip)
  ord <- order(d$flip)
  d <- d[ord, ]
  plot(d$rg, d$rr, pch = 16, cex = 0.3,
       col = ifelse(d$flip, adjustcolor(flip_col, alpha.f = 0.4),
                   adjustcolor(keep_col, alpha.f = 0.15)),
       xlab = expression(r[global]), ylab = expression(r[residualized]),
       main = title, xlim = xlim, ylim = ylim)
  grid(col = "grey92", lty = 1)
  abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)
  abline(h = 0, col = "grey80")
  abline(v = 0, col = "grey80")
  legend("topleft", legend = c(sprintf("sign flip: %.1f%%", 100 * rate), "no flip"),
         pch = 16, col = c(flip_col, keep_col), cex = 0.75, bty = "n")
  if (!is.null(panel)) mtext(panel, side = 3, line = 1.0, adj = 0, font = 2, cex = 1.25)
}

pdf_out <- file.path(out_dir, "Figure_SingleCell_GlobalVsResid.pdf")
png_out <- file.path(out_dir, "Figure_SingleCell_GlobalVsResid.png")
grDevices::cairo_pdf(pdf_out, width = 11, height = 4.5)
par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))
draw_panel(d_pbmc, "PBMC scRNA (Coarse Cell Types)", c(-0.7, 0.8), c(-0.5, 0.6), "A")
draw_panel(d_cite, "CITE-seq ADT-RNA (Joint Context)", c(-0.8, 1.0), c(-0.4, 0.7), "B")
invisible(dev.off())
grDevices::png(png_out, width = 11, height = 4.5, units = "in", res = 300,
               type = "cairo-png")
par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))
draw_panel(d_pbmc, "PBMC scRNA (Coarse Cell Types)", c(-0.7, 0.8), c(-0.5, 0.6), "A")
draw_panel(d_cite, "CITE-seq ADT-RNA (Joint Context)", c(-0.8, 1.0), c(-0.4, 0.7), "B")
invisible(dev.off())

df <- d_cite[order(d_cite$flip), ]
df$label <- factor(ifelse(df$flip, sprintf("sign flip: %.1f%%", 100 * rate_cite), "no flip"),
                   levels = c(sprintf("sign flip: %.1f%%", 100 * rate_cite), "no flip"))
p <- ggplot(df, aes(rg, rr, color = label)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_point(size = 0.45, alpha = 0.5) +
  scale_color_manual(values = setNames(c("#E67E22", "#7F8C8D"), levels(df$label))) +
  labs(title = "CITE-seq: Protein vs RNA", x = expression(italic(r)[global]),
       y = expression(italic(r)[resid])) +
  theme_bw(base_size = 10) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.title = element_blank(), legend.key = element_blank(),
        legend.background = element_blank(), legend.text = element_text(size = 7.5),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0, "pt"),
        axis.title = element_text(size = 11), axis.text = element_text(size = 9),
        plot.title = element_text(size = 10.5, face = "bold"),
        panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "Figure_CITEseq_GlobalVsResid_panel.pdf"), p,
       width = 2.7, height = 2.7, device = cairo_pdf)

cat(sprintf("PBMC: %d cells, %.6f%% flips. CITE-seq: %d cells, %.6f%% flips.\n",
            sum(ar$context %in% names(table(ar$context))[table(ar$context) >= 200]),
            100 * rate_pbmc, sum(ac$context_defs$Joint_RNA_ADT %in% retained_cite), 100 * rate_cite))
