#!/usr/bin/env Rscript

# Genome Biology Brief Report -- Figure 2D.
#
# CD14 protein (ADT) vs VCAN (RNA): how the RESIDUALIZED within-context
# correlation depends on the clustering used to define "context". Built
# directly from the cached pair statistics so the plotted residualized
# values exactly match the cache:
#   results/cache/sc_200_citeseq_protein_rna_contexts/pair_stats_*.rds
#     - r_global is identical across definitions (0.92)
#     - r_resid: -0.03 (RNA-defined), 0.50 (ADT-defined), -0.03 (joint WNN)
#
# Design fixes (vs. the old baked image):
#   * single clean 3-facet comparison (RNA / ADT / joint), one title each,
#     NO nested A/B/C panel letters (this sits as panel D of Figure 2).
#   * restrained, perceptually-ordered palette (grey bars + a single accent),
#     not the previous rainbow per-cluster palette.
#
# Run from repo root: Rscript scripts/gb_fig2D_citeseq_contextdef.R

source("00-paths.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cache_dir <- file.path("results", "cache", "sc_200_citeseq_protein_rna_contexts")

fx <- "CD14-TotalSeqB"; fy <- "VCAN"
defs <- c(RNA_defined = "RNA-defined",
          ADT_defined = "ADT-defined",
          Joint_RNA_ADT = "Joint (WNN)")

rows <- list()
summ <- list()
for (cd in names(defs)) {
  ps <- readRDS(file.path(cache_dir, paste0("pair_stats_", cd, ".rds")))
  sp <- ps$summary_pairs
  idx <- which(sp$feature_x == fx & sp$feature_y == fy)
  stopifnot(length(idx) == 1)
  rr <- ps$r_by_context[idx, ]
  nn <- ps$n_by_context[idx, ]
  ok <- is.finite(rr)
  rr <- rr[ok]; nn <- nn[ok]
  ord <- order(rr, decreasing = TRUE)
  rr <- rr[ord]; nn <- nn[ord]
  cl <- factor(seq_along(rr))   # anonymous, ordered cluster index
  rows[[cd]] <- data.frame(
    def = factor(defs[[cd]], levels = unname(defs)),
    cluster = cl, r = as.numeric(rr), n = as.integer(nn))
  summ[[cd]] <- data.frame(
    def = factor(defs[[cd]], levels = unname(defs)),
    r_global = sp$r_global[idx], r_resid = sp$r_resid[idx])
  message(sprintf("%s: r_global=%.3f r_resid=%.3f (%d clusters)",
                  defs[[cd]], sp$r_global[idx], sp$r_resid[idx], length(rr)))
}
df <- do.call(rbind, rows)
sm <- do.call(rbind, summ)

# Facet titles carry the residualized value (the quantity that changes).
sm$lab <- sprintf("%s\nresidualized r = %+.2f", sm$def, sm$r_resid)
lab_map <- setNames(sm$lab, as.character(sm$def))
df$facet <- factor(lab_map[as.character(df$def)], levels = lab_map[as.character(sm$def)])

# Sign-of-correlation palette: bars colored by whether the within-cluster r is
# positive or negative. This makes the story visible -- ADT-defined clusters are
# consistently positive (a real within-cluster correlation), while RNA-defined
# and joint clusters sit near zero (compositional).
pos_col <- "#1F77B4"  # blue   -> within-cluster r > 0
neg_col <- "#E8601C"  # orange -> within-cluster r < 0
accent <- "#C0392B"   # global reference line
df$sign <- factor(ifelse(df$r >= 0, "within-cluster r > 0", "within-cluster r < 0"),
                  levels = c("within-cluster r > 0", "within-cluster r < 0"))

p <- ggplot(df, aes(x = cluster, y = r, fill = sign)) +
  geom_col(width = 0.95) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  geom_hline(data = sm, aes(yintercept = r_global),
             color = accent, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ facet, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c("within-cluster r > 0" = pos_col,
                               "within-cluster r < 0" = neg_col),
                    name = NULL, drop = FALSE) +
  labs(x = "Cluster (ordered by within-cluster r)",
       y = "Within-cluster correlation",
       title = sprintf("%s vs %s: global r = %.2f (dashed); context defines the residual",
                       "CD14 protein", fy, sm$r_global[1])) +
  coord_cartesian(ylim = c(-0.15, 1)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0.4, "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 9.5),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-6, 0, 0, 0),
        legend.text = element_text(size = 8),
        strip.text = element_text(size = 9, lineheight = 0.95),
        strip.background = element_rect(fill = "grey92", color = NA))

W_in <- 180 / 25.4
H_in <- 70 / 25.4
ggsave(file.path(out_dir, "Figure_SingleCell_CITEseq_ContextDefComparison.pdf"),
       p, width = W_in, height = H_in, units = "in", device = cairo_pdf)
message("Wrote Figure_SingleCell_CITEseq_ContextDefComparison.pdf")
