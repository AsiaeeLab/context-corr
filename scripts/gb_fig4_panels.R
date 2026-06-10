#!/usr/bin/env Rscript

# Genome Biology Brief Report -- Figure 4 single panels (no baked-in letters).
#
# Produces four clean single-panel PDFs in doc/genomebiology/figures/:
#   gb_fig4B_brca_mixsign.pdf  Sign-mixedness rate vs effect-size threshold (PAM50 vs IHC)
#   gb_fig4C_brca_i2.pdf       I^2 distribution by context (PAM50 vs IHC), non-overlapping medians
#   gb_fig4D_brca_pam50.pdf    PAM50 exemplar EN1 x hsa-miR-577 scatter
#   gb_fig4E_brca_ihc.pdf      IHC exemplar SLFN11 x hsa-miR-99a-5p scatter
#
# Per-context r values are placed in the legend (clean) rather than a cramped subtitle.
#
# Run from repo root: Rscript scripts/gb_fig4_panels.R

source("00-paths.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pal_pam <- c(Basal = "#E41A1C", Her2 = "#377EB8", LumA = "#4DAF4A",
             LumB = "#984EA3", Normal = "#FF7F00")
pal_ihc <- c(`HER2+` = "#377EB8", `HR+` = "#4DAF4A", TNBC = "#E41A1C")
pal_context <- c(PAM50 = "#1B9E77", IHC = "#D95F02")

theme_panel <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 9),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5),
        legend.margin = margin(0, 0, 0, 0))

res_pam50 <- readRDS(file.path("results", "cache", "brca_subtype_PAM50.rds"))
res_ihc   <- readRDS(file.path("results", "cache", "brca_subtype_IHC.rds"))
rates     <- read.csv(file.path("results", "brca_subtype_rates.csv"), stringsAsFactors = FALSE)
pl_pam    <- readRDS(file.path("results", "cache", "brca_top_reversal_PAM50.rds"))
pl_ihc    <- readRDS(file.path("results", "cache", "brca_top_reversal_IHC.rds"))

W_in <- 85 / 25.4
H_in <- 60 / 25.4

# ---------------------------------------------------------------------------
# Panel B: sign-mixedness rate vs effect-size threshold
# ---------------------------------------------------------------------------
ratesA <- subset(rates, eps == 0.05 & threshold %in% c(0, 0.05, 0.1, 0.15, 0.2))
if (!any(ratesA$threshold == 0.15)) {
  ratesA <- subset(rates, eps == 0.05 & threshold %in% c(0, 0.05, 0.1, 0.2, 0.3))
}
ratesA$context <- factor(ratesA$context, levels = c("PAM50", "IHC"))

pB <- ggplot(ratesA, aes(x = threshold, y = 100 * frac_mixed,
                         color = context, group = context)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_color_manual(values = pal_context, name = "Context") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title = "Sign-mixedness vs effect-size threshold",
       x = expression(paste("Threshold on ", max[c]*"|"*r[c]*"|")),
       y = "% mixed-sign pairs") +
  theme_panel
ggsave(file.path(out_dir, "gb_fig4B_brca_mixsign.pdf"), pB,
       width = W_in, height = H_in, units = "in", device = cairo_pdf)

# ---------------------------------------------------------------------------
# Panel C: I^2 distribution; medians stacked in a legend-style annotation box
# ---------------------------------------------------------------------------
make_i2_df <- function(res, label) {
  v <- as.vector(res$I2); ok <- is.finite(v)
  data.frame(context = label, I2 = v[ok], stringsAsFactors = FALSE)
}
i2_df <- rbind(make_i2_df(res_pam50, "PAM50"), make_i2_df(res_ihc, "IHC"))
set.seed(1)
if (nrow(i2_df) > 80000) i2_df <- i2_df[sample(nrow(i2_df), 80000), ]
i2_df$context <- factor(i2_df$context, levels = c("PAM50", "IHC"))

med_pam <- median(i2_df$I2[i2_df$context == "PAM50"], na.rm = TRUE)
med_ihc <- median(i2_df$I2[i2_df$context == "IHC"], na.rm = TRUE)

# Single non-overlapping annotation in the upper area, two stacked lines.
ann <- data.frame(
  lab = c(sprintf("median PAM50 = %.2f", med_pam),
          sprintf("median IHC = %.2f", med_ihc)),
  col = c(pal_context["PAM50"], pal_context["IHC"]),
  yrel = c(0.96, 0.86))

pC <- ggplot(i2_df, aes(x = I2, fill = context, color = context)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_vline(xintercept = med_pam, color = pal_context["PAM50"],
             linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = med_ihc, color = pal_context["IHC"],
             linetype = "dashed", linewidth = 0.5) +
  annotate("label", x = 0.97, y = c(Inf, Inf),
           label = c(sprintf("median PAM50 = %.2f", med_pam),
                     sprintf("median IHC = %.2f", med_ihc)),
           hjust = 1, vjust = c(1.4, 3.0), size = 2.3,
           color = c(pal_context["PAM50"], pal_context["IHC"]),
           label.size = 0, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = pal_context, name = "Context") +
  scale_color_manual(values = pal_context, name = "Context") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = expression(paste("Distribution of ", I^2, " by context")),
       x = expression(I^2), y = "Density") +
  theme_panel
ggsave(file.path(out_dir, "gb_fig4C_brca_i2.pdf"), pC,
       width = W_in, height = H_in, units = "in", device = cairo_pdf)

# ---------------------------------------------------------------------------
# Panels D and E: per-pair scatters; per-context r in the legend labels
# ---------------------------------------------------------------------------
make_scatter <- function(pl, palette, ctx_levels, title_main, ctx_label) {
  d <- pl$data
  d$context <- factor(d$context, levels = ctx_levels)
  d <- d[!is.na(d$context), ]
  r_global <- pl$summary$r_global

  per_ctx <- do.call(rbind, lapply(levels(d$context), function(L) {
    sub <- d[d$context == L, ]
    r <- if (nrow(sub) < 3) NA_real_ else suppressWarnings(stats::cor(sub$mRNA, sub$miR))
    data.frame(context = L, r = r, stringsAsFactors = FALSE)
  }))
  # Legend labels carry the per-context r value.
  lab_map <- setNames(
    sprintf("%s (r = %+.2f)", per_ctx$context, per_ctx$r),
    per_ctx$context)
  d$context_lab <- factor(lab_map[as.character(d$context)],
                          levels = lab_map[ctx_levels])
  pal_lab <- setNames(palette[ctx_levels], lab_map[ctx_levels])

  ggplot(d, aes(x = mRNA, y = miR, color = context_lab)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "grey40", linetype = "dashed",
                linewidth = 0.6, formula = y ~ x) +
    geom_smooth(aes(group = context_lab), method = "lm", se = FALSE,
                linewidth = 0.7, formula = y ~ x) +
    scale_color_manual(values = pal_lab,
                       name = sprintf("%s  (global r = %+.2f)", ctx_label, r_global),
                       drop = FALSE) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.6), ncol = 1,
                                title.position = "top")) +
    labs(title = title_main,
         x = sprintf("%s (log2 expression)", pl$gene),
         y = sprintf("%s (log2 RPM)", pl$mir)) +
    theme_panel +
    theme(legend.title = element_text(size = 7))
}

pD <- make_scatter(pl_pam, pal_pam,
                   ctx_levels = c("Basal", "Her2", "LumA", "LumB", "Normal"),
                   title_main = sprintf("PAM50 exemplar: %s x %s", pl_pam$gene, pl_pam$mir),
                   ctx_label = "PAM50")
ggsave(file.path(out_dir, "gb_fig4D_brca_pam50.pdf"), pD,
       width = W_in, height = H_in, units = "in", device = cairo_pdf)

pE <- make_scatter(pl_ihc, pal_ihc,
                   ctx_levels = c("HER2+", "HR+", "TNBC"),
                   title_main = sprintf("IHC exemplar: %s x %s", pl_ihc$gene, pl_ihc$mir),
                   ctx_label = "IHC")
ggsave(file.path(out_dir, "gb_fig4E_brca_ihc.pdf"), pE,
       width = W_in, height = H_in, units = "in", device = cairo_pdf)

cat("\n=== SUMMARY ===\n")
cat(sprintf("median I2 PAM50 = %.3f ; IHC = %.3f\n", med_pam, med_ihc))
cat(sprintf("PAM50 exemplar %s x %s: r_global = %+.3f\n",
            pl_pam$gene, pl_pam$mir, pl_pam$summary$r_global))
cat(sprintf("IHC exemplar %s x %s: r_global = %+.3f\n",
            pl_ihc$gene, pl_ihc$mir, pl_ihc$summary$r_global))
cat("Wrote gb_fig4B/C/D/E.\n")
