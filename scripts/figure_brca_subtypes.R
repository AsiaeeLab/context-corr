#!/usr/bin/env Rscript

# Supplementary Figure S1 (Genome Biology Brief Report):
# BRCA molecular-subtype analysis -- 4-panel composite.
#
# Panels:
#   A. Sign-mixedness rate vs effect-size threshold (PAM50 vs IHC), eps = 0.05
#   B. I^2 distribution across pairs (PAM50 vs IHC, overlaid densities)
#   C. PAM50 exemplar: EN1 x hsa-miR-577 scatter
#   D. IHC exemplar:   SLFN11 x hsa-miR-99a-5p scatter
#
# Outputs:
#   doc/genomebiology/figures/Figure_BRCA_Subtypes.pdf
#   doc/genomebiology/figures/Figure_BRCA_Subtypes.png
#
# Run from repo root: Rscript scripts/figure_brca_subtypes.R

source("00-paths.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------
out_dir <- file.path("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(out_dir, "Figure_BRCA_Subtypes.pdf")
out_png <- file.path(out_dir, "Figure_BRCA_Subtypes.png")

# ---------------------------------------------------------------------------
# Color palettes (colorblind-friendly Set1-based, matching 420-brca-subtypes.Rmd)
# ---------------------------------------------------------------------------
pal_pam <- c(Basal  = "#E41A1C",  # red
             Her2   = "#377EB8",  # blue
             LumA   = "#4DAF4A",  # green
             LumB   = "#984EA3",  # purple
             Normal = "#FF7F00")  # orange

pal_ihc <- c(`HER2+` = "#377EB8",
             `HR+`   = "#4DAF4A",
             TNBC    = "#E41A1C")

pal_context <- c(PAM50 = "#1B9E77", IHC = "#D95F02")  # Dark2

theme_panel <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7.5),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5))

# ---------------------------------------------------------------------------
# Load cached results
# ---------------------------------------------------------------------------
res_pam50 <- readRDS(file.path("results", "cache", "brca_subtype_PAM50.rds"))
res_ihc   <- readRDS(file.path("results", "cache", "brca_subtype_IHC.rds"))
rates     <- read.csv(file.path("results", "brca_subtype_rates.csv"),
                      stringsAsFactors = FALSE)
het       <- read.csv(file.path("results", "brca_subtype_heterogeneity.csv"),
                      stringsAsFactors = FALSE)

pl_pam <- readRDS(file.path("results", "cache", "brca_top_reversal_PAM50.rds"))
pl_ihc <- readRDS(file.path("results", "cache", "brca_top_reversal_IHC.rds"))

cat("Loaded cached results.\n")
cat("PAM50 exemplar: ", pl_pam$gene, "x", pl_pam$mir,
    "| r_global =", round(pl_pam$summary$r_global, 3), "\n")
cat("IHC   exemplar: ", pl_ihc$gene, "x", pl_ihc$mir,
    "| r_global =", round(pl_ihc$summary$r_global, 3), "\n")

# ---------------------------------------------------------------------------
# Panel A: sign-mixedness rate vs effect-size threshold (eps = 0.05)
# ---------------------------------------------------------------------------
ratesA <- subset(rates, eps == 0.05 &
                        threshold %in% c(0, 0.05, 0.1, 0.15, 0.2))
# threshold = 0.15 may be missing in cache; if so, fall back to existing grid.
if (!any(ratesA$threshold == 0.15)) {
  ratesA <- subset(rates, eps == 0.05 &
                          threshold %in% c(0, 0.05, 0.1, 0.2, 0.3))
}
ratesA$context <- factor(ratesA$context, levels = c("PAM50", "IHC"))

pA <- ggplot(ratesA, aes(x = threshold, y = 100 * frac_mixed,
                         color = context, group = context)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  scale_color_manual(values = pal_context, name = "Context") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title = "A.  Sign-mixedness rate vs effect-size threshold",
       x = expression(paste("Effect-size threshold ", epsilon[t], " on ", max[c]*"|"*r[c]*"|")),
       y = "% mixed-sign pairs") +
  theme_panel

# ---------------------------------------------------------------------------
# Panel B: I^2 distribution (overlaid densities, PAM50 vs IHC)
# ---------------------------------------------------------------------------
make_i2_df <- function(res, label) {
  v <- as.vector(res$I2)
  ok <- is.finite(v)
  data.frame(context = label, I2 = v[ok], stringsAsFactors = FALSE)
}
i2_df <- rbind(make_i2_df(res_pam50, "PAM50"),
               make_i2_df(res_ihc,  "IHC"))
# Subsample for plotting speed
set.seed(1)
if (nrow(i2_df) > 80000) {
  i2_df <- i2_df[sample(nrow(i2_df), 80000), ]
}
i2_df$context <- factor(i2_df$context, levels = c("PAM50", "IHC"))

med_pam <- median(i2_df$I2[i2_df$context == "PAM50"], na.rm = TRUE)
med_ihc <- median(i2_df$I2[i2_df$context == "IHC"],   na.rm = TRUE)

pB <- ggplot(i2_df, aes(x = I2, fill = context, color = context)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  geom_vline(xintercept = med_pam, color = pal_context["PAM50"],
             linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = med_ihc, color = pal_context["IHC"],
             linetype = "dashed", linewidth = 0.5) +
  annotate("text",
           x = med_pam, y = Inf,
           label = sprintf("median PAM50 = %.2f", med_pam),
           vjust = 1.6, hjust = -0.05, size = 2.4,
           color = pal_context["PAM50"]) +
  annotate("text",
           x = med_ihc, y = Inf,
           label = sprintf("median IHC = %.2f", med_ihc),
           vjust = 3.0, hjust = -0.05, size = 2.4,
           color = pal_context["IHC"]) +
  scale_fill_manual(values = pal_context, name = "Context") +
  scale_color_manual(values = pal_context, name = "Context") +
  labs(title = expression(bold("B.  Distribution of ") * bolditalic(I)^bold("2") *
                          bold(" across pairs")),
       x = expression(I^2), y = "Density") +
  theme_panel

# ---------------------------------------------------------------------------
# Panels C and D: per-pair scatters with per-context regression lines
# ---------------------------------------------------------------------------
make_scatter <- function(pl, palette, ctx_levels, panel_letter,
                         title_main, ctx_label) {
  d <- pl$data
  d$context <- factor(d$context, levels = ctx_levels)
  d <- d[!is.na(d$context), ]

  r_global <- pl$summary$r_global

  # Per-context regressions
  per_ctx <- do.call(rbind, lapply(levels(d$context), function(L) {
    sub <- d[d$context == L, ]
    if (nrow(sub) < 3) {
      data.frame(context = L, r = NA_real_, n = nrow(sub))
    } else {
      data.frame(context = L,
                 r = suppressWarnings(stats::cor(sub$mRNA, sub$miR)),
                 n = nrow(sub))
    }
  }))

  # Subtitle annotation
  per_ctx_lab <- paste(
    sprintf("%s: %+.2f", per_ctx$context, per_ctx$r),
    collapse = " | "
  )
  sub_text <- sprintf("%s = %+.2f  |  %s",
                      "r[global]", r_global, per_ctx_lab)
  sub_text_plain <- sprintf("r_global = %+.2f  |  %s",
                            r_global, per_ctx_lab)

  p <- ggplot(d, aes(x = mRNA, y = miR, color = context)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "grey40", linetype = "dashed",
                linewidth = 0.6, formula = y ~ x) +
    geom_smooth(aes(group = context), method = "lm", se = FALSE,
                linewidth = 0.7, formula = y ~ x) +
    scale_color_manual(values = palette, name = ctx_label,
                       drop = FALSE) +
    labs(title = paste0(panel_letter, ".  ", title_main),
         subtitle = sub_text_plain,
         x = paste(pl$gene, "(log2 expression)"),
         y = paste(pl$mir,  "(log2 RPM)")) +
    theme_panel

  p
}

pC <- make_scatter(pl_pam, pal_pam,
                   ctx_levels = c("Basal", "Her2", "LumA", "LumB", "Normal"),
                   panel_letter = "C",
                   title_main = paste("PAM50 exemplar:", pl_pam$gene, "x", pl_pam$mir),
                   ctx_label = "PAM50")

pD <- make_scatter(pl_ihc, pal_ihc,
                   ctx_levels = c("HER2+", "HR+", "TNBC"),
                   panel_letter = "D",
                   title_main = paste("IHC exemplar:", pl_ihc$gene, "x", pl_ihc$mir),
                   ctx_label = "IHC")

# ---------------------------------------------------------------------------
# Compose 2 x 2 layout
# ---------------------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))

# Target: full-page-width ~170 mm x 120 mm
W_in <- 170 / 25.4
H_in <- 120 / 25.4

ggsave(out_pdf, fig, width = W_in, height = H_in, units = "in",
       device = cairo_pdf)
ggsave(out_png, fig, width = W_in, height = H_in, units = "in",
       dpi = 300)

cat("Wrote:\n  ", out_pdf, "\n  ", out_png, "\n")
