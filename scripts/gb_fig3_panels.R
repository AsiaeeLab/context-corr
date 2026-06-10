#!/usr/bin/env Rscript

# Genome Biology Brief Report -- Figure 3 single panels (no baked-in letters).
#
# Produces:
#   doc/genomebiology/figures/gb_fig3A_gtex.pdf   GTEx mRNA-mRNA Simpson example
#                                                 (RTKN2 vs RNASE3): global scatter
#                                                 colored by tissue + global regression line.
#   doc/genomebiology/figures/gb_fig3B_pbmc.pdf   PBMC scRNA LYZ vs FTH1: global scatter
#                                                 colored by cell type with per-cell-type
#                                                 regression lines diverging from the global line.
#
# Run from repo root: Rscript scripts/gb_fig3_panels.R

source("00-paths.R")
source("R/contextcorr.R")
source("R/sc_datasets.R")
source("R/sc_context_pipeline.R")

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- file.path("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_panel <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "right",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8.5),
        legend.margin = margin(0, 0, 0, 0))

# ===========================================================================
# Panel A: GTEx RTKN2 vs RNASE3
# ===========================================================================
make_gtex_panel <- function() {
  f2 <- file.path(paths$mirtcga, "xena-gtex.Rda")
  f1 <- "/home/amir/datasets/mirTCGA/clean/xena-gtex.Rda"
  f <- if (file.exists(f1)) f1 else f2
  message("Loading GTEx: ", f)
  e <- new.env()
  load(f, envir = e)
  stopifnot(is.list(e$gtex), "mrna" %in% names(e$gtex))
  X <- as.matrix(e$gtex$mrna)  # genes (Ensembl) x samples, log2(TPM+0.001)
  message("GTEx matrix: ", nrow(X), " genes x ", ncol(X), " samples")

  # True tissue comes from the GTEx v8 Sample Attributes (SMTSD), joined on
  # sample ID. The 4-digit field in the sample barcode is NOT tissue (it is an
  # aliquot/collection-order number that spans many tissues). Cache a copy of
  # the annotation locally.
  attr_local <- "/home/amir/datasets/context/clean/GTEx_SampleAttributes_v8.txt"
  attr_url <- paste0("https://storage.googleapis.com/adult-gtex/annotations/v8/",
                     "metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
  if (!file.exists(attr_local)) {
    dir.create(dirname(attr_local), recursive = TRUE, showWarnings = FALSE)
    if (file.exists("/tmp/gtex_attr.txt")) {
      file.copy("/tmp/gtex_attr.txt", attr_local)
    } else {
      utils::download.file(attr_url, attr_local, quiet = TRUE)
    }
  }
  attr <- utils::read.delim(attr_local, stringsAsFactors = FALSE, quote = "")
  samp_dash <- gsub("\\.", "-", colnames(X))
  smtsd <- attr$SMTSD[match(samp_dash, attr$SAMPID)]

  tab <- sort(table(smtsd[!is.na(smtsd)]), decreasing = TRUE)
  min_n <- 150
  keep_tissues <- names(tab)[tab >= min_n]
  keep_tissues <- keep_tissues[seq_len(min(20, length(keep_tissues)))]

  # Rows are versioned Ensembl IDs; match by base ID.
  base_id <- sub("\\..*$", "", rownames(X))
  # CEACAM3 vs EN1: a clean Simpson example under TRUE tissue grouping --
  # strongly negative global correlation, positive within-tissue in a 15/20
  # majority, so mean-residualization flips the sign.
  g1 <- "CEACAM3"; g2 <- "EN1"
  ens1 <- "ENSG00000170956"; ens2 <- "ENSG00000163064"
  i1 <- which(base_id == ens1); i2 <- which(base_id == ens2)
  stopifnot(length(i1) == 1, length(i2) == 1)

  x <- as.numeric(X[i1, ])
  y <- as.numeric(X[i2, ])
  tt <- smtsd
  in_keep <- tt %in% keep_tissues & is.finite(x) & is.finite(y)
  x <- x[in_keep]; y <- y[in_keep]; tt <- factor(tt[in_keep], levels = keep_tissues)
  tt <- droplevels(tt)

  r_global <- suppressWarnings(stats::cor(x, y))
  # per-tissue r
  r_by <- tapply(seq_along(x), tt, function(idx) suppressWarnings(stats::cor(x[idx], y[idx])))
  message(sprintf("GTEx %s vs %s: global r = %.3f", g1, g2, r_global))
  message(sprintf("  per-tissue r range: %.2f .. %.2f ; positive in %d/%d tissues",
                  min(r_by, na.rm = TRUE), max(r_by, na.rm = TRUE),
                  sum(r_by > 0, na.rm = TRUE), sum(is.finite(r_by))))

  df <- data.frame(x = x, y = y, tissue = tt)
  # Non-highlighted tissues stay muted grey; the two highlighted tissues get
  # clearly distinct accent colors (blue vs orange/red) so they read apart.
  cols <- setNames(rep("grey75", nlevels(df$tissue)), levels(df$tissue))

  # Annotate the two tissues whose within-tissue line diverges most from the
  # global line, BY NAME (SMTSD): the strongest opposite-sign (positive) tissue
  # and the most negative tissue.
  ord <- order(r_by, na.last = NA)
  most_neg <- names(r_by)[ord[1]]                       # most negative tissue
  most_pos <- names(r_by)[ord[length(ord)]]             # strongest positive (opposite to global)
  hi_pos_col <- "#1F77B4"  # blue   -> strongest positive (e.g. Thyroid)
  hi_neg_col <- "#E8601C"  # orange -> most negative (e.g. Esophagus - Mucosa)
  cols[most_pos] <- hi_pos_col
  cols[most_neg] <- hi_neg_col
  ann_codes <- c(most_pos, most_neg)

  # Position the strongest-positive (Thyroid) label in the upper-right of the
  # plot area so it is not lost mid-cloud; keep the most-negative label near its
  # own tissue cloud where it is already legible.
  x_rng <- range(df$x); y_rng <- range(df$y)
  pos_xy <- c(x = x_rng[1] + 0.92 * diff(x_rng),
              y = y_rng[1] + 0.95 * diff(y_rng))
  ann_df <- do.call(rbind, lapply(ann_codes, function(cc) {
    sel <- df$tissue == cc
    if (cc == most_pos) {
      ax <- pos_xy["x"]; ay <- pos_xy["y"]
    } else {
      ax <- stats::median(df$x[sel]); ay <- stats::median(df$y[sel])
    }
    data.frame(x = ax, y = ay, tissue = cc,
               lab = sprintf("%s (r = %+.2f)", cc, r_by[cc]),
               stringsAsFactors = FALSE)
  }))

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = tissue), size = 0.5, alpha = 0.40) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "black", linewidth = 0.9, formula = y ~ x) +
    geom_smooth(data = df[df$tissue %in% ann_codes, ],
                aes(group = tissue, color = tissue), method = "lm", se = FALSE,
                linewidth = 0.9, formula = y ~ x) +
    scale_color_manual(values = cols, guide = "none") +
    ggrepel::geom_label_repel(data = ann_df, aes(x = x, y = y, label = lab, color = tissue),
                              size = 2.6, label.size = 0.25, fill = "white",
                              alpha = 0.95, min.segment.length = 0, seed = 1,
                              show.legend = FALSE) +
    labs(title = sprintf("GTEx: %s vs %s (global r = %+.2f)", g1, g2, r_global),
         x = sprintf("%s (log2 expression)", g1),
         y = sprintf("%s (log2 expression)", g2)) +
    theme_panel + theme(legend.position = "none")
  list(plot = p, r_global = r_global, r_by = r_by)
}

# ===========================================================================
# Panel B: PBMC LYZ vs FTH1 -- compositional view
# ===========================================================================
make_pbmc_panel <- function() {
  fx <- "LYZ"; fy <- "FTH1"
  pbmc <- load_pbmc_scrna(paths = paths, auto_download = FALSE, verbose = FALSE)
  seu <- pbmc$seurat
  seu <- sc_preprocess(seu, assay = "RNA", min_features = 200, max_mt = 20,
                       n_hvg = 2000, normalize_mode = "log", residualize = FALSE,
                       vars_to_regress = NULL, seed = 100, verbose = FALSE)
  seu <- sc_define_coarse_context(seu, context_col = "context_coarse",
                                  prefer_metadata = TRUE, unknown_label = "Other")
  ctx <- as.character(seu$context_coarse)
  mat <- sc_get_expression_matrix(seu, features = c(fx, fy), assay = "RNA", slot = "data")
  x <- as.numeric(mat[fx, ]); y <- as.numeric(mat[fy, ])

  ok <- is.finite(x) & is.finite(y) & !is.na(ctx) & ctx != "Other"
  x <- x[ok]; y <- y[ok]; ctx <- factor(ctx[ok])

  r_global <- suppressWarnings(stats::cor(x, y))
  r_by <- tapply(seq_along(x), ctx, function(idx) suppressWarnings(stats::cor(x[idx], y[idx])))
  n_by <- table(ctx)
  message(sprintf("PBMC %s vs %s: global r = %.3f", fx, fy, r_global))
  for (L in levels(ctx)) message(sprintf("  %s: r=%.2f (n=%d)", L, r_by[L], n_by[L]))

  df <- data.frame(x = x, y = y, celltype = ctx)
  cmap <- sc_context_color_map(levels(ctx))

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = celltype), size = 0.5, alpha = 0.35) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "black", linewidth = 0.9, formula = y ~ x) +
    geom_smooth(aes(group = celltype, color = celltype), method = "lm", se = FALSE,
                linewidth = 0.7, formula = y ~ x) +
    scale_color_manual(values = cmap, name = "Cell type") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.6), ncol = 1)) +
    labs(title = sprintf("PBMC: %s vs %s (global r = %+.2f)", fx, fy, r_global),
         x = sprintf("%s (log-normalized)", fx),
         y = sprintf("%s (log-normalized)", fy)) +
    theme_panel
  list(plot = p, r_global = r_global, r_by = r_by)
}

# ===========================================================================
# Render
# ===========================================================================
W_in <- 90 / 25.4
H_in <- 68 / 25.4

gt <- make_gtex_panel()
ggsave(file.path(out_dir, "gb_fig3A_gtex.pdf"), gt$plot,
       width = W_in, height = H_in, units = "in", device = cairo_pdf)
message("Wrote gb_fig3A_gtex.pdf")

pb <- make_pbmc_panel()
ggsave(file.path(out_dir, "gb_fig3B_pbmc.pdf"), pb$plot,
       width = W_in, height = H_in, units = "in", device = cairo_pdf)
message("Wrote gb_fig3B_pbmc.pdf")

cat("\n=== SUMMARY ===\n")
cat(sprintf("GTEx RTKN2 vs RNASE3 global r = %+.3f\n", gt$r_global))
cat(sprintf("PBMC LYZ vs FTH1 global r = %+.3f\n", pb$r_global))
