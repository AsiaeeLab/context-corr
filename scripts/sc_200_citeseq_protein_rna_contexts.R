#!/usr/bin/env Rscript

source("00-paths.R")
source(file.path("R", "contextcorr.R"))
source(file.path("R", "sc_datasets.R"))
source(file.path("R", "sc_context_pipeline.R"))

sc_require_pkg("ggplot2")

set.seed(200)

experiment_id <- "sc_200_citeseq_protein_rna_contexts"
out_dir <- file.path(paths$results, "sc", experiment_id)
plot_dir <- file.path(out_dir, "plots")
cache_dir <- file.path(paths$results, "cache", experiment_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading PBMC CITE-seq dataset...\n")
ds <- load_pbmc_citeseq(paths = paths, auto_download = FALSE, verbose = TRUE)
seu <- ds$seurat

cat("Preprocessing RNA and ADT...\n")
seu <- sc_preprocess(
  seu,
  assay = "RNA",
  min_features = 200,
  min_counts = NULL,
  max_mt = 20,
  n_hvg = 1800,
  normalize_mode = "log",
  residualize = FALSE,
  vars_to_regress = NULL,
  seed = 200,
  verbose = FALSE
)
Seurat::DefaultAssay(seu) <- "ADT"
seu <- Seurat::NormalizeData(seu, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = FALSE)

cat("Constructing context definitions (RNA / ADT / joint)...\n")
seu <- sc_cluster_context(seu, mode = "rna", context_col = "context_rna", dims_rna = 1:20, resolution = 0.6, seed = 200, verbose = FALSE)
seu <- sc_cluster_context(seu, mode = "adt", context_col = "context_adt", dims_adt = 1:18, resolution = 0.6, seed = 201, verbose = FALSE)
seu <- sc_cluster_context(seu, mode = "joint", context_col = "context_joint", dims_rna = 1:20, dims_adt = 1:18, resolution = 0.6, seed = 202, verbose = FALSE)

context_defs <- list(
  RNA_defined = as.character(seu$context_rna),
  ADT_defined = as.character(seu$context_adt),
  Joint_RNA_ADT = as.character(seu$context_joint)
)

rna_feat <- sc_select_rna_features(
  seu,
  assay = "RNA",
  n_hvg = 1500,
  marker_genes = sc_default_marker_panel(),
  min_detect_frac = 0.01,
  max_features = 1300
)
seu <- rna_feat$seurat

adt_features <- sc_select_adt_features(seu, assay = "ADT", include_isotype = FALSE)
if (length(adt_features) > 80) adt_features <- adt_features[seq_len(80)]

expr_rna <- sc_get_expression_matrix(seu, features = rna_feat$features, assay = "RNA", slot = "data")
expr_adt <- sc_get_expression_matrix(seu, features = adt_features, assay = "ADT", slot = "data")

pair_file <- file.path(cache_dir, "pairs.rds")
if (file.exists(pair_file)) {
  pairs <- readRDS(pair_file)
} else {
  pairs <- sc_make_feature_pairs(rownames(expr_adt), rownames(expr_rna), n_pairs = 200000, seed = 200)
  saveRDS(pairs, pair_file)
}

all_summary <- list()
all_curves <- list()
all_stats <- list()
all_exemplars <- list()

for (nm in names(context_defs)) {
  cat("Running context stats for:", nm, "\n")
  ctx <- context_defs[[nm]]

  stats <- sc_compute_pair_statistics(
    expr_x = expr_adt,
    expr_y = expr_rna,
    pairs = pairs,
    context = ctx,
    eps = 0.05,
    majority = "count",
    n_min = 200,
    method = "pearson",
    cache_file = file.path(cache_dir, paste0("pair_stats_", nm, ".rds")),
    force = FALSE,
    progress_every = 500
  )

  summary_df <- stats$summary_pairs
  summary_df$context_definition <- nm
  curve_df <- sc_reversal_rate_curve(
    r_global = summary_df$r_global,
    r_ctx_mat = stats$r_by_context,
    thresholds = seq(0, 0.5, by = 0.05),
    eps_grid = c(0, 0.05, 0.1)
  )
  curve_df$context_definition <- nm

  ex <- sc_rank_exemplars(summary_df, stats$r_by_context, stats$n_by_context, n_top = 30)

  all_summary[[nm]] <- summary_df
  all_curves[[nm]] <- curve_df
  all_stats[[nm]] <- stats
  all_exemplars[[nm]] <- ex
}

summary_pairs <- do.call(rbind, all_summary)
curve_df <- do.call(rbind, all_curves)

inputs_obj <- list(
  experiment_id = experiment_id,
  dataset = ds[c("dataset_id", "source", "metadata_fields", "notes", "n_cells", "n_features")],
  context_definitions = lapply(context_defs, function(v) sort(table(v), decreasing = TRUE)),
  feature_sets = list(adt_features = rownames(expr_adt), rna_features = rownames(expr_rna)),
  n_pairs = nrow(pairs),
  params = list(
    n_min = 200,
    method = "pearson",
    majority = "count",
    eps = 0.05,
    eps_grid = c(0, 0.05, 0.1)
  )
)

saveRDS(inputs_obj, file.path(out_dir, "inputs.rds"))
utils::write.csv(summary_pairs, file.path(out_dir, "summary_pairs.csv"), row.names = FALSE)
saveRDS(all_exemplars, file.path(out_dir, "top_exemplars.rds"))
utils::write.csv(curve_df, file.path(out_dir, "reversal_rate_curve.csv"), row.names = FALSE)

analysis_ready <- list(
  inputs = inputs_obj,
  pairs = pairs,
  summary_pairs = summary_pairs,
  curve = curve_df,
  context_defs = context_defs,
  per_definition_stats = all_stats,
  top_exemplars = all_exemplars
)
saveRDS(analysis_ready, file.path(out_dir, "analysis_ready.rds"))

p_curve <- ggplot2::ggplot(curve_df, ggplot2::aes(x = abs_r_global_threshold, y = reversal_rate, color = factor(eps))) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 1.2) +
  ggplot2::facet_wrap(~context_definition) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::scale_color_brewer(palette = "Dark2", name = "eps") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(title = "Reversal rates by context definition", x = "|r_global| threshold", y = "Reversal rate")

sc_save_plot_dual(p_curve, file.path(plot_dir, "reversal_rate_vs_threshold"), width = 11, height = 4)

i2_plot_df <- summary_pairs[is.finite(summary_pairs$I2), c("I2", "context_definition")]
p_i2 <- ggplot2::ggplot(i2_plot_df, ggplot2::aes(x = I2)) +
  ggplot2::stat_ecdf(linewidth = 1) +
  ggplot2::facet_wrap(~context_definition) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(title = "I2 ECDF by context definition", x = "I2", y = "ECDF")

sc_save_plot_dual(p_i2, file.path(plot_dir, "I2_distribution_ecdf"), width = 11, height = 4)

scatter_df <- summary_pairs[, c("r_global", "r_resid", "context_definition")]
scatter_df <- scatter_df[is.finite(scatter_df$r_global) & is.finite(scatter_df$r_resid), , drop = FALSE]
scatter_df$sign_flip <- with(scatter_df, sign(r_global) != 0 & sign(r_resid) != 0 & sign(r_global) != sign(r_resid))

p_sc <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = r_global, y = r_resid, color = sign_flip)) +
  ggplot2::geom_point(alpha = 0.25, size = 0.8) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggplot2::facet_wrap(~context_definition) +
  ggplot2::scale_color_manual(values = c("FALSE" = "#1f78b4", "TRUE" = "#d95f02"), name = "Sign flip") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(title = "Global vs residualized correlations", x = "r_global", y = "r_resid")

sc_save_plot_dual(p_sc, file.path(plot_dir, "r_global_vs_r_resid"), width = 11, height = 4)

# Exemplar panel: one top pair per context definition.
ex_counter <- 0
for (nm in names(all_exemplars)) {
  ex <- all_exemplars[[nm]]
  if (length(ex$indices) == 0) next
  idx <- ex$indices[1]
  pair_row <- all_summary[[nm]][idx, , drop = FALSE]
  p_ex <- sc_plot_exemplar_pair(
    expr_x = expr_adt,
    expr_y = expr_rna,
    context = context_defs[[nm]],
    pair_row = pair_row,
    r_ctx = all_stats[[nm]]$r_by_context[idx, , drop = FALSE],
    n_ctx = all_stats[[nm]]$n_by_context[idx, , drop = FALSE],
    seed = 300 + ex_counter
  )
  sc_save_plot_dual(
    p_ex,
    file.path(plot_dir, sprintf("exemplar_%s_%s__%s", nm, pair_row$feature_x, pair_row$feature_y)),
    width = 10,
    height = 4.5
  )
  ex_counter <- ex_counter + 1
  if (ex_counter >= 3) break
}

writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))

cat("Done. Outputs in:\n", out_dir, "\n", sep = "")
