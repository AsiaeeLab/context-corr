#!/usr/bin/env Rscript

source("00-paths.R")
source(file.path("R", "contextcorr.R"))
source(file.path("R", "sc_datasets.R"))
source(file.path("R", "sc_context_pipeline.R"))

set.seed(100)

experiment_id <- "sc_100_pbmc_coarse_celltypes"
cache_root <- file.path(paths$results, "cache")

cat("Loading PBMC scRNA dataset...\n")
ds <- load_pbmc_scrna(paths = paths, auto_download = FALSE, verbose = TRUE)
seu <- ds$seurat

cat("Preprocessing...\n")
seu <- sc_preprocess(
  seu,
  assay = "RNA",
  min_features = 200,
  min_counts = NULL,
  max_mt = 20,
  n_hvg = 2000,
  normalize_mode = "log",
  residualize = FALSE,
  vars_to_regress = NULL,
  seed = 100,
  verbose = FALSE
)

cat("Defining coarse cell-type contexts...\n")
seu <- sc_define_coarse_context(
  seu,
  context_col = "context_coarse",
  prefer_metadata = TRUE,
  marker_sets = sc_marker_gene_sets()$coarse,
  assay = "RNA",
  slot = "data",
  unknown_label = "Other",
  min_score_delta = 0.01
)
ctx <- as.character(seu$context_coarse)

feat_info <- sc_select_rna_features(
  seu,
  assay = "RNA",
  n_hvg = 2000,
  marker_genes = sc_default_marker_panel(),
  min_detect_frac = 0.01,
  max_features = 1800
)
seu <- feat_info$seurat

cat("Running context correlation pipeline...\n")
res <- sc_run_context_experiment(
  experiment_id = experiment_id,
  seurat_obj = seu,
  context = ctx,
  feature_x = feat_info$features,
  feature_y = feat_info$features,
  assay_x = "RNA",
  assay_y = "RNA",
  slot_x = "data",
  slot_y = "data",
  n_pairs = 200000,
  seed = 100,
  n_min = 200,
  method = "pearson",
  majority = "count",
  eps = 0.05,
  eps_grid = c(0, 0.05, 0.1),
  thresholds = seq(0, 0.5, by = 0.05),
  results_root = paths$results,
  cache_root = cache_root,
  top_n = 40,
  force = FALSE,
  robustness = list(
    equalize_context_sizes = FALSE,
    equalize_target_n = NULL,
    run_spearman = FALSE,
    spearman_n_pairs = 20000
  ),
  inputs_extra = list(
    dataset = ds[c("dataset_id", "source", "metadata_fields", "notes", "n_cells", "n_features")],
    feature_selection = list(
      n_hvg = 2000,
      max_features = 1800,
      marker_panel_size = length(sc_default_marker_panel())
    ),
    context_definition = "broad immune cell types (metadata-mapped or marker-scored)"
  ),
  exemplar_plot_n = 3
)

cat("Done. Outputs in:\n", res$out_dir, "\n", sep = "")
