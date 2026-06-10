# Reusable single-cell context-correlation pipeline helpers.
#
# Core correlation definitions are delegated to R/contextcorr.R:
# - context_correlations()
# - heterogeneity_Q_I2()
# - mean_residualized_cor()
# - simpson_flag()

sc_require_contextcorr <- function() {
  needed <- c("context_correlations", "heterogeneity_Q_I2", "mean_residualized_cor", "simpson_flag")
  has_all <- all(vapply(needed, exists, logical(1), inherits = TRUE))
  if (!has_all) {
    if (file.exists(file.path("R", "contextcorr.R"))) {
      source(file.path("R", "contextcorr.R"))
    }
  }
  has_all <- all(vapply(needed, exists, logical(1), inherits = TRUE))
  if (!has_all) {
    stop("Missing core context-correlation helpers. Ensure R/contextcorr.R is sourced.", call. = FALSE)
  }
}

sc_require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package `", pkg, "` is required.", call. = FALSE)
  }
}

sc_marker_gene_sets <- function() {
  list(
    coarse = list(
      T = c("CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "LTB", "MALAT1"),
      B = c("MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "CD37"),
      NK = c("NKG7", "GNLY", "PRF1", "CTSW", "KLRD1", "TRAC"),
      Mono = c("LYZ", "S100A8", "S100A9", "FCN1", "CTSS", "LST1", "FCGR3A"),
      DC = c("FCER1A", "CST3", "CLEC10A", "ITM2B", "HLA-DRA", "LILRA4"),
      Platelet = c("PPBP", "PF4", "SDPR", "GNG11", "NRGN")
    ),
    t_subtypes = list(
      Naive = c("CCR7", "IL7R", "SELL", "TCF7", "LEF1", "MALAT1"),
      Memory = c("LTB", "IL7R", "AQP3", "MALAT1", "MTRNR2L12", "MTRNR2L8"),
      Cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "CTSW", "FGFBP2"),
      Treg = c("FOXP3", "IL2RA", "TIGIT", "CTLA4", "IKZF2", "LAG3"),
      Activated = c("IFNG", "TNFRSF4", "HLA-DRA", "CD69", "FOS", "JUN"),
      Proliferating = c("MKI67", "TOP2A", "TYMS", "STMN1", "UBE2C")
    ),
    adt_coarse = list(
      T = c("CD3", "CD4", "CD8"),
      B = c("CD19", "CD20", "MS4A1", "CD79"),
      NK = c("CD56", "NKG2", "CD16"),
      Mono = c("CD14", "CD16", "CD11b", "CD68"),
      DC = c("CD1C", "CD11c", "CLEC", "CD303")
    )
  )
}

sc_default_marker_panel <- function() {
  sets <- sc_marker_gene_sets()
  unique(unlist(c(sets$coarse, sets$t_subtypes), use.names = FALSE))
}

sc_preprocess <- function(
  seu,
  assay = "RNA",
  min_features = 200,
  min_counts = NULL,
  max_mt = NULL,
  n_hvg = 2000,
  normalize_mode = c("log", "sctransform"),
  residualize = FALSE,
  vars_to_regress = NULL,
  seed = 1,
  verbose = TRUE
) {
  sc_require_pkg("Seurat")
  normalize_mode <- match.arg(normalize_mode)

  if (!inherits(seu, "Seurat")) stop("`seu` must be a Seurat object.", call. = FALSE)
  if (!(assay %in% names(seu@assays))) stop("Assay `", assay, "` not found.", call. = FALSE)

  set.seed(seed)

  Seurat::DefaultAssay(seu) <- assay
  md <- seu@meta.data

  if (!("percent.mt" %in% colnames(md))) {
    feats <- rownames(seu[[assay]])
    mt_pat <- ifelse(any(grepl("^MT-", feats)), "^MT-", "^mt-")
    mt_hits <- sum(grepl(mt_pat, feats))
    if (mt_hits > 0) {
      seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = mt_pat, assay = assay)
    }
  }

  keep <- rep(TRUE, ncol(seu))
  md <- seu@meta.data

  if (!is.null(min_features) && "nFeature_RNA" %in% colnames(md)) {
    keep <- keep & md$nFeature_RNA >= min_features
  }
  if (!is.null(min_counts) && "nCount_RNA" %in% colnames(md)) {
    keep <- keep & md$nCount_RNA >= min_counts
  }
  if (!is.null(max_mt) && "percent.mt" %in% colnames(md)) {
    keep <- keep & md$percent.mt <= max_mt
  }

  if (!all(keep)) {
    seu <- subset(seu, cells = colnames(seu)[which(keep)])
  }

  if (normalize_mode == "log") {
    seu <- Seurat::NormalizeData(seu, assay = assay, normalization.method = "LogNormalize", verbose = verbose)
    seu <- Seurat::FindVariableFeatures(seu, assay = assay, selection.method = "vst", nfeatures = n_hvg, verbose = verbose)
    if (isTRUE(residualize) || (!is.null(vars_to_regress) && length(vars_to_regress) > 0)) {
      feats <- unique(c(Seurat::VariableFeatures(seu), sc_default_marker_panel()))
      feats <- intersect(feats, rownames(seu[[assay]]))
      seu <- Seurat::ScaleData(
        seu,
        assay = assay,
        features = feats,
        vars.to.regress = vars_to_regress,
        verbose = verbose
      )
    }
  } else {
    seu <- Seurat::SCTransform(
      seu,
      assay = assay,
      variable.features.n = n_hvg,
      vars.to.regress = vars_to_regress,
      verbose = verbose
    )
  }

  seu
}

sc_map_coarse_label <- function(x) {
  lx <- tolower(as.character(x))
  if (grepl("cd4|cd8|\\bt cell\\b|\\bt-?cell\\b|\\btreg\\b|\\bt memory\\b|\\bt naive\\b", lx)) return("T")
  if (grepl("\\bb cell\\b|\\bb-?cell\\b|plasma|ms4a1|cd79", lx)) return("B")
  if (grepl("\\bnk\\b|natural killer|nkg7|gnly", lx)) return("NK")
  if (grepl("mono|monocyte|macroph", lx)) return("Mono")
  if (grepl("dendritic|\\bdc\\b|c\\s*dc|p\\s*dc", lx)) return("DC")
  if (grepl("platelet|megakary", lx)) return("Platelet")
  "Other"
}

sc_find_existing_celltype_col <- function(seu) {
  cols <- colnames(seu@meta.data)
  priority <- c(
    "celltype", "cell_type", "CellType", "celltype.l1", "celltype.l2",
    "seurat_annotations", "annotation", "predicted.celltype.l1", "predicted.celltype.l2"
  )
  hit <- intersect(priority, cols)
  if (length(hit) > 0) return(hit[1])

  guess <- cols[grepl("cell|annot|type|lineage", cols, ignore.case = TRUE)]
  if (length(guess) > 0) return(guess[1])
  NULL
}

sc_score_markers <- function(seu, marker_sets, assay = "RNA", slot = "data") {
  mat <- tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = slot),
    error = function(e) Seurat::GetAssayData(seu, assay = assay, slot = slot)
  )
  labs <- names(marker_sets)
  scores <- matrix(NA_real_, nrow = ncol(seu), ncol = length(labs), dimnames = list(colnames(seu), labs))

  for (lab in labs) {
    feats <- intersect(marker_sets[[lab]], rownames(mat))
    if (length(feats) == 0) {
      scores[, lab] <- NA_real_
    } else {
      scores[, lab] <- Matrix::colMeans(mat[feats, , drop = FALSE])
    }
  }

  scores
}

sc_assign_from_scores <- function(scores, unknown_label = "Other", min_delta = 0) {
  if (!is.matrix(scores) || nrow(scores) == 0 || ncol(scores) == 0) {
    return(rep(unknown_label, nrow(scores)))
  }
  max_idx <- max.col(scores, ties.method = "first")
  best <- scores[cbind(seq_len(nrow(scores)), max_idx)]

  tmp <- scores
  tmp[cbind(seq_len(nrow(tmp)), max_idx)] <- -Inf
  second <- apply(tmp, 1, max, na.rm = TRUE)
  second[!is.finite(second)] <- -Inf

  labs <- colnames(scores)[max_idx]
  labs[!is.finite(best)] <- unknown_label
  labs[(best - second) < min_delta] <- unknown_label
  labs
}

sc_define_coarse_context <- function(
  seu,
  context_col = "context_coarse",
  prefer_metadata = TRUE,
  marker_sets = NULL,
  assay = "RNA",
  slot = "data",
  unknown_label = "Other",
  min_score_delta = 0.01
) {
  sc_require_pkg("Seurat")
  if (is.null(marker_sets)) marker_sets <- sc_marker_gene_sets()$coarse

  assigned <- NULL
  source_tag <- "marker"

  if (isTRUE(prefer_metadata)) {
    col <- sc_find_existing_celltype_col(seu)
    if (!is.null(col)) {
      raw <- as.character(seu@meta.data[[col]])
      mapped <- vapply(raw, sc_map_coarse_label, character(1))
      if (length(unique(mapped)) >= 3) {
        assigned <- mapped
        source_tag <- paste0("metadata:", col)
      }
    }
  }

  if (is.null(assigned)) {
    scores <- sc_score_markers(seu, marker_sets = marker_sets, assay = assay, slot = slot)
    assigned <- sc_assign_from_scores(scores, unknown_label = unknown_label, min_delta = min_score_delta)
  }

  seu[[context_col]] <- assigned
  attr(seu[[context_col]], "source") <- source_tag
  seu
}

sc_map_t_subtype <- function(x) {
  lx <- tolower(as.character(x))
  if (grepl("naive|tn", lx)) return("Naive")
  if (grepl("memory|tm", lx)) return("Memory")
  if (grepl("cytotoxic|effector|teff|ctl|cd8", lx)) return("Cytotoxic")
  if (grepl("treg|regulatory", lx)) return("Treg")
  if (grepl("activ|ifn|stim", lx)) return("Activated")
  if (grepl("prolif|cycling|mki67", lx)) return("Proliferating")
  "Other"
}

sc_define_t_subtypes <- function(
  seu,
  context_col = "context_t_subtype",
  prefer_metadata = TRUE,
  marker_sets = NULL,
  assay = "RNA",
  slot = "data",
  min_score_delta = 0.01
) {
  if (is.null(marker_sets)) marker_sets <- sc_marker_gene_sets()$t_subtypes

  assigned <- NULL
  source_tag <- "marker"

  if (isTRUE(prefer_metadata)) {
    col <- sc_find_existing_celltype_col(seu)
    if (!is.null(col)) {
      raw <- as.character(seu@meta.data[[col]])
      mapped <- vapply(raw, sc_map_t_subtype, character(1))
      if (length(unique(mapped)) >= 3) {
        assigned <- mapped
        source_tag <- paste0("metadata:", col)
      }
    }
  }

  if (is.null(assigned)) {
    scores <- sc_score_markers(seu, marker_sets = marker_sets, assay = assay, slot = slot)
    assigned <- sc_assign_from_scores(scores, unknown_label = "Other", min_delta = min_score_delta)
  }

  seu[[context_col]] <- assigned
  attr(seu[[context_col]], "source") <- source_tag
  seu
}

sc_subset_to_context <- function(seu, context_col, keep_values) {
  if (!(context_col %in% colnames(seu@meta.data))) {
    stop("Context column `", context_col, "` not found in metadata.", call. = FALSE)
  }
  keep <- as.character(seu@meta.data[[context_col]]) %in% keep_values
  subset(seu, cells = colnames(seu)[which(keep)])
}

sc_cluster_context <- function(
  seu,
  mode = c("rna", "adt", "joint"),
  context_col,
  dims_rna = 1:20,
  dims_adt = 1:18,
  resolution = 0.6,
  seed = 1,
  verbose = TRUE
) {
  sc_require_pkg("Seurat")
  mode <- match.arg(mode)
  set.seed(seed)

  if (mode %in% c("rna", "joint")) {
    if (!("RNA" %in% names(seu@assays))) stop("RNA assay required for RNA/joint contexts.", call. = FALSE)
    Seurat::DefaultAssay(seu) <- "RNA"
    if (length(Seurat::VariableFeatures(seu)) == 0) {
      seu <- Seurat::FindVariableFeatures(seu, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = verbose)
    }
    if (!("pca" %in% names(seu@reductions))) {
      seu <- Seurat::ScaleData(seu, assay = "RNA", features = Seurat::VariableFeatures(seu), verbose = verbose)
      rna_feats <- Seurat::VariableFeatures(seu)
      npcs_rna <- min(max(dims_rna), max(1, length(rna_feats) - 1))
      seu <- Seurat::RunPCA(seu, assay = "RNA", features = rna_feats, npcs = npcs_rna, verbose = verbose)
    }
  }

  if (mode %in% c("adt", "joint")) {
    if (!("ADT" %in% names(seu@assays))) stop("ADT assay required for ADT/joint contexts.", call. = FALSE)
    Seurat::DefaultAssay(seu) <- "ADT"
    seu <- Seurat::NormalizeData(seu, assay = "ADT", normalization.method = "CLR", margin = 2, verbose = verbose)
    seu <- Seurat::ScaleData(seu, assay = "ADT", verbose = verbose)
    adt_feats <- rownames(seu[["ADT"]])
    npcs_adt <- min(max(dims_adt), max(1, length(adt_feats) - 1))
    seu <- Seurat::RunPCA(
      seu,
      assay = "ADT",
      features = adt_feats,
      reduction.name = "apca",
      npcs = npcs_adt,
      verbose = verbose
    )
  }

  if (mode == "rna") {
    n_pcs <- ncol(Seurat::Embeddings(seu, "pca"))
    dims_use <- dims_rna[dims_rna <= n_pcs]
    if (length(dims_use) == 0) stop("No valid RNA PCA dimensions available for FindNeighbors().", call. = FALSE)
    seu <- Seurat::FindNeighbors(seu, reduction = "pca", dims = dims_use, verbose = verbose)
    seu <- Seurat::FindClusters(seu, resolution = resolution, random.seed = seed, verbose = verbose)
    seu[[context_col]] <- as.character(Seurat::Idents(seu))
  } else if (mode == "adt") {
    n_pcs <- ncol(Seurat::Embeddings(seu, "apca"))
    dims_use <- dims_adt[dims_adt <= n_pcs]
    if (length(dims_use) == 0) stop("No valid ADT PCA dimensions available for FindNeighbors().", call. = FALSE)
    seu <- Seurat::FindNeighbors(seu, reduction = "apca", dims = dims_use, verbose = verbose)
    seu <- Seurat::FindClusters(seu, resolution = resolution, random.seed = seed, verbose = verbose)
    seu[[context_col]] <- as.character(Seurat::Idents(seu))
  } else {
    n_pcs_rna <- ncol(Seurat::Embeddings(seu, "pca"))
    n_pcs_adt <- ncol(Seurat::Embeddings(seu, "apca"))
    dims_rna_use <- dims_rna[dims_rna <= n_pcs_rna]
    dims_adt_use <- dims_adt[dims_adt <= n_pcs_adt]
    if (length(dims_rna_use) == 0) stop("No valid RNA PCA dimensions available for FindMultiModalNeighbors().", call. = FALSE)
    if (length(dims_adt_use) == 0) stop("No valid ADT PCA dimensions available for FindMultiModalNeighbors().", call. = FALSE)
    seu <- Seurat::FindMultiModalNeighbors(
      seu,
      reduction.list = list("pca", "apca"),
      dims.list = list(dims_rna_use, dims_adt_use),
      verbose = verbose
    )
    seu <- Seurat::FindClusters(seu, graph.name = "wsnn", resolution = resolution, random.seed = seed, verbose = verbose)
    seu[[context_col]] <- as.character(Seurat::Idents(seu))
  }

  seu
}

sc_select_rna_features <- function(
  seu,
  assay = "RNA",
  n_hvg = 2000,
  marker_genes = NULL,
  min_detect_frac = 0.01,
  max_features = NULL
) {
  sc_require_pkg("Seurat")

  if (is.null(marker_genes)) marker_genes <- sc_default_marker_panel()
  Seurat::DefaultAssay(seu) <- assay

  hvgs <- Seurat::VariableFeatures(seu)
  if (length(hvgs) == 0) {
    seu <- Seurat::FindVariableFeatures(seu, assay = assay, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
    hvgs <- Seurat::VariableFeatures(seu)
  }
  hvgs <- hvgs[seq_len(min(length(hvgs), n_hvg))]

  counts <- tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = "counts"),
    error = function(e) Seurat::GetAssayData(seu, assay = assay, slot = "counts")
  )
  detect_frac <- Matrix::rowMeans(counts > 0)

  hvgs <- hvgs[hvgs %in% names(detect_frac)[detect_frac >= min_detect_frac]]
  markers_used <- intersect(unique(marker_genes), rownames(counts))

  feat <- unique(c(markers_used, hvgs))
  if (!is.null(max_features) && length(feat) > max_features) {
    keep_hvg <- setdiff(feat, markers_used)
    feat <- unique(c(markers_used, keep_hvg[seq_len(max(0, max_features - length(markers_used)))]))
  }

  list(
    features = feat,
    hvg = hvgs,
    marker_genes_used = markers_used,
    detect_frac = detect_frac[feat],
    seurat = seu
  )
}

sc_select_adt_features <- function(seu, assay = "ADT", include_isotype = FALSE) {
  if (!(assay %in% names(seu@assays))) stop("Assay `", assay, "` not found.", call. = FALSE)
  feats <- rownames(seu[[assay]])

  if (!include_isotype) {
    drop_pat <- "isotype|ctrl|control|igg"
    feats <- feats[!grepl(drop_pat, feats, ignore.case = TRUE)]
  }

  feats
}

sc_make_feature_pairs <- function(features_x, features_y = NULL, n_pairs = 200000, seed = 1) {
  set.seed(seed)
  if (is.null(features_y)) features_y <- features_x

  nx <- length(features_x)
  ny <- length(features_y)
  if (nx < 1 || ny < 1) stop("Feature lists cannot be empty.", call. = FALSE)

  symmetric <- identical(features_x, features_y)

  if (symmetric) {
    if (nx < 2) stop("Need at least 2 features for symmetric pair sampling.", call. = FALSE)
    total <- as.integer(nx * (nx - 1) / 2)

    if (n_pairs >= total && total <= 3000000) {
      cmb <- utils::combn(nx, 2)
      idx_x <- cmb[1, ]
      idx_y <- cmb[2, ]
    } else {
      target <- min(n_pairs, total)
      out <- data.frame(idx_x = integer(0), idx_y = integer(0))
      while (nrow(out) < target) {
        chunk <- max(10000, target * 2)
        ix <- sample.int(nx, size = chunk, replace = TRUE)
        iy <- sample.int(nx, size = chunk, replace = TRUE)
        ok <- ix < iy
        if (!any(ok)) next
        cand <- unique(data.frame(idx_x = ix[ok], idx_y = iy[ok]))
        out <- unique(rbind(out, cand))
      }
      out <- out[seq_len(target), , drop = FALSE]
      idx_x <- out$idx_x
      idx_y <- out$idx_y
    }
  } else {
    total <- as.integer(nx * ny)
    if (n_pairs >= total && total <= 3000000) {
      grid <- expand.grid(idx_x = seq_len(nx), idx_y = seq_len(ny), KEEP.OUT.ATTRS = FALSE)
      idx_x <- grid$idx_x
      idx_y <- grid$idx_y
    } else {
      target <- min(n_pairs, total)
      lin <- sample.int(total, size = target, replace = FALSE)
      idx_x <- ((lin - 1L) %% nx) + 1L
      idx_y <- ((lin - 1L) %/% nx) + 1L
    }
  }

  data.frame(
    pair_id = seq_along(idx_x),
    idx_x = as.integer(idx_x),
    idx_y = as.integer(idx_y),
    feature_x = features_x[idx_x],
    feature_y = features_y[idx_y],
    stringsAsFactors = FALSE
  )
}

sc_get_expression_matrix <- function(seu, features, assay = "RNA", slot = "data") {
  sc_require_pkg("Seurat")
  mat <- tryCatch(
    Seurat::GetAssayData(seu, assay = assay, layer = slot),
    error = function(e) Seurat::GetAssayData(seu, assay = assay, slot = slot)
  )
  feats <- intersect(features, rownames(mat))
  if (length(feats) == 0) {
    stop("No requested features found in assay=", assay, ", slot=", slot, ".", call. = FALSE)
  }
  as.matrix(mat[feats, , drop = FALSE])
}

sc_equalize_context_sizes <- function(context, seed = 1, target_n = NULL) {
  set.seed(seed)
  fctx <- as.factor(context)
  tab <- table(fctx)
  if (length(tab) < 2) return(seq_along(context))
  if (is.null(target_n)) target_n <- min(tab)

  idx <- integer(0)
  for (lv in names(tab)) {
    ii <- which(fctx == lv)
    if (length(ii) <= target_n) {
      idx <- c(idx, ii)
    } else {
      idx <- c(idx, sample(ii, size = target_n))
    }
  }
  sort(idx)
}

sc_row_zscore_matrix <- function(mat) {
  stopifnot(is.matrix(mat))
  n <- ncol(mat)
  if (n < 2) stop("Need at least 2 columns to compute correlations.", call. = FALSE)
  denom <- n - 1

  mu <- rowMeans(mat)
  centered <- mat - mu
  s2 <- rowSums(centered ^ 2) / denom
  s <- sqrt(s2)
  z <- centered / s
  z[!is.finite(z)] <- NA_real_
  z
}

sc_pairwise_cor_from_z <- function(zx, zy, idx_x, idx_y) {
  stopifnot(is.matrix(zx), is.matrix(zy), ncol(zx) == ncol(zy))
  denom <- ncol(zx) - 1
  if (denom <= 0) return(rep(NA_real_, length(idx_x)))
  rowSums(zx[idx_x, , drop = FALSE] * zy[idx_y, , drop = FALSE]) / denom
}

sc_residualize_matrix_by_context <- function(mat, context_factor) {
  stopifnot(is.matrix(mat), is.factor(context_factor), length(context_factor) == ncol(mat))
  out <- mat
  lev <- levels(context_factor)
  for (lv in lev) {
    idx <- which(context_factor == lv)
    if (length(idx) < 1) next
    mu <- rowMeans(out[, idx, drop = FALSE])
    out[, idx] <- out[, idx, drop = FALSE] - mu
  }
  out
}

sc_compute_pair_statistics_slow <- function(
  expr_x,
  expr_y,
  pairs,
  context,
  eps = 0.05,
  majority = "count",
  n_min = 200,
  method = "pearson",
  cache_file = NULL,
  force = FALSE,
  progress_every = 1000
) {
  sc_require_contextcorr()

  if (!is.null(cache_file) && file.exists(cache_file) && !isTRUE(force)) {
    return(readRDS(cache_file))
  }

  if (ncol(expr_x) != ncol(expr_y)) stop("expr_x and expr_y must have matching cells.", call. = FALSE)
  if (length(context) != ncol(expr_x)) stop("context length must equal number of cells.", call. = FALSE)

  context <- as.factor(context)
  context_levels <- levels(context)

  n_pairs <- nrow(pairs)
  k_ctx <- length(context_levels)

  r_global <- rep(NA_real_, n_pairs)
  r_resid <- rep(NA_real_, n_pairs)
  flag_simpson <- rep(FALSE, n_pairs)
  majority_sign <- rep(0L, n_pairs)
  n_pos <- rep(0L, n_pairs)
  n_neg <- rep(0L, n_pairs)
  n_used <- rep(0L, n_pairs)
  Q <- rep(NA_real_, n_pairs)
  p_heterogeneity <- rep(NA_real_, n_pairs)
  I2 <- rep(NA_real_, n_pairs)
  k_contexts_used <- rep(0L, n_pairs)
  min_n_context <- rep(NA_real_, n_pairs)
  median_n_context <- rep(NA_real_, n_pairs)
  min_abs_r_context <- rep(NA_real_, n_pairs)
  max_abs_r_context <- rep(NA_real_, n_pairs)

  r_ctx_mat <- matrix(NA_real_, nrow = n_pairs, ncol = k_ctx, dimnames = list(NULL, context_levels))
  n_ctx_mat <- matrix(NA_integer_, nrow = n_pairs, ncol = k_ctx, dimnames = list(NULL, context_levels))

  for (i in seq_len(n_pairs)) {
    x <- as.numeric(expr_x[pairs$idx_x[i], ])
    y <- as.numeric(expr_y[pairs$idx_y[i], ])

    cc <- context_correlations(x, y, context = context, method = method)
    r_ctx <- cc$r_by_context
    n_ctx <- cc$n_by_context

    use_ctx <- n_ctx >= n_min
    r_ctx_use <- r_ctx
    r_ctx_use[!use_ctx] <- NA_real_

    sf <- simpson_flag(cc$r_global, r_ctx_use, eps = eps, majority = majority)
    ht <- heterogeneity_Q_I2(r_ctx_use, n_ctx)

    ctx_ok <- as.character(context)
    if (any(!use_ctx)) {
      ctx_ok[!(ctx_ok %in% names(use_ctx)[use_ctx])] <- NA_character_
    }
    r_res <- mean_residualized_cor(x, y, context = ctx_ok, method = method)

    used <- is.finite(r_ctx_use)

    r_global[i] <- cc$r_global
    r_resid[i] <- r_res
    flag_simpson[i] <- isTRUE(sf$flag)
    majority_sign[i] <- as.integer(sf$majority_sign)
    n_pos[i] <- as.integer(sf$n_pos)
    n_neg[i] <- as.integer(sf$n_neg)
    n_used[i] <- as.integer(sf$n_used)

    Q[i] <- ht$Q
    p_heterogeneity[i] <- ht$p
    I2[i] <- ht$I2

    k_contexts_used[i] <- sum(used)
    if (any(used)) {
      min_n_context[i] <- min(n_ctx[used])
      median_n_context[i] <- stats::median(n_ctx[used])
      min_abs_r_context[i] <- min(abs(r_ctx_use[used]))
      max_abs_r_context[i] <- max(abs(r_ctx_use[used]))
    }

    r_ctx_mat[i, ] <- as.numeric(r_ctx_use[context_levels])
    n_ctx_mat[i, ] <- as.integer(n_ctx[context_levels])

    if (i == 1L || (i %% progress_every) == 0L || i == n_pairs) {
      message(sprintf("Processed %d/%d pairs", i, n_pairs))
    }
  }

  summary_pairs <- data.frame(
    feature_x = pairs$feature_x,
    feature_y = pairs$feature_y,
    pair_id = pairs$pair_id,
    idx_x = pairs$idx_x,
    idx_y = pairs$idx_y,
    r_global = r_global,
    r_resid = r_resid,
    flag_simpson = flag_simpson,
    majority_sign = majority_sign,
    n_pos = n_pos,
    n_neg = n_neg,
    n_used = n_used,
    eps = eps,
    Q = Q,
    p_heterogeneity = p_heterogeneity,
    I2 = I2,
    k_contexts_used = k_contexts_used,
    min_n_context = min_n_context,
    median_n_context = median_n_context,
    min_abs_r_context = min_abs_r_context,
    max_abs_r_context = max_abs_r_context,
    stringsAsFactors = FALSE
  )

  out <- list(
    summary_pairs = summary_pairs,
    pairs = pairs,
    r_by_context = r_ctx_mat,
    n_by_context = n_ctx_mat,
    context_levels = context_levels,
    params = list(eps = eps, majority = majority, n_min = n_min, method = method)
  )

  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cache_file)
  }

  out
}

sc_compute_pair_statistics_fast_pearson <- function(
  expr_x,
  expr_y,
  pairs,
  context,
  eps = 0.05,
  majority = "count",
  n_min = 200,
  cache_file = NULL,
  force = FALSE,
  progress_every = 20000,
  chunk_size = 2000,
  validate_n = 25,
  validate_tol = 1e-6
) {
  sc_require_contextcorr()

  if (!is.null(cache_file) && file.exists(cache_file) && !isTRUE(force)) {
    return(readRDS(cache_file))
  }

  if (!is.matrix(expr_x) || !is.matrix(expr_y)) stop("expr_x and expr_y must be matrices (features x cells).", call. = FALSE)
  if (ncol(expr_x) != ncol(expr_y)) stop("expr_x and expr_y must have matching cells.", call. = FALSE)
  if (length(context) != ncol(expr_x)) stop("context length must equal number of cells.", call. = FALSE)

  context <- as.factor(context)
  context_levels <- levels(context)
  n_by_context <- setNames(as.integer(table(context)[context_levels]), context_levels)
  k_ctx <- length(context_levels)
  n_pairs <- nrow(pairs)

  # Preallocate.
  r_global <- rep(NA_real_, n_pairs)
  r_resid <- rep(NA_real_, n_pairs)
  r_ctx_mat <- matrix(NA_real_, nrow = n_pairs, ncol = k_ctx, dimnames = list(NULL, context_levels))
  n_ctx_mat <- matrix(NA_integer_, nrow = n_pairs, ncol = k_ctx, dimnames = list(NULL, context_levels))
  for (j in seq_len(k_ctx)) n_ctx_mat[, j] <- n_by_context[j]

  # Global correlations.
  zx <- sc_row_zscore_matrix(expr_x)
  zy <- if (identical(expr_x, expr_y)) zx else sc_row_zscore_matrix(expr_y)

  idx_x <- pairs$idx_x
  idx_y <- pairs$idx_y

  n_chunks <- ceiling(n_pairs / chunk_size)
  for (ch in seq_len(n_chunks)) {
    a <- (ch - 1L) * chunk_size + 1L
    b <- min(n_pairs, ch * chunk_size)
    ii <- a:b
    r_global[ii] <- sc_pairwise_cor_from_z(zx, zy, idx_x[ii], idx_y[ii])
    if (ch == 1L || (b %% progress_every) == 0L || b == n_pairs) {
      message(sprintf("Global correlations: %d/%d pairs", b, n_pairs))
    }
  }

  # Within-context correlations with n_min enforcement.
  for (j in seq_len(k_ctx)) {
    lv <- context_levels[j]
    cells <- which(context == lv)
    n_ctx <- length(cells)
    if (n_ctx < n_min || n_ctx < 3) {
      next
    }

    zx_c <- sc_row_zscore_matrix(expr_x[, cells, drop = FALSE])
    zy_c <- if (identical(expr_x, expr_y)) zx_c else sc_row_zscore_matrix(expr_y[, cells, drop = FALSE])

    for (ch in seq_len(n_chunks)) {
      a <- (ch - 1L) * chunk_size + 1L
      b <- min(n_pairs, ch * chunk_size)
      ii <- a:b
      r_ctx_mat[ii, j] <- sc_pairwise_cor_from_z(zx_c, zy_c, idx_x[ii], idx_y[ii])
      if (j == 1L && (ch == 1L || (b %% progress_every) == 0L || b == n_pairs)) {
        message(sprintf("Within-context correlations (example context '%s'): %d/%d pairs", lv, b, n_pairs))
      }
    }
  }

  # Residualized correlations: exclude contexts with n < n_min (matches mean_residualized_cor usage in slow path).
  good_levels <- context_levels[n_by_context >= n_min]
  keep_cells <- which(as.character(context) %in% good_levels)
  if (length(keep_cells) >= 4) {
    ctx_keep <- droplevels(context[keep_cells])
    x_keep <- expr_x[, keep_cells, drop = FALSE]
    y_keep <- expr_y[, keep_cells, drop = FALSE]

    x_res <- sc_residualize_matrix_by_context(x_keep, ctx_keep)
    y_res <- if (identical(expr_x, expr_y)) x_res else sc_residualize_matrix_by_context(y_keep, ctx_keep)

    zx_res <- sc_row_zscore_matrix(x_res)
    zy_res <- if (identical(expr_x, expr_y)) zx_res else sc_row_zscore_matrix(y_res)

    for (ch in seq_len(n_chunks)) {
      a <- (ch - 1L) * chunk_size + 1L
      b <- min(n_pairs, ch * chunk_size)
      ii <- a:b
      r_resid[ii] <- sc_pairwise_cor_from_z(zx_res, zy_res, idx_x[ii], idx_y[ii])
      if (ch == 1L || (b %% progress_every) == 0L || b == n_pairs) {
        message(sprintf("Residualized correlations: %d/%d pairs", b, n_pairs))
      }
    }
  }

  # Simpson majority/sign counts (vectorized).
  if (majority == "count") {
    n_pos <- rowSums(is.finite(r_ctx_mat) & (r_ctx_mat >= eps))
    n_neg <- rowSums(is.finite(r_ctx_mat) & (r_ctx_mat <= -eps))
    n_used <- n_pos + n_neg
    majority_sign <- integer(n_pairs)
    has_used <- n_used > 0
    majority_sign[has_used] <- sign(n_pos[has_used] - n_neg[has_used])
  } else {
    # Median of within-context r values after eps thresholding.
    r_use <- r_ctx_mat
    r_use[!(is.finite(r_use) & abs(r_use) >= eps)] <- NA_real_
    n_used <- rowSums(is.finite(r_use))
    majority_sign <- apply(r_use, 1, function(v) {
      v <- v[is.finite(v)]
      if (length(v) == 0) return(0L)
      as.integer(sign(stats::median(v)))
    })
    n_pos <- rowSums(is.finite(r_use) & (r_use > 0))
    n_neg <- rowSums(is.finite(r_use) & (r_use < 0))
  }

  global_sign <- sign(r_global)
  flag_simpson <- is.finite(r_global) & (global_sign != 0) & (majority_sign != 0) & (global_sign != majority_sign)
  k_contexts_used <- rowSums(is.finite(r_ctx_mat))

  # Vectorized heterogeneity Q/I2 (matches heterogeneity_Q_I2 algebra).
  # Only contexts with n>=4 and finite r contribute; r_ctx_mat already NA for n<n_min.
  w <- pmax(n_by_context - 3L, 0L)
  w[n_by_context < 4L] <- 0L
  w_mat <- matrix(rep(as.numeric(w), each = n_pairs), nrow = n_pairs, ncol = k_ctx)

  r_for_z <- r_ctx_mat
  r_for_z[!is.finite(r_for_z)] <- NA_real_
  z <- contextcorr_fisher_z(r_for_z)
  valid <- is.finite(z) & (w_mat > 0)
  z[!valid] <- NA_real_

  sum_w <- rowSums(w_mat * valid)
  sum_wz <- rowSums(w_mat * z, na.rm = TRUE)
  sum_wz2 <- rowSums(w_mat * (z ^ 2), na.rm = TRUE)
  k_valid <- rowSums(valid)
  df <- pmax(k_valid - 1L, 0L)

  Q <- rep(NA_real_, n_pairs)
  p_heterogeneity <- rep(NA_real_, n_pairs)
  I2 <- rep(NA_real_, n_pairs)

  ok <- df >= 1L & is.finite(sum_w) & (sum_w > 0)
  Q[ok] <- sum_wz2[ok] - (sum_wz[ok] ^ 2) / sum_w[ok]
  p_heterogeneity[ok] <- stats::pchisq(Q[ok], df = df[ok], lower.tail = FALSE)
  I2[ok] <- ifelse(is.finite(Q[ok]) & Q[ok] > 0, pmax(0, (Q[ok] - df[ok]) / Q[ok]), NA_real_)

  # Context-size summaries (based on contexts passing n_min).
  use_n <- as.numeric(n_by_context[n_by_context >= n_min])
  min_n_context <- if (length(use_n) > 0) rep(min(use_n), n_pairs) else rep(NA_real_, n_pairs)
  median_n_context <- if (length(use_n) > 0) rep(stats::median(use_n), n_pairs) else rep(NA_real_, n_pairs)

  # min/max abs r across used contexts (after n_min filtering).
  abs_r <- abs(r_ctx_mat)
  abs_r_min <- abs_r
  abs_r_min[!is.finite(abs_r_min)] <- Inf
  min_abs_r_context <- abs_r_min[, 1]
  if (k_ctx >= 2) {
    for (j in 2:k_ctx) min_abs_r_context <- pmin(min_abs_r_context, abs_r_min[, j])
  }
  min_abs_r_context[is.infinite(min_abs_r_context)] <- NA_real_

  abs_r_max <- abs_r
  abs_r_max[!is.finite(abs_r_max)] <- -Inf
  max_abs_r_context <- abs_r_max[, 1]
  if (k_ctx >= 2) {
    for (j in 2:k_ctx) max_abs_r_context <- pmax(max_abs_r_context, abs_r_max[, j])
  }
  max_abs_r_context[!is.finite(max_abs_r_context)] <- NA_real_

  # Validate against the canonical helper functions on a small random subset.
  if (n_pairs > 0 && validate_n > 0) {
    set.seed(1)
    test_idx <- sample.int(n_pairs, size = min(validate_n, n_pairs))
    for (i in test_idx) {
      x <- as.numeric(expr_x[idx_x[i], ])
      y <- as.numeric(expr_y[idx_y[i], ])
      cc <- context_correlations(x, y, context = context, method = "pearson")
      r_ctx_use <- cc$r_by_context
      r_ctx_use[cc$n_by_context < n_min] <- NA_real_

      if (is.finite(cc$r_global) && is.finite(r_global[i])) {
        if (abs(cc$r_global - r_global[i]) > validate_tol) {
          stop("Fast/slow mismatch on r_global at pair ", i, call. = FALSE)
        }
      }

      fast_vec <- as.numeric(r_ctx_mat[i, ])
      slow_vec <- as.numeric(r_ctx_use[context_levels])
      diff <- abs(fast_vec - slow_vec)
      ok_diff <- (is.finite(diff) & diff <= validate_tol) | (!is.finite(fast_vec) & !is.finite(slow_vec))
      if (!all(ok_diff)) {
        stop("Fast/slow mismatch on r_by_context at pair ", i, call. = FALSE)
      }

      sf <- simpson_flag(cc$r_global, r_ctx_use, eps = eps, majority = majority)
      if (!identical(isTRUE(sf$flag), isTRUE(flag_simpson[i]))) stop("Fast/slow mismatch on simpson_flag at pair ", i, call. = FALSE)

      ht <- heterogeneity_Q_I2(r_ctx_use, cc$n_by_context)
      if (is.finite(ht$Q) && is.finite(Q[i]) && abs(ht$Q - Q[i]) > 1e-5) stop("Fast/slow mismatch on Q at pair ", i, call. = FALSE)
      if (is.finite(ht$I2) && is.finite(I2[i]) && abs(ht$I2 - I2[i]) > 1e-5) stop("Fast/slow mismatch on I2 at pair ", i, call. = FALSE)

      ctx_ok <- as.character(context)
      ctx_ok[!(ctx_ok %in% names(cc$n_by_context)[cc$n_by_context >= n_min])] <- NA_character_
      r_res_slow <- mean_residualized_cor(x, y, context = ctx_ok, method = "pearson")
      if (is.finite(r_res_slow) && is.finite(r_resid[i]) && abs(r_res_slow - r_resid[i]) > 1e-5) {
        stop("Fast/slow mismatch on r_resid at pair ", i, call. = FALSE)
      }
    }
  }

  summary_pairs <- data.frame(
    feature_x = pairs$feature_x,
    feature_y = pairs$feature_y,
    pair_id = pairs$pair_id,
    idx_x = idx_x,
    idx_y = idx_y,
    r_global = r_global,
    r_resid = r_resid,
    flag_simpson = flag_simpson,
    majority_sign = as.integer(majority_sign),
    n_pos = as.integer(n_pos),
    n_neg = as.integer(n_neg),
    n_used = as.integer(n_used),
    eps = eps,
    Q = Q,
    p_heterogeneity = p_heterogeneity,
    I2 = I2,
    k_contexts_used = as.integer(k_contexts_used),
    min_n_context = min_n_context,
    median_n_context = median_n_context,
    min_abs_r_context = min_abs_r_context,
    max_abs_r_context = max_abs_r_context,
    stringsAsFactors = FALSE
  )

  out <- list(
    summary_pairs = summary_pairs,
    pairs = pairs,
    r_by_context = r_ctx_mat,
    n_by_context = n_ctx_mat,
    context_levels = context_levels,
    params = list(eps = eps, majority = majority, n_min = n_min, method = "pearson")
  )

  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cache_file)
  }

  out
}

sc_compute_pair_statistics <- function(
  expr_x,
  expr_y,
  pairs,
  context,
  eps = 0.05,
  majority = "count",
  n_min = 200,
  method = "pearson",
  cache_file = NULL,
  force = FALSE,
  progress_every = 1000
) {
  method <- match.arg(method, c("pearson", "spearman"))
  majority <- match.arg(majority, c("count", "median"))
  if (method == "pearson") {
    return(sc_compute_pair_statistics_fast_pearson(
      expr_x = expr_x,
      expr_y = expr_y,
      pairs = pairs,
      context = context,
      eps = eps,
      majority = majority,
      n_min = n_min,
      cache_file = cache_file,
      force = force,
      progress_every = max(1000, progress_every)
    ))
  }
  sc_compute_pair_statistics_slow(
    expr_x = expr_x,
    expr_y = expr_y,
    pairs = pairs,
    context = context,
    eps = eps,
    majority = majority,
    n_min = n_min,
    method = method,
    cache_file = cache_file,
    force = force,
    progress_every = progress_every
  )
}

sc_reversal_rate_curve <- function(
  r_global,
  r_ctx_mat,
  thresholds = seq(0, 0.5, by = 0.05),
  eps_grid = c(0, 0.05, 0.1)
) {
  out <- data.frame()
  global_sign <- sign(r_global)

  for (eps in eps_grid) {
    pos <- rowSums(r_ctx_mat >= eps, na.rm = TRUE)
    neg <- rowSums(r_ctx_mat <= -eps, na.rm = TRUE)
    maj <- sign(pos - neg)

    for (thr in thresholds) {
      elig <- is.finite(r_global) & abs(r_global) >= thr & global_sign != 0 & maj != 0
      frac <- if (sum(elig) > 0) mean(global_sign[elig] != maj[elig]) else NA_real_
      out <- rbind(
        out,
        data.frame(
          eps = eps,
          abs_r_global_threshold = thr,
          eligible_pairs = sum(elig),
          reversal_rate = frac,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  out
}

sc_rank_exemplars <- function(summary_pairs, r_ctx_mat, n_ctx_mat, n_top = 30) {
  score <- with(
    summary_pairs,
    ifelse(flag_simpson, 100, 0) +
      20 * pmin(abs(r_global), 1) +
      10 * ifelse(is.finite(I2), I2, 0) +
      pmin(abs(r_global - r_resid), 1)
  )
  ord <- order(score, decreasing = TRUE, na.last = TRUE)
  keep <- ord[seq_len(min(length(ord), n_top))]

  exemplars <- lapply(keep, function(i) {
    list(
      pair = summary_pairs[i, , drop = FALSE],
      per_context = data.frame(
        context = colnames(r_ctx_mat),
        r_c = as.numeric(r_ctx_mat[i, ]),
        n_c = as.integer(n_ctx_mat[i, ]),
        stringsAsFactors = FALSE
      )
    )
  })

  list(indices = keep, exemplars = exemplars)
}

sc_plot_reversal_curve <- function(curve_df, title = "Reversal rate vs |r_global| threshold") {
  sc_require_pkg("ggplot2")
  ggplot2::ggplot(curve_df, ggplot2::aes(x = abs_r_global_threshold, y = reversal_rate, color = factor(eps))) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_brewer(palette = "Dark2", name = "eps") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = "|r_global| threshold", y = "Simpson reversal rate")
}

sc_plot_I2_distribution <- function(summary_pairs, title = "I2 distribution") {
  sc_require_pkg("ggplot2")
  dat <- summary_pairs[is.finite(summary_pairs$I2), , drop = FALSE]
  ggplot2::ggplot(dat, ggplot2::aes(x = I2)) +
    ggplot2::stat_ecdf(linewidth = 1) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = "I2", y = "ECDF")
}

sc_plot_global_vs_resid <- function(summary_pairs, title = "Global vs residualized correlation") {
  sc_require_pkg("ggplot2")
  dat <- summary_pairs
  dat$sign_flip_global_resid <- with(dat, sign(r_global) != 0 & sign(r_resid) != 0 & sign(r_global) != sign(r_resid))
  frac <- mean(dat$sign_flip_global_resid, na.rm = TRUE)

  ggplot2::ggplot(dat, ggplot2::aes(x = r_global, y = r_resid, color = sign_flip_global_resid)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("FALSE" = "#1f78b4", "TRUE" = "#d95f02"), name = "Sign flip") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = sprintf("%s (flip fraction=%.3f)", title, frac),
      x = "r_global",
      y = "r_resid"
    )
}

sc_plot_exemplar_pair <- function(
  expr_x,
  expr_y,
  context,
  pair_row,
  r_ctx,
  n_ctx,
  max_points = 5000,
  seed = 1
) {
  sc_require_pkg("ggplot2")

  x <- as.numeric(expr_x[pair_row$idx_x, ])
  y <- as.numeric(expr_y[pair_row$idx_y, ])
  ctx <- as.character(context)

  df <- data.frame(x = x, y = y, context = ctx, stringsAsFactors = FALSE)
  if (nrow(df) > max_points) {
    set.seed(seed)
    df <- df[sample(seq_len(nrow(df)), size = max_points), , drop = FALSE]
  }

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = context)) +
    ggplot2::geom_point(alpha = 0.45, size = 0.8) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = sprintf("%s vs %s", pair_row$feature_x, pair_row$feature_y),
      x = pair_row$feature_x,
      y = pair_row$feature_y,
      color = "context"
    )

  dfc <- data.frame(context = colnames(r_ctx), r_c = as.numeric(r_ctx), n_c = as.numeric(n_ctx), stringsAsFactors = FALSE)
  p2 <- ggplot2::ggplot(dfc, ggplot2::aes(x = stats::reorder(context, r_c), y = r_c, fill = n_c)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(title = "Within-context correlations", x = "context", y = "r_c", fill = "n_c")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p1 + p2 + patchwork::plot_layout(widths = c(2, 1))
  } else {
    p1
  }
}

sc_context_color_map <- function(context_levels) {
  stopifnot(is.character(context_levels))
  k <- length(context_levels)
  if (k == 0) return(setNames(character(), character()))

  if (k <= 12) {
    pal <- grDevices::hcl.colors(k, palette = "Dark 3")
  } else {
    pal <- grDevices::hcl.colors(k, palette = "Dynamic")
  }
  setNames(pal, context_levels)
}

sc_pick_highlight_contexts <- function(r_by_context, n_by_context, eps = 0.05, top_k = 6, mode = c("majority", "extreme")) {
  mode <- match.arg(mode)

  r <- as.numeric(r_by_context)
  n <- as.numeric(n_by_context)
  nm <- names(r_by_context)
  if (is.null(nm)) nm <- names(n_by_context)
  if (is.null(nm)) nm <- as.character(seq_along(r))
  names(r) <- nm
  names(n) <- nm

  valid <- is.finite(r) & is.finite(n) & n >= 3
  r <- r[valid]
  n <- n[valid]
  if (length(r) == 0) return(character())

  r_use <- r
  r_use[abs(r_use) < eps] <- NA_real_

  if (mode == "majority") {
    sf <- simpson_flag(r_global = 1, r_by_context = r_use, eps = eps, majority = "count")
    maj <- sf$majority_sign
    if (!is.finite(maj) || maj == 0) {
      mode <- "extreme"
    } else {
      keep <- names(r_use)[is.finite(r_use) & (sign(r_use) == maj)]
      if (length(keep) == 0) mode <- "extreme"
      else {
        ord <- order(abs(r_use[keep]), decreasing = TRUE, na.last = TRUE)
        return(keep[ord][seq_len(min(length(ord), top_k))])
      }
    }
  }

  keep <- names(r_use)[is.finite(r_use)]
  if (length(keep) == 0) keep <- names(r)
  ord <- order(abs(r[keep]), decreasing = TRUE, na.last = TRUE)
  keep[ord][seq_len(min(length(ord), top_k))]
}

sc_plot_exemplar_paperstyle <- function(
  x,
  y,
  context,
  feature_x = "x",
  feature_y = "y",
  eps = 0.05,
  n_min = 200,
  highlight_top_k = 6,
  highlight_mode = c("majority", "extreme"),
  max_points = 20000,
  seed = 1
) {
  sc_require_contextcorr()
  sc_require_pkg("ggplot2")
  sc_require_pkg("patchwork")

  highlight_mode <- match.arg(highlight_mode)

  stopifnot(length(x) == length(y), length(x) == length(context))
  ok <- is.finite(x) & is.finite(y) & !is.na(context)
  x <- x[ok]
  y <- y[ok]
  context <- as.factor(context[ok])

  cc <- context_correlations(x, y, context = context, method = "pearson")
  r_ctx <- cc$r_by_context
  n_ctx <- cc$n_by_context
  r_ctx_use <- r_ctx
  r_ctx_use[n_ctx < n_min] <- NA_real_

  ctx_ok <- as.character(context)
  ctx_ok[!(ctx_ok %in% names(n_ctx)[n_ctx >= n_min])] <- NA_character_
  r_resid <- mean_residualized_cor(x, y, context = ctx_ok, method = "pearson")

  means_x <- tapply(x, context, mean)
  means_y <- tapply(y, context, mean)
  keep_means <- names(n_ctx)[n_ctx >= n_min]
  means_x <- means_x[keep_means]
  means_y <- means_y[keep_means]
  r_means <- suppressWarnings(stats::cor(as.numeric(means_x), as.numeric(means_y), method = "pearson"))

  highlight <- sc_pick_highlight_contexts(r_ctx_use, n_ctx, eps = eps, top_k = highlight_top_k, mode = highlight_mode)
  cmap <- sc_context_color_map(highlight)

  df <- data.frame(x = x, y = y, context = as.character(context), stringsAsFactors = FALSE)
  if (nrow(df) > max_points) {
    set.seed(seed)
    df <- df[sample.int(nrow(df), size = max_points), , drop = FALSE]
  }
  df$context_hl <- ifelse(df$context %in% highlight, df$context, "Other")

  # Panel A: global scatter.
  pA <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(color = "grey70", alpha = 0.35, size = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = sprintf("Global r = %.3f", cc$r_global),
      x = feature_x,
      y = feature_y
    )

  # Panel B: context correlation bars.
  dbar <- data.frame(
    context = names(r_ctx_use),
    r_c = as.numeric(r_ctx_use),
    n_c = as.numeric(n_ctx[names(r_ctx_use)]),
    stringsAsFactors = FALSE
  )
  dbar$highlight <- dbar$context %in% highlight
  dbar <- dbar[order(dbar$r_c, decreasing = TRUE), , drop = FALSE]
  dbar$context <- factor(dbar$context, levels = dbar$context)
  dbar$fill <- ifelse(dbar$highlight, dbar$context, "Other")
  fill_map <- c(cmap, Other = "grey75")

  pB <- ggplot2::ggplot(dbar, ggplot2::aes(x = context, y = r_c, fill = fill)) +
    ggplot2::geom_col(color = "white", linewidth = 0.2, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = cc$r_global, color = "black", linewidth = 0.7) +
    ggplot2::scale_fill_manual(values = fill_map, guide = "none") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggplot2::labs(title = "Context correlations", x = NULL, y = "r_c")

  # Panel C: context means scatter.
  dmeans <- data.frame(
    context = names(means_x),
    x_mean = as.numeric(means_x),
    y_mean = as.numeric(means_y),
    n_c = as.numeric(n_ctx[names(means_x)]),
    stringsAsFactors = FALSE
  )
  dmeans$context_hl <- ifelse(dmeans$context %in% highlight, dmeans$context, "Other")
  col_map <- c(cmap, Other = "grey65")

  pC <- ggplot2::ggplot(dmeans, ggplot2::aes(x = x_mean, y = y_mean, color = context_hl)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    ggplot2::scale_color_manual(values = col_map, guide = "none") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = sprintf("Context means r = %.3f", r_means),
      x = sprintf("%s (mean)", feature_x),
      y = sprintf("%s (mean)", feature_y)
    )

  # Panel D: filtered overlay.
  df_hl <- df[df$context %in% highlight, , drop = FALSE]
  pD <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(color = "grey80", alpha = 0.25, size = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
    ggplot2::geom_point(
      data = df_hl,
      mapping = ggplot2::aes(color = context),
      alpha = 0.5,
      size = 0.8
    ) +
    ggplot2::geom_smooth(
      data = df_hl,
      mapping = ggplot2::aes(color = context),
      method = "lm",
      se = FALSE,
      linewidth = 0.6
    ) +
    ggplot2::scale_color_manual(values = cmap, guide = "none") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = sprintf("Filtered overlay (mean-resid r = %.3f)", r_resid),
      x = feature_x,
      y = feature_y
    )

  (pA | pB) / (pC | pD) +
    patchwork::plot_annotation(
      title = sprintf("%s vs %s (eps=%.2f, n_min=%d)", feature_x, feature_y, eps, n_min)
    )
}

sc_save_plot_dual <- function(plot_obj, out_base, width = 7, height = 5) {
  sc_require_pkg("ggplot2")
  ggplot2::ggsave(filename = paste0(out_base, ".pdf"), plot = plot_obj, width = width, height = height)
  ggplot2::ggsave(filename = paste0(out_base, ".png"), plot = plot_obj, width = width, height = height, dpi = 180)
}

sc_run_context_experiment <- function(
  experiment_id,
  seurat_obj,
  context,
  feature_x,
  feature_y = NULL,
  assay_x = "RNA",
  assay_y = NULL,
  slot_x = "data",
  slot_y = NULL,
  n_pairs = 200000,
  seed = 1,
  n_min = 200,
  method = "pearson",
  majority = "count",
  eps = 0.05,
  eps_grid = c(0, 0.05, 0.1),
  thresholds = seq(0, 0.5, by = 0.05),
  results_root,
  cache_root,
  top_n = 40,
  force = FALSE,
  robustness = list(
    equalize_context_sizes = FALSE,
    equalize_target_n = NULL,
    run_spearman = FALSE,
    spearman_n_pairs = 20000
  ),
  inputs_extra = list(),
  exemplar_plot_n = 3
) {
  sc_require_contextcorr()
  sc_require_pkg("Seurat")

  if (is.null(assay_y)) assay_y <- assay_x
  if (is.null(slot_y)) slot_y <- slot_x

  out_dir <- file.path(results_root, "sc", experiment_id)
  plot_dir <- file.path(out_dir, "plots")
  cache_dir <- file.path(cache_root, experiment_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(context) != ncol(seurat_obj)) {
    stop("`context` length must equal number of cells in `seurat_obj`.", call. = FALSE)
  }

  keep <- !is.na(context)
  if (!all(keep)) {
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[which(keep)])
    context <- context[keep]
  }

  if (isTRUE(robustness$equalize_context_sizes)) {
    eq_idx <- sc_equalize_context_sizes(context, seed = seed, target_n = robustness$equalize_target_n)
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[eq_idx])
    context <- context[eq_idx]
  }

  expr_x <- sc_get_expression_matrix(seurat_obj, features = feature_x, assay = assay_x, slot = slot_x)
  if (is.null(feature_y)) feature_y <- feature_x
  expr_y <- sc_get_expression_matrix(seurat_obj, features = feature_y, assay = assay_y, slot = slot_y)

  feature_x2 <- rownames(expr_x)
  feature_y2 <- rownames(expr_y)

  pair_file <- file.path(cache_dir, "pairs.rds")
  if (file.exists(pair_file) && !isTRUE(force)) {
    pairs <- readRDS(pair_file)
  } else {
    pairs <- sc_make_feature_pairs(feature_x2, feature_y2, n_pairs = n_pairs, seed = seed)
    saveRDS(pairs, pair_file)
  }

  stat_cache <- file.path(cache_dir, paste0("pair_stats_", method, ".rds"))
  stats <- sc_compute_pair_statistics(
    expr_x = expr_x,
    expr_y = expr_y,
    pairs = pairs,
    context = context,
    eps = eps,
    majority = majority,
    n_min = n_min,
    method = method,
    cache_file = stat_cache,
    force = force
  )

  summary_pairs <- stats$summary_pairs
  curve_df <- sc_reversal_rate_curve(
    r_global = summary_pairs$r_global,
    r_ctx_mat = stats$r_by_context,
    thresholds = thresholds,
    eps_grid = eps_grid
  )

  exemplars <- sc_rank_exemplars(summary_pairs, stats$r_by_context, stats$n_by_context, n_top = top_n)

  spearman_df <- NULL
  if (isTRUE(robustness$run_spearman)) {
    n_sp <- min(nrow(pairs), as.integer(robustness$spearman_n_pairs))
    sp_pairs <- pairs[seq_len(n_sp), , drop = FALSE]
    sp <- sc_compute_pair_statistics(
      expr_x = expr_x,
      expr_y = expr_y,
      pairs = sp_pairs,
      context = context,
      eps = eps,
      majority = majority,
      n_min = n_min,
      method = "spearman",
      cache_file = file.path(cache_dir, "pair_stats_spearman.rds"),
      force = force,
      progress_every = 500
    )
    spearman_df <- sp$summary_pairs
  }

  inputs_obj <- list(
    experiment_id = experiment_id,
    n_cells = ncol(seurat_obj),
    context_levels = sort(unique(as.character(context))),
    context_table = sort(table(context), decreasing = TRUE),
    assay_x = assay_x,
    assay_y = assay_y,
    slot_x = slot_x,
    slot_y = slot_y,
    feature_x = feature_x2,
    feature_y = feature_y2,
    n_pairs = nrow(pairs),
    params = list(
      n_min = n_min,
      method = method,
      majority = majority,
      eps = eps,
      eps_grid = eps_grid,
      thresholds = thresholds,
      seed = seed,
      robustness = robustness
    ),
    extra = inputs_extra
  )

  saveRDS(inputs_obj, file.path(out_dir, "inputs.rds"))
  utils::write.csv(summary_pairs, file.path(out_dir, "summary_pairs.csv"), row.names = FALSE)
  saveRDS(exemplars, file.path(out_dir, "top_exemplars.rds"))
  utils::write.csv(curve_df, file.path(out_dir, "reversal_rate_curve.csv"), row.names = FALSE)

  analysis_ready <- list(
    inputs = inputs_obj,
    summary_pairs = summary_pairs,
    curve = curve_df,
    pairs = pairs,
    r_by_context = stats$r_by_context,
    n_by_context = stats$n_by_context,
    context = as.character(context),
    context_levels = stats$context_levels,
    spearman_subset = spearman_df
  )
  saveRDS(analysis_ready, file.path(out_dir, "analysis_ready.rds"))

  p_curve <- sc_plot_reversal_curve(curve_df)
  p_i2 <- sc_plot_I2_distribution(summary_pairs)
  p_scatter <- sc_plot_global_vs_resid(summary_pairs)

  sc_save_plot_dual(p_curve, file.path(plot_dir, "reversal_rate_vs_threshold"), width = 7, height = 4.5)
  sc_save_plot_dual(p_i2, file.path(plot_dir, "I2_distribution_ecdf"), width = 6.5, height = 4.5)
  sc_save_plot_dual(p_scatter, file.path(plot_dir, "r_global_vs_r_resid"), width = 6.5, height = 5)

  n_ex_plot <- min(exemplar_plot_n, length(exemplars$indices))
  if (n_ex_plot > 0) {
    for (i in seq_len(n_ex_plot)) {
      idx <- exemplars$indices[i]
      pair_row <- summary_pairs[idx, , drop = FALSE]
      p_ex <- sc_plot_exemplar_pair(
        expr_x = expr_x,
        expr_y = expr_y,
        context = context,
        pair_row = pair_row,
        r_ctx = stats$r_by_context[idx, , drop = FALSE],
        n_ctx = stats$n_by_context[idx, , drop = FALSE],
        seed = seed + i
      )
      base <- file.path(plot_dir, sprintf("exemplar_%02d_%s__%s", i, pair_row$feature_x, pair_row$feature_y))
      sc_save_plot_dual(p_ex, base, width = 10, height = 4.5)
    }
  }

  if (!is.null(spearman_df)) {
    utils::write.csv(spearman_df, file.path(out_dir, "summary_pairs_spearman_subset.csv"), row.names = FALSE)
  }

  writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))

  list(
    out_dir = out_dir,
    inputs = inputs_obj,
    summary_pairs = summary_pairs,
    curve = curve_df,
    exemplars = exemplars,
    analysis_ready = analysis_ready
  )
}
