# Single-cell dataset adapters for context-correlation experiments.
#
# Contract: each loader returns a list with
# - seurat: Seurat object
# - dataset_id: stable dataset identifier
# - source: where data came from
# - metadata_fields: candidate metadata columns useful for context definitions
# - notes: short text notes

sc_require_paths <- function(paths = NULL) {
  if (!is.null(paths)) return(paths)
  if (!exists("paths", inherits = TRUE)) {
    stop("`paths` is not available. Run source('00-paths.R') before calling dataset loaders.", call. = FALSE)
  }
  get("paths", inherits = TRUE)
}

sc_inputs_dir <- function(paths = NULL) {
  paths <- sc_require_paths(paths)
  p <- file.path(paths$scratch, "sc_inputs")
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

sc_dataset_result <- function(seurat_obj, dataset_id, source, metadata_fields = character(), notes = NULL) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Dataset loader did not produce a Seurat object.", call. = FALSE)
  }
  list(
    seurat = seurat_obj,
    dataset_id = dataset_id,
    source = source,
    metadata_fields = metadata_fields,
    notes = notes,
    n_cells = ncol(seurat_obj),
    n_features = nrow(seurat_obj)
  )
}

sc_infer_metadata_fields <- function(seurat_obj) {
  cols <- colnames(seurat_obj@meta.data)
  pat <- c(
    "celltype", "cell_type", "annotation", "cluster", "seurat", "lineage",
    "subtype", "predicted", "l1", "l2", "tissue", "organ"
  )
  keep <- cols[vapply(cols, function(x) any(grepl(paste(pat, collapse = "|"), x, ignore.case = TRUE)), logical(1))]
  unique(keep)
}

sc_try_load_local_rds <- function(input_dir, filenames) {
  for (nm in filenames) {
    fp <- if (grepl("^/", nm)) nm else file.path(input_dir, nm)
    if (!file.exists(fp)) next
    obj <- tryCatch(readRDS(fp), error = function(e) e)
    if (inherits(obj, "Seurat")) {
      return(list(object = obj, source = paste0("local_rds:", normalizePath(fp))))
    }
    if (is.list(obj) && !is.null(obj$seurat) && inherits(obj$seurat, "Seurat")) {
      return(list(object = obj$seurat, source = paste0("local_rds_list:", normalizePath(fp))))
    }
  }
  NULL
}

sc_try_load_seuratdata <- function(dataset_candidates, auto_download = FALSE, verbose = TRUE) {
  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    return(NULL)
  }

  installed <- tryCatch(SeuratData::InstalledData(), error = function(e) NULL)
  installed_names <- character()
  if (is.data.frame(installed)) {
    installed_names <- unique(c(rownames(installed), installed$Dataset))
  }

  for (ds in dataset_candidates) {
    if (!(ds %in% installed_names) && auto_download) {
      if (isTRUE(verbose)) message("Installing SeuratData dataset: ", ds)
      try(SeuratData::InstallData(ds), silent = TRUE)
      installed <- tryCatch(SeuratData::InstalledData(), error = function(e) NULL)
      if (is.data.frame(installed)) {
        installed_names <- unique(c(rownames(installed), installed$Dataset))
      }
    }

    if (!(ds %in% installed_names)) next

    obj <- tryCatch(SeuratData::LoadData(ds), error = function(e) NULL)
    if (inherits(obj, "Seurat")) {
      return(list(object = obj, source = paste0("SeuratData:", ds)))
    }

    # Some datasets expose named entries.
    obj_alt <- tryCatch(SeuratData::LoadData(ds, type = ds), error = function(e) NULL)
    if (inherits(obj_alt, "Seurat")) {
      return(list(object = obj_alt, source = paste0("SeuratData:", ds, ":type")))
    }
  }

  NULL
}

sc_missing_dataset_error <- function(dataset_label, easy_sources, input_dir, expected_files, expected_assays = "RNA") {
  file_lines <- paste0("  - ", file.path(input_dir, expected_files), collapse = "\n")
  easy_lines <- paste0("  - ", easy_sources, collapse = "\n")
  stop(
    paste0(
      "Could not load ", dataset_label, ".\n",
      "Attempted easy sources:\n", easy_lines, "\n",
      "Place a Seurat object (.rds) at one of:\n", file_lines, "\n",
      "Expected assays: ", paste(expected_assays, collapse = ", "), "\n",
      "Data root for single-cell inputs: ", input_dir
    ),
    call. = FALSE
  )
}

load_pbmc_scrna <- function(paths = NULL, auto_download = FALSE, verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package `Seurat` is required.", call. = FALSE)
  }

  paths <- sc_require_paths(paths)
  input_dir <- sc_inputs_dir(paths)

  easy_candidates <- c("pbmc3k", "pbmcsca")
  local_files <- c(
    "pbmc_scrna.rds",
    "pbmc3k.rds",
    "pbmc_scRNA.rds",
    "pbmc_seurat.rds"
  )

  out <- sc_try_load_seuratdata(easy_candidates, auto_download = auto_download, verbose = verbose)
  if (is.null(out)) {
    out <- sc_try_load_local_rds(input_dir, local_files)
  }

  if (is.null(out)) {
    sc_missing_dataset_error(
      dataset_label = "PBMC scRNA-seq dataset",
      easy_sources = paste0("SeuratData::LoadData('", easy_candidates, "')"),
      input_dir = input_dir,
      expected_files = local_files,
      expected_assays = "RNA"
    )
  }

  seu <- out$object
  if (!("RNA" %in% names(seu@assays))) {
    stop("PBMC scRNA loader found data but assay `RNA` is missing.", call. = FALSE)
  }

  sc_dataset_result(
    seurat_obj = seu,
    dataset_id = "pbmc_scrna",
    source = out$source,
    metadata_fields = sc_infer_metadata_fields(seu),
    notes = "PBMC scRNA-seq canonical/teaching dataset"
  )
}

load_tabula_sapiens_immune <- function(
  paths = NULL,
  auto_download = FALSE,
  max_cells = 50000,
  seed = 1,
  tissue_include = NULL,
  celltype_include = NULL,
  verbose = TRUE
) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package `Seurat` is required.", call. = FALSE)
  }

  paths <- sc_require_paths(paths)
  input_dir <- sc_inputs_dir(paths)

  # There is no guaranteed SeuratData Tabula Sapiens immune dataset.
  # Still attempt common names first for convenience.
  easy_candidates <- c("tabulasapiens", "tabula.sapiens", "tabula_sapiens_immune")
  local_files <- c(
    "tabula_sapiens_immune.rds",
    "tabula_sapiens_blood_immune.rds",
    "tabula_sapiens.rds"
  )

  out <- sc_try_load_seuratdata(easy_candidates, auto_download = auto_download, verbose = verbose)
  if (is.null(out)) {
    out <- sc_try_load_local_rds(input_dir, local_files)
  }

  if (is.null(out)) {
    sc_missing_dataset_error(
      dataset_label = "Tabula Sapiens immune subset",
      easy_sources = c(
        paste0("SeuratData::LoadData('", easy_candidates, "')"),
        "local pre-converted Seurat .rds"
      ),
      input_dir = input_dir,
      expected_files = local_files,
      expected_assays = "RNA"
    )
  }

  seu <- out$object
  if (!("RNA" %in% names(seu@assays))) {
    stop("Tabula Sapiens loader found data but assay `RNA` is missing.", call. = FALSE)
  }

  md <- seu@meta.data

  if (!is.null(tissue_include)) {
    tissue_cols <- intersect(c("tissue", "organ_tissue", "organ", "tissue_general"), colnames(md))
    if (length(tissue_cols) > 0) {
      keep <- as.character(md[[tissue_cols[1]]]) %in% tissue_include
      seu <- subset(seu, cells = colnames(seu)[which(keep)])
      md <- seu@meta.data
    }
  }

  if (!is.null(celltype_include)) {
    ct_cols <- intersect(c("cell_type", "celltype", "cell_ontology_class", "celltype.l1", "annotation"), colnames(md))
    if (length(ct_cols) > 0) {
      keep <- as.character(md[[ct_cols[1]]]) %in% celltype_include
      seu <- subset(seu, cells = colnames(seu)[which(keep)])
    }
  }

  if (!is.null(max_cells) && ncol(seu) > max_cells) {
    set.seed(seed)
    keep_cells <- sample(colnames(seu), size = max_cells)
    seu <- subset(seu, cells = keep_cells)
  }

  sc_dataset_result(
    seurat_obj = seu,
    dataset_id = "tabula_sapiens_immune",
    source = out$source,
    metadata_fields = sc_infer_metadata_fields(seu),
    notes = "Tabula Sapiens immune/blood subset with optional downsampling"
  )
}

load_pbmc_citeseq <- function(paths = NULL, auto_download = FALSE, verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package `Seurat` is required.", call. = FALSE)
  }

  paths <- sc_require_paths(paths)
  input_dir <- sc_inputs_dir(paths)

  easy_candidates <- c("bmcite", "pbmcmultimodal")
  local_files <- c(
    "pbmc_citeseq.rds",
    "bmcite.rds",
    "pbmc_multimodal_citeseq.rds"
  )

  out <- sc_try_load_seuratdata(easy_candidates, auto_download = auto_download, verbose = verbose)
  if (is.null(out)) {
    out <- sc_try_load_local_rds(input_dir, local_files)
  }

  if (is.null(out)) {
    sc_missing_dataset_error(
      dataset_label = "PBMC CITE-seq dataset",
      easy_sources = paste0("SeuratData::LoadData('", easy_candidates, "')"),
      input_dir = input_dir,
      expected_files = local_files,
      expected_assays = c("RNA", "ADT")
    )
  }

  seu <- out$object
  assays <- names(seu@assays)
  if (!("RNA" %in% assays && "ADT" %in% assays)) {
    stop(
      paste0("PBMC CITE-seq loader requires both RNA and ADT assays; found assays: ", paste(assays, collapse = ", ")),
      call. = FALSE
    )
  }

  sc_dataset_result(
    seurat_obj = seu,
    dataset_id = "pbmc_citeseq",
    source = out$source,
    metadata_fields = sc_infer_metadata_fields(seu),
    notes = "PBMC CITE-seq dataset with RNA+ADT assays"
  )
}
