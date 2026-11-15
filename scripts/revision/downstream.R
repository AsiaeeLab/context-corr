#!/usr/bin/env Rscript

# Downstream network inference on pooled versus context-aware input matrices.
#
# Compare WGCNA and GENIE3 rankings under alternative correlation inputs,
# with STRING and TRRUST as exploratory external references.
# The analysis uses the 2,000 most variable mRNAs of the matched TCGA data
# (8,890 samples, 31 cancer types) with three inputs over the same genes and
# samples:
#   (i)   POOLED:        the raw log2(10 + RSEM) matrix
#   (ii)  RESIDUALIZED:  residualize_matrix(), per-cohort mean removed
#   (iii) RE-ADJACENCY:  random_effects_adjacency() r_RE, the DerSimonian-Laird
#                        pooled within-cohort correlation, used directly as the
#                        co-expression similarity
# Tools: WGCNA (soft threshold picked on the pooled matrix and reused, signed
# hybrid adjacency, TOM, average-linkage clustering, dynamicTreeCut modules) and
# GENIE3 (random-forest regulatory network with a fixed set of 100 TF
# regulators; GENIE3 needs an expression matrix, so it runs on (i) and (ii)
# only). External reference: STRING v12.0 human physical links with combined
# score >= 700 (a proxy for functional relatedness, not causal regulation), plus
# TRRUST v2 literature-curated TF-target pairs for the GENIE3 edges only.
#
# Run from the repo root:  Rscript scripts/revision/downstream.R
# Smoke test (small gene set, few trees, scratch output dir):
#   DOWNSTREAM_QUICK=1 DOWNSTREAM_OUT=/some/dir Rscript scripts/revision/downstream.R
#
# Downloads are cached under results/revision/downstream/cache/ and reused.
#
# Outputs (results/revision/downstream/):
#   genes.csv                   the 2,000 genes, variance rank, TF flag, STRING
#                               mapping, WGCNA module labels for the 3 inputs
#   soft_threshold.csv          pickSoftThreshold fit indices (pooled matrix)
#   edge_overlap.csv            ranking agreement between inputs at the
#                               similarity, WGCNA TOM, and GENIE3 levels
#   module_agreement.csv        adjusted Rand index between module assignments
#   edge_classification.csv     classify_pair labels for the top-10,000 edges
#                               of each input and a random subsample of pairs
#   edge_classification_summary.csv  counts and fractions per class and edge set
#   benchmark_precision.csv     STRING / TRRUST precision of top-k edges,
#                               random baselines, precision by pair class
#   Figure_Downstream.pdf       three-panel figure
#   edges_all.rds               all 1,999,000 pairs with r_pooled, r_resid,
#                               r_RE, I2, TOM weights, ranks, STRING support
#   genie3_weights.rds          GENIE3 weight matrices (pooled, residualized, and
#                               the pooled seed-replicate subset). Doubles as a cache:
#                               a rerun with identical parameters reuses it instead of
#                               recomputing (the two full runs take ~80 min). Delete it
#                               to force GENIE3 to run again.
#   runtime.txt                 wall times and sessionInfo()

t_start <- Sys.time()
section_times <- list()
tic <- function() Sys.time()
toc <- function(label, t0) {
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  section_times[[label]] <<- dt
  message(sprintf("[%s] %.1f s", label, dt))
}

source("00-paths.R")
source(file.path("R", "contextcorr.R"))
suppressPackageStartupMessages(library(data.table))
options(timeout = 900, stringsAsFactors = FALSE)

QUICK <- identical(Sys.getenv("DOWNSTREAM_QUICK"), "1")
out_dir <- Sys.getenv("DOWNSTREAM_OUT", unset = file.path("results", "revision", "downstream"))
cache_dir <- file.path("results", "revision", "downstream", "cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
n_cores <- max(1L, min(64L, parallel::detectCores() - 4L))

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
G <- if (QUICK) 300L else 2000L            # top mRNAs by global variance
k_top <- c(1000L, 5000L, 10000L)           # top-k edge sets for precision
k_small <- 1000L                           # "top-1,000 pooled edges ..."
k_big <- 10000L                            # "... outside the top-10,000 residualized edges"
string_min_score <- 700L                   # STRING combined score threshold
n_regulators <- if (QUICK) 20L else 100L   # GENIE3 regulators (most variable TFs)
genie3_ntrees <- if (QUICK) 50L else 500L  # GENIE3 default is 1000; 500 keeps the two full runs within ~1 h on 64 cores
genie3_replicate_targets <- if (QUICK) 60L else 500L   # seed-replicate control (pooled, second seed) on this many random targets
wgcna_type <- "signed hybrid"              # adjacency = r^beta for r > 0, else 0
wgcna_rsq_cut <- 0.85
wgcna_powers <- c(1:10, seq(12, 20, 2))
wgcna_fallback_power <- 12                 # WGCNA FAQ default for signed networks, n > 60
min_module_size <- 30L
deep_split <- 2L
n_scatter <- if (QUICK) 5000L else 30000L  # random pairs classified for panel A
n_random_draws <- 200L                     # Monte Carlo random-pair baselines
eps <- 0.05                                # classify_pair conventions (paper)
switch_thresh <- 0.3
alpha <- 0.05
seed_scatter <- 20260905L
seed_genie3 <- 20260906L
seed_genie3_replicate <- 20260916L
seed_random <- 20260907L

string_version <- "v12.0"
string_base <- "https://stringdb-downloads.org/download"
string_files <- c(
  links = "protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz",
  info = "protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
  aliases = "protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
)
tf_url <- "http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"
trrust_url <- "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"

fetch <- function(url, dest) {
  if (!file.exists(dest) || file.size(dest) == 0) {
    message("downloading ", url)
    ok <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE) == 0,
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) { unlink(dest); return(FALSE) }
  }
  TRUE
}
pkg_version <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not installed"
fmt <- function(x, d = 3) ifelse(is.finite(x), formatC(x, digits = d, format = "f"), "NA")
pct <- function(x, d = 1) ifelse(is.finite(x), paste0(formatC(100 * x, digits = d, format = "f"), "%"), "NA")
fnum <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
key_pair <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "|")
rank_desc <- function(x) {
  o <- order(x, decreasing = TRUE, na.last = TRUE)
  r <- integer(length(x)); r[o] <- seq_along(x); r
}

# ---------------------------------------------------------------------------
# 1. Data and gene set
# ---------------------------------------------------------------------------
t0 <- tic()
f1 <- file.path(paths$mirtcga, "matchedNoGBM.Rda")
f2 <- file.path(paths$mirtcga, "matchedNoGBM.rda")
if (file.exists(f1)) load(f1) else if (file.exists(f2)) load(f2) else stop("Cannot find matchedNoGBM.rda in ", paths$mirtcga)
stopifnot(ncol(mRNAdata) == length(cancerType), !anyNA(cancerType))
cohorts <- levels(cancerType)
N_c <- as.integer(table(cancerType)[cohorts]); names(N_c) <- cohorts
N_total <- ncol(mRNAdata)
k_ctx <- length(cohorts)

gene_var <- apply(mRNAdata, 1, var)
var_order <- order(gene_var, decreasing = TRUE)
genes <- rownames(mRNAdata)[var_order][seq_len(G)]
gene_var_sel <- gene_var[genes]
X_pooled <- mRNAdata[genes, , drop = FALSE]
rm(mRNAdata, miRdata); invisible(gc())
toc("load_data", t0)

# ---------------------------------------------------------------------------
# 2. The three inputs and the per-cohort correlations
# ---------------------------------------------------------------------------
t0 <- tic()
X_resid <- residualize_matrix(X_pooled, cancerType)
S_pooled <- stats::cor(t(X_pooled))
S_resid <- stats::cor(t(X_resid))
toc("pooled_resid_cor", t0)

t0 <- tic()
RE <- random_effects_adjacency(X_pooled, cancerType, min_n = 4)
S_RE <- RE$r_RE
toc("re_adjacency", t0)

t0 <- tic()
lt <- which(lower.tri(S_pooled), arr.ind = TRUE)      # pair m: gene i = lt[m, 1] (> j), gene j = lt[m, 2]
n_pairs <- nrow(lt)
r_ctx <- matrix(NA_real_, nrow = k_ctx, ncol = n_pairs, dimnames = list(cohorts, NULL))
for (ci in seq_len(k_ctx)) {
  idx <- which(cancerType == cohorts[ci])
  Rc <- suppressWarnings(stats::cor(t(X_pooled[, idx, drop = FALSE])))
  r_ctx[ci, ] <- Rc[lower.tri(Rc)]
}
rm(Rc)
toc("per_cohort_cor", t0)

edges <- data.table(
  pair = seq_len(n_pairs), i = lt[, 1], j = lt[, 2],
  gene_i = genes[lt[, 1]], gene_j = genes[lt[, 2]],
  r_pooled = S_pooled[lower.tri(S_pooled)],
  r_resid = S_resid[lower.tri(S_resid)],
  r_RE = S_RE[lower.tri(S_RE)],
  I2 = RE$I2[lower.tri(RE$I2)],
  k_RE = RE$k[lower.tri(RE$k)]
)
n_na_RE <- sum(is.na(edges$r_RE))
edges[, rank_pooled := rank_desc(abs(r_pooled))]
edges[, rank_resid := rank_desc(abs(r_resid))]
edges[, rank_RE := rank_desc(abs(r_RE))]

# ---------------------------------------------------------------------------
# 3. WGCNA: soft threshold on the pooled matrix, reused for all three inputs
# ---------------------------------------------------------------------------
t0 <- tic()
have_wgcna <- requireNamespace("WGCNA", quietly = TRUE) && requireNamespace("dynamicTreeCut", quietly = TRUE)
if (!have_wgcna) stop("WGCNA / dynamicTreeCut are not installed; install them (BiocManager::install('WGCNA')) and rerun.")
suppressPackageStartupMessages(library(WGCNA))
invisible(capture.output(suppressMessages(try(enableWGCNAThreads(nThreads = min(8L, n_cores)), silent = TRUE))))
sft <- pickSoftThreshold(t(X_pooled), networkType = wgcna_type, RsquaredCut = wgcna_rsq_cut,
                         powerVector = wgcna_powers, corFnc = cor, corOptions = list(use = "p"), verbose = 0)
beta <- sft$powerEstimate
if (is.na(beta)) {
  beta <- wgcna_fallback_power
  power_source <- sprintf("no power reached scale-free R^2 >= %.2f; fallback to the WGCNA FAQ default (%d) for signed networks with more than 60 samples",
                          wgcna_rsq_cut, wgcna_fallback_power)
} else {
  power_source <- sprintf("pickSoftThreshold powerEstimate on the pooled matrix (lowest power with scale-free R^2 >= %.2f, networkType = '%s')",
                          wgcna_rsq_cut, wgcna_type)
}
sft_tab <- sft$fitIndices
sft_tab$chosen <- sft_tab$Power == beta
utils::write.csv(sft_tab, file.path(out_dir, "soft_threshold.csv"), row.names = FALSE)

build_network <- function(S) {
  S <- as.matrix(S)
  S[is.na(S)] <- 0
  diag(S) <- 1
  A <- adjacency.fromSimilarity(S, type = wgcna_type, power = beta)
  TOM <- TOMsimilarity(A, TOMType = "unsigned", verbose = 0)   # adjacency is non-negative, so this equals the signed TOM
  dissTOM <- 1 - TOM
  tree <- stats::hclust(stats::as.dist(dissTOM), method = "average")
  mods <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = deep_split,
                                        pamRespectsDendro = FALSE, minClusterSize = min_module_size, verbose = 0)
  list(n_adj_nonzero = sum(A[lower.tri(A)] > 0), TOM_lt = TOM[lower.tri(TOM)], modules = as.integer(mods))
}
net <- list(pooled = build_network(S_pooled), resid = build_network(S_resid), RE = build_network(S_RE))
for (nm in names(net)) {
  edges[, paste0("TOM_", nm) := net[[nm]]$TOM_lt]
  edges[, paste0("rank_TOM_", nm) := rank_desc(net[[nm]]$TOM_lt)]
}
n_nonzero_adj <- vapply(net, function(z) z$n_adj_nonzero, numeric(1))
n_nonzero_TOM <- vapply(net, function(z) sum(z$TOM_lt > 0), numeric(1))
toc("wgcna", t0)

# ---------------------------------------------------------------------------
# 4. GENIE3 on (i) and (ii) with 100 fixed TF regulators
# ---------------------------------------------------------------------------
t0 <- tic()
tf_file <- file.path(cache_dir, "TF_names_v_1.01.txt")
have_tf <- fetch(tf_url, tf_file)
tf_names <- if (have_tf) readLines(tf_file) else character(0)
tfs_in_set <- genes[genes %in% tf_names]                 # genes are in decreasing-variance order
regulators <- head(tfs_in_set, n_regulators)
have_genie3 <- requireNamespace("GENIE3", quietly = TRUE) && requireNamespace("doRNG", quietly = TRUE) &&
  have_tf && length(regulators) >= min(n_regulators, 10L)
genie3_note <- if (!requireNamespace("GENIE3", quietly = TRUE)) "GENIE3 not installed" else
  if (!requireNamespace("doRNG", quietly = TRUE)) "doRNG not installed (needed for parallel GENIE3)" else
  if (!have_tf) "TF list download failed" else if (length(regulators) < min(n_regulators, 10L)) "too few TFs in the gene set" else "run"
gl <- NULL
gl_rep <- NULL
genie3_from_cache <- FALSE
genie3_compute_secs <- NA_real_
if (have_genie3) {
  run_genie3 <- function(M, seed, targets = NULL) {
    set.seed(seed)
    GENIE3::GENIE3(M, regulators = regulators, targets = targets, treeMethod = "RF", K = "sqrt", nTrees = genie3_ntrees,
                   nCores = n_cores, returnMatrix = TRUE, verbose = FALSE)
  }
  # The two full GENIE3 runs dominate the runtime, so a completed result is reused when every
  # parameter it depends on is unchanged. Delete genie3_weights.rds to force a recomputation.
  # The cached object records its own replicate target set, which is authoritative on reuse:
  # GENIE3 runs through doRNG and can leave the RNG in a different kind, so redrawing the
  # targets outside the original call order would not reproduce the same sample.
  cache_file <- file.path(out_dir, "genie3_weights.rds")
  cached <- if (file.exists(cache_file)) tryCatch(readRDS(cache_file), error = function(e) NULL) else NULL
  usable <- !is.null(cached) && identical(cached$regulators, regulators) && identical(cached$nTrees, genie3_ntrees) &&
    identical(cached$seed, seed_genie3) && identical(cached$seed_replicate, seed_genie3_replicate) &&
    setequal(rownames(cached$pooled), regulators) && setequal(colnames(cached$pooled), genes) &&
    identical(dimnames(cached$pooled), dimnames(cached$resid)) &&
    !is.null(cached$pooled_replicate_subset) && all(colnames(cached$pooled_replicate_subset) %in% genes)
  if (isTRUE(usable)) {
    W_pooled <- cached$pooled; W_resid <- cached$resid; W_pooled_rep <- cached$pooled_replicate_subset
    rep_targets <- colnames(W_pooled_rep)
    genie3_compute_secs <- if (is.null(cached$compute_secs)) NA_real_ else cached$compute_secs
    genie3_from_cache <- TRUE
    message("reusing cached GENIE3 weights from ", cache_file)
  } else {
    t_g3 <- Sys.time()
    W_pooled <- run_genie3(X_pooled, seed_genie3)
    W_resid <- run_genie3(X_resid, seed_genie3)
    stopifnot(identical(dimnames(W_pooled), dimnames(W_resid)))
    # seed-replicate control: the pooled matrix again with a different seed, on a random subset of targets,
    # so random-forest noise can be separated from the pooled-vs-residualized difference
    set.seed(seed_genie3_replicate)
    rep_targets <- sort(sample(genes, min(genie3_replicate_targets, G)))
    W_pooled_rep <- run_genie3(X_pooled, seed_genie3_replicate, targets = rep_targets)
    genie3_compute_secs <- as.numeric(difftime(Sys.time(), t_g3, units = "secs"))
    saveRDS(list(pooled = W_pooled, resid = W_resid, pooled_replicate_subset = W_pooled_rep, regulators = regulators,
                 nTrees = genie3_ntrees, K = "sqrt", seed = seed_genie3, seed_replicate = seed_genie3_replicate,
                 compute_secs = genie3_compute_secs),
            cache_file)
  }
  gl <- data.table(regulator = rep(rownames(W_pooled), times = ncol(W_pooled)),
                   target = rep(colnames(W_pooled), each = nrow(W_pooled)),
                   w_pooled = as.vector(W_pooled), w_resid = as.vector(W_resid))
  gl <- gl[regulator != target]
  gl[, rank_pooled := rank_desc(w_pooled)]
  gl[, rank_resid := rank_desc(w_resid)]
  Wp_sub <- W_pooled[, rep_targets, drop = FALSE]
  Wr_sub <- W_resid[, rep_targets, drop = FALSE]
  gl_rep <- data.table(regulator = rep(rownames(Wp_sub), times = ncol(Wp_sub)),
                       target = rep(colnames(Wp_sub), each = nrow(Wp_sub)),
                       w_pooled = as.vector(Wp_sub), w_resid = as.vector(Wr_sub), w_pooled_rep = as.vector(W_pooled_rep[, rep_targets]))
  gl_rep <- gl_rep[regulator != target]
  gl_rep[, rank_pooled := rank_desc(w_pooled)]
  gl_rep[, rank_resid := rank_desc(w_resid)]
  gl_rep[, rank_pooled_rep := rank_desc(w_pooled_rep)]
}
toc("genie3", t0)

# ---------------------------------------------------------------------------
# 5. External reference: STRING v12.0 physical links (and TRRUST for GENIE3)
# ---------------------------------------------------------------------------
t0 <- tic()
string_paths <- file.path(cache_dir, basename(string_files))
names(string_paths) <- names(string_files)
have_string <- all(mapply(function(u, d) fetch(u, d), file.path(string_base, string_files), string_paths))
if (!have_string) stop("STRING download failed; cannot benchmark.")
string_md5 <- tools::md5sum(string_paths)

info <- fread(string_paths["info"], sep = "\t", quote = "", header = TRUE)
setnames(info, 1, "string_protein_id")
ali <- fread(string_paths["aliases"], sep = "\t", quote = "", header = TRUE)
setnames(ali, 1, "string_protein_id")

ensp_of_gene <- setNames(rep(NA_character_, G), genes)
map_via <- setNames(rep("unmapped", G), genes)
pref <- setNames(info$string_protein_id, info$preferred_name)
direct <- genes %in% names(pref)
ensp_of_gene[direct] <- pref[genes[direct]]
map_via[direct] <- "preferred_name"
miss <- genes[!direct]
ali_hit <- ali[alias %in% miss, .(n_prot = uniqueN(string_protein_id), ensp = string_protein_id[1]), by = alias]
uniq <- ali_hit[n_prot == 1]
ensp_of_gene[uniq$alias] <- uniq$ensp
map_via[uniq$alias] <- "alias"
n_ambiguous_alias <- sum(ali_hit$n_prot > 1)
# a protein may serve only one gene in the set: keep the preferred_name match, drop alias-derived duplicates
dup <- !is.na(ensp_of_gene) & duplicated(ensp_of_gene) & map_via == "alias"
ensp_of_gene[dup] <- NA_character_; map_via[dup] <- "unmapped_duplicate_protein"
mapped <- !is.na(ensp_of_gene)
n_mapped_genes <- sum(mapped)
ensp2gene <- setNames(names(ensp_of_gene)[mapped], ensp_of_gene[mapped])
rm(ali); invisible(gc())

pl <- fread(string_paths["links"], header = TRUE)
n_links_rows <- nrow(pl)
pl <- pl[combined_score >= string_min_score]
n_links_rows_700 <- nrow(pl)
pl[, g1 := ensp2gene[protein1]]
pl[, g2 := ensp2gene[protein2]]
pl <- pl[!is.na(g1) & !is.na(g2) & g1 != g2]
ref_keys <- unique(key_pair(pl$g1, pl$g2))
rm(pl); invisible(gc())

edges[, mappable := mapped[i] & mapped[j]]
edges[, string_support := NA]
edges[mappable == TRUE, string_support := key_pair(gene_i, gene_j) %in% ref_keys]
n_mappable_pairs <- sum(edges$mappable)
n_ref_pairs <- sum(edges$string_support, na.rm = TRUE)
baseline_string <- n_ref_pairs / n_mappable_pairs

trrust_file <- file.path(cache_dir, "trrust_rawdata.human.tsv")
have_trrust <- fetch(trrust_url, trrust_file)
trrust_keys <- character(0)
if (have_trrust) {
  tr <- fread(trrust_file, header = FALSE)
  trrust_keys <- unique(paste(tr$V1, tr$V2, sep = ">"))
}
if (!is.null(gl)) {
  gl[, mappable := mapped[regulator] & mapped[target]]
  gl[, string_support := NA]
  gl[mappable == TRUE, string_support := key_pair(regulator, target) %in% ref_keys]
  gl[, trrust_support := if (have_trrust) paste(regulator, target, sep = ">") %in% trrust_keys else NA]
  baseline_genie3_string <- sum(gl$string_support, na.rm = TRUE) / sum(gl$mappable)
  baseline_genie3_trrust <- if (have_trrust) mean(gl$trrust_support) else NA_real_
}
toc("string_reference", t0)

# ---------------------------------------------------------------------------
# 6. Ranking agreement (edge_overlap.csv) and module agreement
# ---------------------------------------------------------------------------
t0 <- tic()
jaccard_topk <- function(ra, rb, k) { a <- ra <= k; b <- rb <= k; sum(a & b) / sum(a | b) }
compare_rankings <- function(level, name_a, name_b, x_a, x_b, ra, rb) {
  data.frame(
    level = level, input_a = name_a, input_b = name_b, n_edges = length(x_a),
    spearman_weights = suppressWarnings(stats::cor(x_a, x_b, method = "spearman", use = "complete.obs")),
    spearman_abs_weights = suppressWarnings(stats::cor(abs(x_a), abs(x_b), method = "spearman", use = "complete.obs")),
    jaccard_top1000 = jaccard_topk(ra, rb, 1000L),
    jaccard_top5000 = jaccard_topk(ra, rb, 5000L),
    jaccard_top10000 = jaccard_topk(ra, rb, 10000L),
    n_top1000_a_outside_top10000_b = sum(ra <= k_small & rb > k_big),
    n_top1000_b_outside_top10000_a = sum(rb <= k_small & ra > k_big),
    stringsAsFactors = FALSE
  )
}
ov <- list(
  compare_rankings("similarity_r", "pooled", "residualized", edges$r_pooled, edges$r_resid, edges$rank_pooled, edges$rank_resid),
  compare_rankings("similarity_r", "pooled", "RE_pooled", edges$r_pooled, edges$r_RE, edges$rank_pooled, edges$rank_RE),
  compare_rankings("similarity_r", "residualized", "RE_pooled", edges$r_resid, edges$r_RE, edges$rank_resid, edges$rank_RE),
  compare_rankings("WGCNA_TOM", "pooled", "residualized", edges$TOM_pooled, edges$TOM_resid, edges$rank_TOM_pooled, edges$rank_TOM_resid),
  compare_rankings("WGCNA_TOM", "pooled", "RE_pooled", edges$TOM_pooled, edges$TOM_RE, edges$rank_TOM_pooled, edges$rank_TOM_RE),
  compare_rankings("WGCNA_TOM", "residualized", "RE_pooled", edges$TOM_resid, edges$TOM_RE, edges$rank_TOM_resid, edges$rank_TOM_RE)
)
if (!is.null(gl)) {
  ov[[length(ov) + 1]] <- compare_rankings("GENIE3_weight", "pooled", "residualized", gl$w_pooled, gl$w_resid, gl$rank_pooled, gl$rank_resid)
  ov[[length(ov) + 1]] <- compare_rankings("GENIE3_weight_subset", "pooled", "pooled_second_seed", gl_rep$w_pooled, gl_rep$w_pooled_rep, gl_rep$rank_pooled, gl_rep$rank_pooled_rep)
  ov[[length(ov) + 1]] <- compare_rankings("GENIE3_weight_subset", "pooled", "residualized", gl_rep$w_pooled, gl_rep$w_resid, gl_rep$rank_pooled, gl_rep$rank_resid)
}
edge_overlap <- do.call(rbind, ov)
edge_overlap$ranking_basis <- ifelse(edge_overlap$level == "similarity_r", "|r| (top-k) and signed r (Spearman)",
                                     ifelse(edge_overlap$level == "WGCNA_TOM", "TOM weight (signed hybrid adjacency, power beta)",
                                            ifelse(edge_overlap$level == "GENIE3_weight", "GENIE3 importance weight, regulator -> target",
                                                   sprintf("GENIE3 importance weight on the %d-target seed-replicate subset", length(if (!is.null(gl_rep)) unique(gl_rep$target) else integer(0))))))
utils::write.csv(edge_overlap, file.path(out_dir, "edge_overlap.csv"), row.names = FALSE)

ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  comb2 <- function(x) x * (x - 1) / 2
  s_ij <- sum(comb2(tab)); s_a <- sum(comb2(rowSums(tab))); s_b <- sum(comb2(colSums(tab)))
  expected <- s_a * s_b / comb2(n)
  (s_ij - expected) / ((s_a + s_b) / 2 - expected)
}
mod <- data.frame(gene = genes, pooled = net$pooled$modules, resid = net$resid$modules, RE = net$RE$modules)
mod_pairs <- list(c("pooled", "resid"), c("pooled", "RE"), c("resid", "RE"))
module_agreement <- do.call(rbind, lapply(mod_pairs, function(p) {
  a <- mod[[p[1]]]; b <- mod[[p[2]]]
  both <- a != 0 & b != 0
  data.frame(
    input_a = ifelse(p[1] == "resid", "residualized", ifelse(p[1] == "RE", "RE_pooled", p[1])),
    input_b = ifelse(p[2] == "resid", "residualized", ifelse(p[2] == "RE", "RE_pooled", p[2])),
    n_genes = G,
    ari_all_genes_unassigned_as_module = ari(a, b),
    n_genes_assigned_in_both = sum(both),
    ari_assigned_in_both = if (sum(both) > 1) ari(a[both], b[both]) else NA_real_,
    n_modules_a = length(setdiff(unique(a), 0L)), n_modules_b = length(setdiff(unique(b), 0L)),
    n_unassigned_a = sum(a == 0L), n_unassigned_b = sum(b == 0L),
    largest_module_a = if (any(a != 0)) max(table(a[a != 0])) else 0L,
    largest_module_b = if (any(b != 0)) max(table(b[b != 0])) else 0L,
    stringsAsFactors = FALSE
  )
}))
module_agreement$power <- beta
module_agreement$network_type <- wgcna_type
module_agreement$min_module_size <- min_module_size
module_agreement$deep_split <- deep_split
utils::write.csv(module_agreement, file.path(out_dir, "module_agreement.csv"), row.names = FALSE)
toc("overlap_modules", t0)

# ---------------------------------------------------------------------------
# 7. classify_pair on the top-10,000 edges of each input and a random subsample
# ---------------------------------------------------------------------------
t0 <- tic()
set.seed(seed_scatter)
scatter_idx <- sort(sample.int(n_pairs, min(n_scatter, n_pairs)))
edges[, in_top10000_pooled := rank_pooled <= k_big]
edges[, in_top10000_resid := rank_resid <= k_big]
edges[, in_top10000_RE := rank_RE <= k_big]
edges[, in_random_subsample := FALSE]
edges[scatter_idx, in_random_subsample := TRUE]
to_classify <- edges[in_top10000_pooled | in_top10000_resid | in_top10000_RE | in_random_subsample, pair]

classify_one <- function(m) {
  cl <- classify_pair(edges$r_pooled[m], r_ctx[, m], N_c, eps = eps, switch_thresh = switch_thresh, alpha = alpha)
  data.table(pair = m, label = cl$label, k_informative = cl$k_informative, n_sig = cl$n_sig,
             n_sig_agree = cl$n_sig_agree, n_sig_disagree = cl$n_sig_disagree,
             n_ctx_pos_eps = cl$simpson$n_pos, n_ctx_neg_eps = cl$simpson$n_neg,
             n_switch_pos = cl$n_switch_pos, n_switch_neg = cl$n_switch_neg)
}
cls_list <- parallel::mclapply(to_classify, classify_one, mc.cores = n_cores)
failed <- vapply(cls_list, inherits, logical(1), what = "try-error")
if (any(failed)) stop("classify_pair failed for ", sum(failed), " pairs")
cls <- rbindlist(cls_list)
rm(cls_list)
edge_class <- merge(edges[pair %in% to_classify,
                          .(pair, gene_i, gene_j, r_pooled, r_resid, r_RE, I2, k_RE, rank_pooled, rank_resid, rank_RE,
                            TOM_pooled, TOM_resid, TOM_RE, in_top10000_pooled, in_top10000_resid, in_top10000_RE,
                            in_random_subsample, mappable, string_support)],
                    cls, by = "pair")
edge_class[, median_r_ctx := apply(r_ctx[, pair, drop = FALSE], 2, stats::median, na.rm = TRUE)]
setorder(edge_class, rank_pooled)
fwrite(edge_class, file.path(out_dir, "edge_classification.csv"))

class_levels <- c("globally_robust", "compositional", "context_switching", "indeterminate")
class_counts <- function(sel) {
  tab <- table(factor(edge_class$label[sel], levels = class_levels))
  data.frame(label = names(tab), n = as.integer(tab), fraction = as.numeric(tab) / sum(tab), stringsAsFactors = FALSE)
}
cc_pooled <- class_counts(edge_class$in_top10000_pooled)
cc_resid <- class_counts(edge_class$in_top10000_resid)
cc_RE <- class_counts(edge_class$in_top10000_RE)
cc_random <- class_counts(edge_class$in_random_subsample)
class_summary <- rbind(cbind(set = "top10000_pooled", cc_pooled), cbind(set = "top10000_residualized", cc_resid),
                       cbind(set = "top10000_RE_pooled", cc_RE), cbind(set = "random_subsample", cc_random))
utils::write.csv(class_summary, file.path(out_dir, "edge_classification_summary.csv"), row.names = FALSE)
# sign disagreement between pooled and residualized r among the top-10,000 pooled edges
top_pooled <- edge_class[in_top10000_pooled == TRUE]
n_sign_flip_resid <- sum(sign(top_pooled$r_pooled) != sign(top_pooled$r_resid))
n_sign_flip_RE <- sum(sign(top_pooled$r_pooled) != sign(top_pooled$r_RE), na.rm = TRUE)
min_abs_r_top_pooled <- min(abs(top_pooled$r_pooled))
# attenuation of the top-1,000 pooled edges once context is removed
lost_sel <- edges$rank_pooled <= k_small & edges$rank_resid > k_big
kept_sel <- edges$rank_pooled <= k_small & edges$rank_resid <= k_big
atten <- function(sel) c(n = sum(sel), med_abs_r_pooled = stats::median(abs(edges$r_pooled[sel])),
                         med_abs_r_resid = stats::median(abs(edges$r_resid[sel])),
                         med_I2 = stats::median(edges$I2[sel], na.rm = TRUE))
atten_lost <- atten(lost_sel); atten_kept <- atten(kept_sel)
toc("classification", t0)

# ---------------------------------------------------------------------------
# 8. Benchmark precision (benchmark_precision.csv)
# ---------------------------------------------------------------------------
t0 <- tic()
prec_row <- function(edge_set, ranking, subset, k, support, baseline, reference) {
  data.frame(edge_set = edge_set, ranking = ranking, subset = subset, reference = reference, k = k,
             n_edges = length(support), n_mappable = sum(!is.na(support)), n_supported = sum(support, na.rm = TRUE),
             precision_mappable = sum(support, na.rm = TRUE) / sum(!is.na(support)),
             precision_all_edges = sum(support, na.rm = TRUE) / length(support),
             baseline_precision = baseline,
             enrichment_over_baseline = (sum(support, na.rm = TRUE) / sum(!is.na(support))) / baseline,
             stringsAsFactors = FALSE)
}
bp <- list()
rank_cols <- c(pooled_abs_r = "rank_pooled", residualized_abs_r = "rank_resid", RE_pooled_abs_r = "rank_RE",
               pooled_TOM = "rank_TOM_pooled", residualized_TOM = "rank_TOM_resid", RE_pooled_TOM = "rank_TOM_RE")
for (rn in names(rank_cols)) for (k in k_top) {
  sel <- edges[[rank_cols[rn]]] <= k
  bp[[length(bp) + 1]] <- prec_row("gene_pairs", rn, "top_k", k, edges$string_support[sel], baseline_string, "STRING_physical_700")
}
# random-pair baselines: exact density plus Monte Carlo draws of k mappable pairs
bp[[length(bp) + 1]] <- data.frame(edge_set = "gene_pairs", ranking = "random_pairs_exact", subset = "all_mappable_pairs",
                                   reference = "STRING_physical_700", k = NA_integer_, n_edges = n_mappable_pairs,
                                   n_mappable = n_mappable_pairs, n_supported = n_ref_pairs, precision_mappable = baseline_string,
                                   precision_all_edges = n_ref_pairs / n_pairs, baseline_precision = baseline_string,
                                   enrichment_over_baseline = 1, stringsAsFactors = FALSE)
set.seed(seed_random)
mappable_idx <- which(edges$mappable)
mc_rows <- list()
for (k in k_top) {
  draws <- replicate(n_random_draws, mean(edges$string_support[sample(mappable_idx, k)]))
  mc_rows[[length(mc_rows) + 1]] <- data.frame(edge_set = "gene_pairs", ranking = "random_pairs_monte_carlo", subset = sprintf("%d draws", n_random_draws),
                                               reference = "STRING_physical_700", k = k, n_edges = k, n_mappable = k,
                                               n_supported = round(mean(draws) * k, 1), precision_mappable = mean(draws),
                                               precision_all_edges = mean(draws), baseline_precision = baseline_string,
                                               enrichment_over_baseline = mean(draws) / baseline_string, stringsAsFactors = FALSE)
  mc_rows[[length(mc_rows)]]$mc_q025 <- as.numeric(stats::quantile(draws, 0.025))
  mc_rows[[length(mc_rows)]]$mc_q975 <- as.numeric(stats::quantile(draws, 0.975))
}
# precision by classify_pair class among the top-10,000 pooled edges
for (lab in class_levels) {
  sel <- top_pooled$label == lab
  if (sum(sel) > 0) bp[[length(bp) + 1]] <- prec_row("gene_pairs", "pooled_abs_r", paste0("top10000_class_", lab), k_big, top_pooled$string_support[sel], baseline_string, "STRING_physical_700")
}
bp[[length(bp) + 1]] <- prec_row("gene_pairs", "pooled_abs_r", "top10000_all_classes", k_big, top_pooled$string_support, baseline_string, "STRING_physical_700")
# edges lost / kept when context is removed (top-1,000 pooled vs top-10,000 residualized or RE), and the converse
lost_sets <- list(
  list("pooled_abs_r", "top1000_lost_outside_top10000_residualized", edges$rank_pooled <= k_small & edges$rank_resid > k_big),
  list("pooled_abs_r", "top1000_kept_inside_top10000_residualized", edges$rank_pooled <= k_small & edges$rank_resid <= k_big),
  list("pooled_abs_r", "top1000_lost_outside_top10000_RE_pooled", edges$rank_pooled <= k_small & edges$rank_RE > k_big),
  list("pooled_abs_r", "top1000_kept_inside_top10000_RE_pooled", edges$rank_pooled <= k_small & edges$rank_RE <= k_big),
  list("residualized_abs_r", "top1000_new_outside_top10000_pooled", edges$rank_resid <= k_small & edges$rank_pooled > k_big),
  list("residualized_abs_r", "top1000_shared_inside_top10000_pooled", edges$rank_resid <= k_small & edges$rank_pooled <= k_big),
  list("RE_pooled_abs_r", "top1000_new_outside_top10000_pooled", edges$rank_RE <= k_small & edges$rank_pooled > k_big),
  list("RE_pooled_abs_r", "top1000_shared_inside_top10000_pooled", edges$rank_RE <= k_small & edges$rank_pooled <= k_big),
  list("pooled_TOM", "top1000_lost_outside_top10000_residualized", edges$rank_TOM_pooled <= k_small & edges$rank_TOM_resid > k_big),
  list("pooled_TOM", "top1000_kept_inside_top10000_residualized", edges$rank_TOM_pooled <= k_small & edges$rank_TOM_resid <= k_big)
)
for (ls in lost_sets) if (sum(ls[[3]]) > 0) bp[[length(bp) + 1]] <- prec_row("gene_pairs", ls[[1]], ls[[2]], k_small, edges$string_support[ls[[3]]], baseline_string, "STRING_physical_700")
# GENIE3 edges
if (!is.null(gl)) {
  for (rn in c("pooled", "residualized")) for (k in k_top) {
    sel <- (if (rn == "pooled") gl$rank_pooled else gl$rank_resid) <= k
    bp[[length(bp) + 1]] <- prec_row("GENIE3_regulator_target", paste0(rn, "_GENIE3_weight"), "top_k", k, gl$string_support[sel], baseline_genie3_string, "STRING_physical_700")
    if (have_trrust) bp[[length(bp) + 1]] <- prec_row("GENIE3_regulator_target", paste0(rn, "_GENIE3_weight"), "top_k", k, gl$trrust_support[sel], baseline_genie3_trrust, "TRRUST_v2")
  }
  bp[[length(bp) + 1]] <- data.frame(edge_set = "GENIE3_regulator_target", ranking = "random_pairs_exact", subset = "all_mappable_regulator_target_pairs",
                                     reference = "STRING_physical_700", k = NA_integer_, n_edges = sum(gl$mappable), n_mappable = sum(gl$mappable),
                                     n_supported = sum(gl$string_support, na.rm = TRUE), precision_mappable = baseline_genie3_string,
                                     precision_all_edges = sum(gl$string_support, na.rm = TRUE) / nrow(gl), baseline_precision = baseline_genie3_string,
                                     enrichment_over_baseline = 1, stringsAsFactors = FALSE)
  if (have_trrust) bp[[length(bp) + 1]] <- data.frame(edge_set = "GENIE3_regulator_target", ranking = "random_pairs_exact", subset = "all_regulator_target_pairs",
                                                      reference = "TRRUST_v2", k = NA_integer_, n_edges = nrow(gl), n_mappable = nrow(gl),
                                                      n_supported = sum(gl$trrust_support), precision_mappable = baseline_genie3_trrust,
                                                      precision_all_edges = baseline_genie3_trrust, baseline_precision = baseline_genie3_trrust,
                                                      enrichment_over_baseline = 1, stringsAsFactors = FALSE)
}
benchmark <- do.call(rbind, bp)
benchmark$mc_q025 <- NA_real_; benchmark$mc_q975 <- NA_real_
benchmark <- rbind(benchmark, do.call(rbind, mc_rows))
rownames(benchmark) <- NULL
benchmark$row <- seq_len(nrow(benchmark))
utils::write.csv(benchmark, file.path(out_dir, "benchmark_precision.csv"), row.names = FALSE)
toc("benchmark", t0)

# gene table
gene_tab <- data.frame(gene = genes, variance_rank = seq_len(G), variance = as.numeric(gene_var_sel),
                       is_TF_Lambert2018 = genes %in% tf_names, genie3_regulator = genes %in% regulators,
                       string_protein_id = as.character(ensp_of_gene), string_mapped_via = as.character(map_via),
                       module_pooled = net$pooled$modules, module_residualized = net$resid$modules, module_RE_pooled = net$RE$modules,
                       stringsAsFactors = FALSE)
utils::write.csv(gene_tab, file.path(out_dir, "genes.csv"), row.names = FALSE)
saveRDS(edges, file.path(out_dir, "edges_all.rds"))

# ---------------------------------------------------------------------------
# 9. Figure
# ---------------------------------------------------------------------------
t0 <- tic()
class_cols <- c(globally_robust = "#0072B2", compositional = "#D55E00", context_switching = "#CC79A7", indeterminate = "grey72")
input_cols <- c(pooled = "#E69F00", residualized = "#009E73", RE_pooled = "#7570B3")
# Panel letter at the top left of the panel. adj = 0 with no `at` anchors it to the left edge
# of the plot region, which is safe on a log axis and in a panel of any width (an `at` computed
# as a fraction of the user range lands outside the device in a wide panel and is clipped).
panel_label <- function(lab) mtext(lab, side = 3, line = 1.2, adj = 0, font = 2, cex = 1.3)

pdf_file <- file.path(out_dir, "Figure_Downstream.pdf")
# Two rows so the figure keeps a legible panel size when scaled to a portrait
# text block: A and B on top, C spanning the bottom.
fig_w <- 10.0; fig_h <- 8.4
grDevices::cairo_pdf(pdf_file, width = fig_w, height = fig_h)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.95))
par(mar = c(4.6, 4.8, 3.4, 1.2), oma = c(0, 0.6, 0, 0), mgp = c(2.8, 0.7, 0))

# A. pooled vs residualized r for the random subsample, colored by class
sc <- edge_class[in_random_subsample == TRUE]
sc <- sc[order(factor(sc$label, levels = c("indeterminate", "globally_robust", "compositional", "context_switching")))]
lim <- range(c(sc$r_pooled, sc$r_resid), na.rm = TRUE); lim <- c(-1, 1) * max(abs(lim))
plot(sc$r_pooled, sc$r_resid, pch = 16, cex = 0.45, col = grDevices::adjustcolor(class_cols[sc$label], alpha.f = 0.55),
     xlim = lim, ylim = lim, xlab = bquote(r[pooled] ~ .(sprintf("(all %s samples)", fnum(N_total)))),
     ylab = expression(r[resid] ~ "(cohort means removed)"), main = "", las = 1)
abline(h = 0, v = 0, col = "grey60", lty = 3)
abline(0, 1, col = "grey40", lty = 2)
title(main = "Pooled vs Residualized Edge Weights", font.main = 2, cex.main = 1.1)
panel_label("A.")
cnt <- table(factor(sc$label, levels = class_levels))
legend("topleft", bty = "n", cex = 0.85, pch = 16, pt.cex = 1.1, col = class_cols[class_levels],
       legend = sprintf("%s (n = %s)", c("Globally robust", "Compositional", "Context switching", "Indeterminate"), fnum(as.integer(cnt))))

# B. STRING precision of the top-k edges
prec_curve <- function(rn) { b <- benchmark[benchmark$ranking == rn & benchmark$subset == "top_k" & benchmark$edge_set == "gene_pairs", ]; b[match(k_top, b$k), "precision_mappable"] }
curves <- list(pooled = prec_curve("pooled_abs_r"), residualized = prec_curve("residualized_abs_r"), RE_pooled = prec_curve("RE_pooled_abs_r"),
               pooled_TOM = prec_curve("pooled_TOM"), residualized_TOM = prec_curve("residualized_TOM"), RE_pooled_TOM = prec_curve("RE_pooled_TOM"))
ymax <- max(unlist(curves), baseline_string, na.rm = TRUE) * 1.15
plot(NA, xlim = c(800, 12500), ylim = c(0, ymax), log = "x", xaxt = "n", las = 1,
     xlab = expression("Top-" * italic(k) * " edges by weight"), ylab = "Precision (fraction with STRING physical link)", main = "")
axis(1, at = k_top, labels = fnum(k_top))
abline(h = baseline_string, col = "grey35", lty = 3, lwd = 1.6)
for (nm in c("pooled", "residualized", "RE_pooled")) {
  lines(k_top, curves[[nm]], col = input_cols[nm], lwd = 2.4); points(k_top, curves[[nm]], col = input_cols[nm], pch = 16, cex = 1.2)
  lines(k_top, curves[[paste0(nm, "_TOM")]], col = input_cols[nm], lwd = 2, lty = 2); points(k_top, curves[[paste0(nm, "_TOM")]], col = input_cols[nm], pch = 17, cex = 1.1)
}
title(main = "STRING Support of Top-k Edges", font.main = 2, cex.main = 1.1)
panel_label("B.")
legend("topright", bty = "n", cex = 0.82, ncol = 1,
       legend = c("Pooled", "Residualized", expression("RE-pooled " * r[RE]), expression("Ranking by |" * r * "|"), "Ranking by WGCNA TOM", "Random gene pairs"),
       col = c(input_cols, "grey20", "grey20", "grey35"), lwd = c(2.4, 2.4, 2.4, 2.4, 2, 1.6), lty = c(1, 1, 1, 1, 2, 3),
       pch = c(16, 16, 16, 16, 17, NA), pt.cex = c(1.2, 1.2, 1.2, 1.2, 1.1, NA))

# C. agreement between pooled and context-aware outputs
get_ov <- function(level, b, col) { r <- edge_overlap[edge_overlap$level == level & edge_overlap$input_a == "pooled" & edge_overlap$input_b == b, ]; if (nrow(r) == 0) NA_real_ else r[[col]] }
get_ari <- function(b) module_agreement$ari_all_genes_unassigned_as_module[module_agreement$input_a == "pooled" & module_agreement$input_b == b]
bar_mat <- rbind(
  residualized = c(get_ov("similarity_r", "residualized", "jaccard_top1000"), get_ov("similarity_r", "residualized", "jaccard_top10000"),
                   get_ov("WGCNA_TOM", "residualized", "jaccard_top1000"), get_ov("WGCNA_TOM", "residualized", "jaccard_top10000"),
                   get_ari("residualized"), get_ov("GENIE3_weight", "residualized", "jaccard_top1000"), get_ov("GENIE3_weight", "residualized", "jaccard_top10000")),
  RE_pooled = c(get_ov("similarity_r", "RE_pooled", "jaccard_top1000"), get_ov("similarity_r", "RE_pooled", "jaccard_top10000"),
                get_ov("WGCNA_TOM", "RE_pooled", "jaccard_top1000"), get_ov("WGCNA_TOM", "RE_pooled", "jaccard_top10000"),
                get_ari("RE_pooled"), NA, NA)
)
colnames(bar_mat) <- c("|r|\ntop 1,000", "|r|\ntop 10,000", "TOM\ntop 1,000", "TOM\ntop 10,000", "Modules\n(ARI)", "GENIE3\ntop 1,000", "GENIE3\ntop 10,000")
par(mar = c(5.2, 4.8, 3.4, 1.2))
bp_x <- barplot(bar_mat, beside = TRUE, col = input_cols[c("residualized", "RE_pooled")], border = NA, ylim = c(0, 1), las = 1,
                names.arg = rep("", ncol(bar_mat)), ylab = "Agreement with the pooled output", main = "")
axis(1, at = colMeans(bp_x), labels = FALSE, tick = FALSE)
mtext(colnames(bar_mat), side = 1, at = colMeans(bp_x), line = 2.4, cex = 0.8)
title(main = "Agreement of Pooled and Context-Aware Outputs", font.main = 2, cex.main = 1.1)
panel_label("C.")
# GENIE3 seed-replicate ceiling: the same pooled matrix rerun with a different seed, on the target subset
ceil <- if (!is.null(gl_rep)) c(get_ov("GENIE3_weight_subset", "pooled_second_seed", "jaccard_top1000"),
                                get_ov("GENIE3_weight_subset", "pooled_second_seed", "jaccard_top10000")) else c(NA, NA)
if (all(is.finite(ceil))) {
  for (m in 1:2) {
    xg <- bp_x[, ncol(bar_mat) - 2 + m]
    segments(min(xg) - 0.45, ceil[m], max(xg) + 0.45, ceil[m], lwd = 2.2, lty = 1, col = "grey20")
  }
}
# Stack the legend at the left: the seed-replicate segments sit at y > 0.8 over the
# GENIE3 bars on the right, and a wide multi-column legend runs straight through them.
legend("topleft", bty = "n", cex = 0.85, ncol = 1, inset = c(0.01, 0),
       fill = c(input_cols[c("residualized", "RE_pooled")], NA), border = NA,
       legend = c("Pooled vs residualized", expression("Pooled vs RE-pooled " * r[RE]), "GENIE3 same-input seed replicate"),
       lty = c(NA, NA, 1), lwd = c(NA, NA, 2.2), col = c(NA, NA, "grey20"), merge = FALSE)
invisible(dev.off())
toc("figure", t0)

# Runtime and environment information
t0 <- tic()
total_s <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
rt <- c(sprintf("downstream.R finished %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        sprintf("QUICK mode: %s", QUICK), sprintf("cores: %d", n_cores), sprintf("G = %d genes, %d pairs", G, n_pairs),
        sprintf("WGCNA power beta = %s (%s)", beta, power_source),
        sprintf("GENIE3: %s, nTrees = %d, regulators = %d, weights %s", genie3_note, genie3_ntrees, length(regulators),
                if (genie3_from_cache) "reused from cached genie3_weights.rds" else "computed in this run"),
        sprintf("total wall time: %.1f s (%.1f min)", total_s, total_s / 60),
        "section wall times (s):", sprintf("  %-20s %8.1f", names(section_times), unlist(section_times)),
        "", capture.output(print(sessionInfo())))
writeLines(rt, file.path(out_dir, "runtime.txt"))
toc("report", t0)
message(sprintf("Done in %.1f min. Outputs in %s", total_s / 60, out_dir))
