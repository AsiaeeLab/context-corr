#!/usr/bin/env Rscript

# Composition-adjusted reversal and heterogeneity rates for the TCGA
# miRNA-mRNA analysis using published cell-fraction estimates, tumor purity,
# ESTIMATE scores, and marker-based covariates.
#
# Run from the repo root:  Rscript scripts/revision/composition.R
# Smoke test (small domain, few bootstrap replicates, scratch output dir):
#   COMP_QUICK=1 COMP_OUT=/some/dir Rscript scripts/revision/composition.R
#
# Question answered: how much of the TCGA Simpson-reversal rate and of the
# between-cohort heterogeneity (I^2) survives adjustment for per-sample
# cell-type composition?
#
# Composition covariates (per sample):
#   Published (Thorsson et al. 2018, Immunity 48:812, TCGA PanCanAtlas immune
#   landscape; files hosted by the GDC at
#   https://gdc.cancer.gov/about-data/publications/panimmune):
#     - leukocyte fraction (LF), DNA-methylation based
#       (TCGA_all_leuk_estimate.masked.20170107.tsv)
#     - CIBERSORT relative fractions of 22 leukocyte subsets, run on the
#       PanCanAtlas Kallisto RNA-seq (TCGA.Kallisto.fullIDs.cibersort.relative.tsv);
#       absolute fraction = relative fraction x LF, as in Thorsson et al.
#   Published (TCGA PanCanAtlas, https://gdc.cancer.gov/about-data/publications/pancanatlas):
#     - ABSOLUTE tumor purity (TCGA_mastercalls.abs_tables_JSedit.fixed.txt)
#   Computed from the expression matrix analyzed in the paper:
#     - ESTIMATE stromal and immune scores (Yoshihara et al. 2013), computed
#       with the algorithm of the `estimate` R package (v1.0.13, R-Forge),
#       platform = "illumina"; the in-memory implementation below is verified
#       against estimate::estimateScore on a subset of samples at run time.
#     - marker-based immune score (mean standardized log expression of PTPRC,
#       CD3E, CD3D, CD2, CD19, MS4A1, CD14, CD68, FCGR3A, NKG7, GZMB), stromal
#       score (COL1A1, COL1A2, COL3A1, DCN, LUM, FAP, PDGFRB, ACTA2), and a
#       purity proxy = -(immune + stromal). This is a marker-based proxy, not
#       a deconvolution; it is included so that the analysis has a fallback
#       that needs no external download, and its agreement with the
#       published estimates is reported.
#
# Adjustment modes (each feature is regressed on the covariates and replaced
# by the part not explained by them):
#   pooled          one regression over all samples, ignoring cohort (what an
#                   analyst gets by adding cell-type fractions as covariates
#                   without stratifying); cohort means change only through
#                   the covariates.
#   within          cohort-specific regressions; residuals plus the cohort
#                   mean, so cohort means are retained and only the
#                   within-cohort composition component is removed. Within-
#                   cohort r_c become partial correlations given composition.
#   within_between  `within` plus a cohort-level regression of the cohort
#                   means on the cohort-mean covariates (weighted by N_c), so
#                   the between-cohort mean shifts that covary with cohort
#                   composition are removed as well. This regression has only k
#                   cohort-level observations, so with many covariates it is
#                   close to saturated and removes between-cohort variance by
#                   construction; the residual degrees of freedom k - 1 - p are
#                   reported as n_cohort_level_residual_df.
#   within_resid    literal within-cohort residuals (cohort means removed).
#                   Its r_c equal those of `within`; its pooled r_global is
#                   the pooled within-cohort partial correlation, so its
#                   reversal rate is not comparable to the headline rate
#                   (the between-cohort component is gone by construction).
#                   Reported for completeness.
#
# Paper conventions reused (scripts/robustness_signflip.R,
# scripts/heterogeneity_I2_FDR.R, R/contextcorr.R):
#   - high-variance domain: top 5,000 mRNAs x top 200 miRNAs by global variance
#     (10^6 pairs); no pair subsampling here, the full domain is used
#   - r_c valid when finite and |r_c| < 1; pair evaluable when >= 10 valid cohorts
#   - sign tolerance eps: s_c = +1 if r_c >= eps, -1 if r_c <= -eps, else 0
#   - mixed-sign: both signs present; Simpson reversal: sign(r_global) !=
#     sign(sum_c s_c), ties and r_global = 0 excluded; reported among
#     |r_global| >= 0.2
#   - Cochran Q on Fisher z with weights n_c - 3; I^2 = max(0, (Q - (k-1))/Q).
#     For within-cohort adjustments with p covariates the weights are
#     n_c - p - 3 (the Fisher-z variance of a partial correlation).
#   - residualized flip: sign(r_global) != sign(r_resid), where r_resid is the
#     correlation after subtracting cohort means from both features
#
# Sample blocks: "all" = the paper's 8,890 samples and 31 cohorts (marker and
# ESTIMATE covariates exist for every sample); "matched" = samples with LF,
# CIBERSORT and ABSOLUTE purity all available (tumor and metastasis samples
# only; DLBC and THYM have no published estimates and drop out). Every
# adjusted rate is reported next to the unadjusted rate computed on the same
# samples and the same 10^6 pairs.
#
# Outputs (results/revision/composition/):
#   covariates_per_sample.csv   all covariates for the 8,890 samples
#   covariates_summary.csv      source, URL, md5, match rate, distribution
#   covariate_correlations.csv  agreement between covariate sources
#   adjusted_rates.csv          one row per (block, covariate set, mode) with
#                               unadjusted-vs-adjusted rates and bootstrap CIs
#   reversal_transitions.csv    persist / vanish / newly appear, per condition
#   meanshift_shrinkage.csv     between-cohort variance of cohort means after
#                               vs before adjustment, per condition
#   exemplars.csv               miR-200c-3p vs HTRA3 and let-7a-5p vs canonical
#                               targets, per condition
#   pair_level_key_conditions.rds  per-pair r_global / I^2 / reversal flags for
#                               the baselines and the headline conditions
#   Figure_CompositionAdjusted.pdf
#   runtime.txt

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

QUICK <- identical(Sys.getenv("COMP_QUICK"), "1")
out_dir <- Sys.getenv("COMP_OUT", unset = file.path("results", "revision", "composition"))
raw_dir <- Sys.getenv("COMP_RAW", unset = file.path(paths$raw, "panimmune"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

n_cores <- max(1L, min(24L, parallel::detectCores() - 4L))

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
G <- if (QUICK) 400 else 5000       # top mRNAs by global variance (paper domain)
M <- if (QUICK) 40 else 200         # top miRNAs by global variance (paper domain)
min_n_per_cohort <- 20              # paper convention
min_k <- 10                         # paper convention
eps_main <- 0.05
eps_grid <- c(0, 0.05, 0.10)
t_main <- 0.2
B <- if (QUICK) 25L else 500L       # pair-level bootstrap replicates
seed_boot <- 20260905
seed_fig <- 20260906

immune_markers <- c("PTPRC", "CD3E", "CD3D", "CD2", "CD19", "MS4A1", "CD14", "CD68", "FCGR3A", "NKG7", "GZMB")
stromal_markers <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRB", "ACTA2")

exemplar_pairs <- rbind(
  data.frame(set = "HTRA3", mir = "hsa-miR-200c-3p", gene = "HTRA3", stringsAsFactors = FALSE),
  data.frame(set = "let-7a", mir = "hsa-let-7a-5p",
             gene = c("KRAS", "NRAS", "HRAS", "HMGA2", "MYC", "LIN28A", "LIN28B", "IGF2BP1", "IGF2BP2"),
             stringsAsFactors = FALSE)
)

# ---------------------------------------------------------------------------
# External covariate files (download if missing, record md5)
# ---------------------------------------------------------------------------
ext_files <- list(
  leuk = list(
    name = "TCGA_all_leuk_estimate.masked.20170107.tsv",
    url = "https://api.gdc.cancer.gov/data/6f75c9d7-5134-4ed1-b8f3-72856c98a4e8",
    page = "https://gdc.cancer.gov/about-data/publications/panimmune",
    source = "Thorsson et al. 2018 (PanCanAtlas immune landscape), DNA-methylation leukocyte fraction"),
  cibersort = list(
    name = "TCGA.Kallisto.fullIDs.cibersort.relative.tsv",
    url = "https://api.gdc.cancer.gov/data/b3df502e-3594-46ef-9f94-d041a20a0b9a",
    page = "https://gdc.cancer.gov/about-data/publications/panimmune",
    source = "Thorsson et al. 2018 (PanCanAtlas immune landscape), CIBERSORT relative fractions (22 subsets)"),
  absolute = list(
    name = "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
    url = "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
    page = "https://gdc.cancer.gov/about-data/publications/pancanatlas",
    source = "TCGA PanCanAtlas ABSOLUTE purity/ploidy calls (Aran et al. 2015 method; Taylor et al. 2018 release)")
)
t0 <- tic()
for (nm in names(ext_files)) {
  f <- file.path(raw_dir, ext_files[[nm]]$name)
  if (!file.exists(f) || file.size(f) < 1000) {
    message("Downloading ", ext_files[[nm]]$name)
    ok <- FALSE
    for (attempt in 1:6) {
      st <- tryCatch(utils::download.file(ext_files[[nm]]$url, f, method = "curl", quiet = TRUE,
                                          extra = "-L --http1.1 --retry 5 --retry-delay 3"),
                     error = function(e) 1L, warning = function(w) 1L)
      if (identical(st, 0L) && file.exists(f) && file.size(f) > 1000) { ok <- TRUE; break }
      Sys.sleep(3)
    }
    if (!ok) stop("Could not download ", ext_files[[nm]]$name, " from ", ext_files[[nm]]$url)
  }
  ext_files[[nm]]$path <- normalizePath(f)
  ext_files[[nm]]$md5 <- unname(tools::md5sum(f))
  ext_files[[nm]]$bytes <- file.size(f)
}
toc("external_files", t0)

# ---------------------------------------------------------------------------
# Load matched TCGA data
# ---------------------------------------------------------------------------
t0 <- tic()
f1 <- file.path(paths$mirtcga, "matchedNoGBM.Rda")
f2 <- file.path(paths$mirtcga, "matchedNoGBM.rda")
if (file.exists(f1)) {
  load(f1)
} else if (file.exists(f2)) {
  load(f2)
} else {
  stop("Cannot find matchedNoGBM.Rda or matchedNoGBM.rda in ", paths$mirtcga)
}
stopifnot(identical(colnames(miRdata), colnames(mRNAdata)))
N_all <- ncol(mRNAdata)
sample_id <- colnames(mRNAdata)                       # e.g. TCGA.OR.A5J1.01
key15 <- gsub(".", "-", sample_id, fixed = TRUE)      # TCGA-OR-A5J1-01
sample_type <- as.character(sampleType)
cohort_all <- as.character(cancerType)

gene_var <- apply(mRNAdata, 1, var)
mir_var <- apply(miRdata, 1, var)
genes <- names(sort(gene_var, decreasing = TRUE))[1:G]
mirs <- names(sort(mir_var, decreasing = TRUE))[1:M]

missing_ex <- c(setdiff(exemplar_pairs$gene, rownames(mRNAdata)), setdiff(exemplar_pairs$mir, rownames(miRdata)))
if (length(missing_ex) > 0) {
  message("Exemplar features not in the data (dropped): ", paste(missing_ex, collapse = ", "))
  exemplar_pairs <- exemplar_pairs[exemplar_pairs$gene %in% rownames(mRNAdata) & exemplar_pairs$mir %in% rownames(miRdata), ]
}
extra_genes <- setdiff(unique(exemplar_pairs$gene), genes)
extra_mirs <- setdiff(unique(exemplar_pairs$mir), mirs)

# Feature matrix: domain genes, domain miRNAs, then exemplar features outside the domain.
E_all <- rbind(mRNAdata[genes, , drop = FALSE], miRdata[mirs, , drop = FALSE],
               mRNAdata[extra_genes, , drop = FALSE], miRdata[extra_mirs, , drop = FALSE])
row_gene <- seq_len(G)
row_mir <- G + seq_len(M)
feature_type <- c(rep("mRNA", G), rep("miRNA", M), rep("mRNA", length(extra_genes)), rep("miRNA", length(extra_mirs)))
toc("load_data", t0)

# ---------------------------------------------------------------------------
# Covariates
# ---------------------------------------------------------------------------
t0 <- tic()
cov <- data.frame(sample = sample_id, key15 = key15, cohort = cohort_all, sample_type = sample_type,
                  stringsAsFactors = FALSE)

# Leukocyte fraction: aliquot barcodes, matched on the first 15 characters (patient + sample type);
# several aliquots of one sample are averaged.
lf <- utils::read.delim(ext_files$leuk$path, header = FALSE, col.names = c("cancer", "aliquot", "LF"),
                        stringsAsFactors = FALSE)
lf$key15 <- substr(lf$aliquot, 1, 15)
lf <- lf[is.finite(lf$LF), ]
lf_agg <- stats::aggregate(LF ~ key15, data = lf, FUN = mean)
cov$LF <- lf_agg$LF[match(cov$key15, lf_agg$key15)]

# CIBERSORT relative fractions (22 subsets), matched on the first 15 characters, aliquots averaged.
cs <- utils::read.delim(ext_files$cibersort$path, stringsAsFactors = FALSE, check.names = FALSE)
cs_cols <- setdiff(colnames(cs), c("SampleID", "CancerType", "P.value", "Correlation", "RMSE"))
stopifnot(length(cs_cols) == 22)
cs$key15 <- gsub(".", "-", substr(cs$SampleID, 1, 15), fixed = TRUE)
cs_agg <- stats::aggregate(cs[, c(cs_cols, "P.value")], by = list(key15 = cs$key15), FUN = mean)
m_cs <- match(cov$key15, cs_agg$key15)
cs_rel <- as.matrix(cs_agg[m_cs, cs_cols])
rownames(cs_rel) <- cov$sample
colnames(cs_rel) <- paste0("CS_", gsub("[^A-Za-z0-9]+", "_", cs_cols))
colnames(cs_rel) <- gsub("_+$", "", colnames(cs_rel))
cov$CS_pvalue <- cs_agg$P.value[m_cs]

# ABSOLUTE purity, matched on the 15-character `array` column.
ab <- utils::read.delim(ext_files$absolute$path, stringsAsFactors = FALSE, check.names = FALSE)
ab <- ab[is.finite(ab$purity), ]
cov$purity <- ab$purity[match(cov$key15, ab$array)]

# Marker-based scores (z-score across all samples, then average).
zrow <- function(x) (x - mean(x)) / stats::sd(x)
Zimm <- t(apply(mRNAdata[immune_markers, , drop = FALSE], 1, zrow))
Zstr <- t(apply(mRNAdata[stromal_markers, , drop = FALSE], 1, zrow))
cov$marker_immune <- colMeans(Zimm)
cov$marker_stromal <- colMeans(Zstr)
cov$marker_purity_proxy <- -(cov$marker_immune + cov$marker_stromal)

# ESTIMATE stromal / immune scores. In-memory re-implementation of
# estimate::estimateScore (v1.0.13) restricted to the package's common genes,
# verified below against the package's own function on a subset of samples.
if (!requireNamespace("estimate", quietly = TRUE)) {
  utils::install.packages("estimate", repos = "https://download.r-forge.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages(library(estimate))
utils::data("common_genes", package = "estimate", envir = environment())
utils::data("SI_geneset", package = "estimate", envir = environment())
estimate_version <- as.character(utils::packageVersion("estimate"))
est_genes <- intersect(common_genes$GeneSymbol, rownames(mRNAdata))
n_common_matched <- length(est_genes)
estimate_scores_inmem <- function(m) {
  # m: genes x samples (expression on any monotone scale). Mirrors estimateScore().
  Ng <- nrow(m)
  Ns <- ncol(m)
  for (j in seq_len(Ns)) m[, j] <- rank(m[, j], ties.method = "average")
  m <- 10000 * m / Ng
  gs <- as.matrix(SI_geneset[, -1])
  gene.names <- rownames(m)
  out <- matrix(NA_real_, nrow = 2, ncol = Ns, dimnames = list(c("StromalScore", "ImmuneScore"), colnames(m)))
  for (gs.i in 1:2) {
    gene.overlap <- intersect(gs[gs.i, ], gene.names)
    gene.set2 <- match(gene.overlap, gene.names)
    Nh <- length(gene.set2)
    Nm <- Ng - Nh
    for (S in seq_len(Ns)) {
      gene.list <- order(m[, S], decreasing = TRUE)
      correl.vector <- abs(m[gene.list, S])^0.25
      TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
      no.TAG <- 1 - TAG
      sum.correl <- sum(correl.vector[TAG == 1])
      F0 <- cumsum(no.TAG / Nm)
      Fn <- cumsum(TAG * correl.vector / sum.correl)
      out[gs.i, S] <- sum(Fn - F0)
    }
  }
  out
}
est <- estimate_scores_inmem(mRNAdata[est_genes, , drop = FALSE])
cov$est_stromal <- est["StromalScore", ]
cov$est_immune <- est["ImmuneScore", ]
cov$est_score <- cov$est_stromal + cov$est_immune

# Verification against the package on the first 60 samples (GCT round trip in a temp dir).
verify_idx <- seq_len(min(60, N_all))
tmpd <- tempfile("estimate_verify_")
dir.create(tmpd)
gct_in <- file.path(tmpd, "input.gct")
gct_out <- file.path(tmpd, "scores.gct")
df_in <- as.data.frame(mRNAdata[est_genes, verify_idx, drop = FALSE])
estimate::outputGCT(df_in, gct_in)
invisible(utils::capture.output(estimate::estimateScore(gct_in, gct_out, platform = "illumina")))
pkg_scores <- utils::read.delim(gct_out, skip = 2, row.names = 1, check.names = FALSE)[, -1]
estimate_verify_maxdiff <- max(abs(as.numeric(as.matrix(pkg_scores[c("StromalScore", "ImmuneScore"), ])) -
                                     as.numeric(est[, verify_idx])))
message(sprintf("ESTIMATE in-memory vs package, max abs difference over %d samples: %.3e", length(verify_idx), estimate_verify_maxdiff))
if (!is.finite(estimate_verify_maxdiff) || estimate_verify_maxdiff > 1e-6) stop("ESTIMATE re-implementation does not match the package.")
unlink(tmpd, recursive = TRUE)

# CIBERSORT absolute fractions = relative x LF (Thorsson et al.).
cs_abs <- cs_rel * cov$LF
colnames(cs_abs) <- sub("^CS_", "CSabs_", colnames(cs_rel))

cov_out <- cbind(cov, as.data.frame(cs_rel), as.data.frame(cs_abs))
utils::write.csv(cov_out, file.path(out_dir, "covariates_per_sample.csv"), row.names = FALSE)

# Summary of covariates: source, match rate, distribution.
is_tumor <- cov$sample_type != "normal"
summ_row <- function(var, x, source, file = NA, url = NA, md5 = NA, note = NA) {
  q <- stats::quantile(x, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
  data.frame(covariate = var, source = source, file = file, url = url, md5 = md5,
             n_available = sum(is.finite(x)), n_samples = N_all, match_rate = sum(is.finite(x)) / N_all,
             n_available_tumor_or_metastasis = sum(is.finite(x) & is_tumor), n_tumor_or_metastasis = sum(is_tumor),
             match_rate_tumor_or_metastasis = sum(is.finite(x) & is_tumor) / sum(is_tumor),
             n_available_normal = sum(is.finite(x) & !is_tumor), n_normal = sum(!is_tumor),
             mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE),
             min = q[1], q25 = q[2], median = q[3], q75 = q[4], max = q[5], note = note,
             stringsAsFactors = FALSE)
}
cov_summary <- rbind(
  summ_row("LF", cov$LF, ext_files$leuk$source, ext_files$leuk$name, ext_files$leuk$url, ext_files$leuk$md5,
           sprintf("file has %d aliquot rows (%d unique 15-char sample keys, %d with finite LF); matched on the first 15 barcode characters; aliquots averaged",
                   nrow(lf), length(unique(lf$key15)), nrow(lf_agg))),
  summ_row("purity", cov$purity, ext_files$absolute$source, ext_files$absolute$name, ext_files$absolute$url, ext_files$absolute$md5,
           sprintf("file has %d rows with finite purity; matched on the 15-char `array` column", nrow(ab))),
  summ_row("CS_pvalue", cov$CS_pvalue, ext_files$cibersort$source, ext_files$cibersort$name, ext_files$cibersort$url, ext_files$cibersort$md5,
           sprintf("CIBERSORT deconvolution p-value; file has %d rows (%d unique 15-char keys); the %d relative fractions sum to 1 per sample; matched on the first 15 characters; aliquots averaged",
                   nrow(cs), nrow(cs_agg), length(cs_cols))),
  summ_row("marker_immune", cov$marker_immune, "computed here from mRNAdata (marker-based proxy)", note = paste("mean z-score of", paste(immune_markers, collapse = ","))),
  summ_row("marker_stromal", cov$marker_stromal, "computed here from mRNAdata (marker-based proxy)", note = paste("mean z-score of", paste(stromal_markers, collapse = ","))),
  summ_row("marker_purity_proxy", cov$marker_purity_proxy, "computed here from mRNAdata (marker-based proxy)", note = "-(marker_immune + marker_stromal)"),
  summ_row("est_stromal", cov$est_stromal, sprintf("ESTIMATE (Yoshihara 2013) via estimate v%s algorithm, platform illumina, computed here from mRNAdata", estimate_version),
           note = sprintf("%d of %d ESTIMATE common genes present; in-memory implementation verified against the package (max abs diff %.2e on %d samples)", n_common_matched, nrow(common_genes), estimate_verify_maxdiff, length(verify_idx))),
  summ_row("est_immune", cov$est_immune, sprintf("ESTIMATE via estimate v%s algorithm, computed here from mRNAdata", estimate_version)),
  summ_row("est_score", cov$est_score, sprintf("ESTIMATE score = stromal + immune, estimate v%s algorithm", estimate_version))
)
for (j in seq_len(ncol(cs_abs))) {
  cov_summary <- rbind(cov_summary, summ_row(colnames(cs_abs)[j], cs_abs[, j], "CIBERSORT relative fraction x LF (Thorsson et al. 2018)",
                                             ext_files$cibersort$name, ext_files$cibersort$url, ext_files$cibersort$md5))
}
utils::write.csv(cov_summary, file.path(out_dir, "covariates_summary.csv"), row.names = FALSE)

# Agreement between sources (complete pairs).
cor_vars <- c("LF", "purity", "marker_immune", "marker_stromal", "marker_purity_proxy", "est_immune", "est_stromal", "est_score")
cc_rows <- list()
for (a in seq_along(cor_vars)) for (b in seq_along(cor_vars)) if (a < b) {
  x <- cov[[cor_vars[a]]]; y <- cov[[cor_vars[b]]]
  ok <- is.finite(x) & is.finite(y)
  cc_rows[[length(cc_rows) + 1]] <- data.frame(var1 = cor_vars[a], var2 = cor_vars[b], n = sum(ok),
                                               pearson = stats::cor(x[ok], y[ok]),
                                               spearman = stats::cor(x[ok], y[ok], method = "spearman"),
                                               stringsAsFactors = FALSE)
}
cov_cor <- do.call(rbind, cc_rows)
utils::write.csv(cov_cor, file.path(out_dir, "covariate_correlations.csv"), row.names = FALSE)
toc("covariates", t0)

# ---------------------------------------------------------------------------
# Blocks, covariate sets, conditions
# ---------------------------------------------------------------------------
matched <- is.finite(cov$LF) & is.finite(cov$purity) & apply(is.finite(cs_abs), 1, all)
blocks <- list(
  all = list(samples = seq_len(N_all), label = "all 8,890 samples"),
  matched = list(samples = which(matched), label = "samples with LF, CIBERSORT and ABSOLUTE purity")
)
# Drop cohorts below min_n_per_cohort inside each block.
for (b in names(blocks)) {
  s <- blocks[[b]]$samples
  tab <- table(cohort_all[s])
  keep_ct <- names(tab)[tab >= min_n_per_cohort]
  blocks[[b]]$samples <- s[cohort_all[s] %in% keep_ct]
  blocks[[b]]$cohorts <- keep_ct
  blocks[[b]]$dropped_cohorts <- setdiff(levels(cancerType), keep_ct)
}

covsets <- list(
  none = list(cols = character(0), label = "unadjusted"),
  leuk = list(cols = "LF", label = "leukocyte fraction (Thorsson)"),
  leuk_purity = list(cols = c("LF", "purity"), label = "leukocyte fraction + ABSOLUTE purity"),
  cibersort = list(cols = colnames(cs_abs), label = "22 CIBERSORT absolute fractions"),
  cibersort_purity = list(cols = c(colnames(cs_abs), "purity"), label = "22 CIBERSORT absolute fractions + ABSOLUTE purity"),
  marker = list(cols = c("marker_immune", "marker_stromal"), label = "marker immune + stromal scores (proxy)"),
  estimate = list(cols = c("est_immune", "est_stromal"), label = "ESTIMATE immune + stromal scores")
)
block_sets <- list(all = c("none", "marker", "estimate"),
                   matched = c("none", "leuk", "leuk_purity", "cibersort", "cibersort_purity", "marker", "estimate"))
modes <- c("pooled", "within", "within_between", "within_resid")

conditions <- list()
for (b in names(block_sets)) for (cs_name in block_sets[[b]]) {
  if (cs_name == "none") {
    conditions[[length(conditions) + 1]] <- list(block = b, covset = "none", mode = "none")
  } else {
    for (md in modes) conditions[[length(conditions) + 1]] <- list(block = b, covset = cs_name, mode = md)
  }
}
cond_id <- vapply(conditions, function(cn) paste(cn$block, cn$covset, cn$mode, sep = "|"), character(1))
names(conditions) <- cond_id
message(length(conditions), " conditions")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
cor_rows <- function(A, Bm, idx) {
  A <- A[, idx, drop = FALSE]
  Bm <- Bm[, idx, drop = FALSE]
  scaled_B <- scale(t(Bm)) / (ncol(Bm) - 1)
  std_A <- t(scale(t(A)))
  std_A %*% scaled_B
}

# Orthonormal basis of the column space of X (rank-aware).
qbasis <- function(X) {
  q <- qr(X)
  qr.Q(q)[, seq_len(q$rank), drop = FALSE]
}

# Regress each row of E on the covariates W (N x p) and return the adjusted matrix.
adjust_matrix <- function(E, ctx, W, mode) {
  p <- ncol(W)
  N <- ncol(E)
  if (mode == "pooled") {
    Q <- qbasis(cbind(1, W))
    A <- E - (E %*% Q) %*% t(Q)
    R2_pooled <- 1 - rowSums(A^2) / rowSums((E - rowMeans(E))^2)
    A <- A + rowMeans(E)
    return(list(A = A, p_within = 0L, R2 = R2_pooled))
  }
  lev <- levels(ctx)
  k <- length(lev)
  A <- E
  Mc <- matrix(NA_real_, nrow = nrow(E), ncol = k, dimnames = list(rownames(E), lev))
  ss_res <- ss_within <- 0
  for (j in seq_len(k)) {
    idx <- which(ctx == lev[j])
    Ec <- E[, idx, drop = FALSE]
    Mc[, j] <- rowMeans(Ec)
    Q <- qbasis(cbind(1, W[idx, , drop = FALSE]))
    Rc <- Ec - (Ec %*% Q) %*% t(Q)
    A[, idx] <- Rc
    ss_res <- ss_res + rowSums(Rc^2)
    ss_within <- ss_within + rowSums((Ec - Mc[, j])^2)
  }
  R2_within <- 1 - ss_res / ss_within
  if (mode == "within_resid") return(list(A = A, p_within = p, R2 = R2_within))
  n_c <- as.numeric(table(ctx)[lev])
  if (mode == "within") {
    for (j in seq_len(k)) {
      idx <- which(ctx == lev[j])
      A[, idx] <- A[, idx, drop = FALSE] + Mc[, j]
    }
    return(list(A = A, p_within = p, R2 = R2_within))
  }
  # within_between: cohort-level WLS (weights N_c) of cohort means on cohort-mean covariates.
  Wbar <- t(vapply(seq_len(k), function(j) colMeans(W[ctx == lev[j], , drop = FALSE]), numeric(p)))
  sw <- sqrt(n_c)
  Xb <- cbind(1, matrix(Wbar, nrow = k)) * sw
  Qb <- qbasis(Xb)
  Mt <- Mc * rep(sw, each = nrow(Mc))
  Rb <- (Mt - (Mt %*% Qb) %*% t(Qb)) / rep(sw, each = nrow(Mc))
  gm <- as.numeric(Mc %*% n_c) / sum(n_c)
  for (j in seq_len(k)) {
    idx <- which(ctx == lev[j])
    A[, idx] <- A[, idx, drop = FALSE] + Rb[, j] + gm
  }
  list(A = A, p_within = p, R2 = R2_within)
}

# Per-pair statistics from the pooled correlations r_global (P) and the
# within-cohort matrix R (P x k), with per-cohort effective sizes n_eff (k).
pair_metrics <- function(r_global, R, n_eff, eps_grid) {
  P <- nrow(R)
  k <- ncol(R)
  w_c <- n_eff - 3
  valid <- is.finite(R) & abs(R) < 1
  valid <- valid & matrix(w_c > 0, nrow = P, ncol = k, byrow = TRUE)
  k_valid <- rowSums(valid)
  max_abs <- rep(0, P)
  for (j in seq_len(k)) {
    a <- abs(R[, j])
    a[!valid[, j]] <- 0
    max_abs <- pmax(max_abs, a)
  }
  W <- matrix(w_c, nrow = P, ncol = k, byrow = TRUE)
  W[!valid] <- 0
  Z <- contextcorr_fisher_z(R)
  Z[!valid] <- 0
  sum_w <- rowSums(W)
  sum_wz <- rowSums(W * Z)
  sum_wz2 <- rowSums(W * Z^2)
  Q <- sum_wz2 - sum_wz^2 / pmax(sum_w, 1e-12)
  df <- pmax(k_valid - 1L, 0L)
  ok <- df >= 1L & is.finite(Q) & Q >= 0
  I2 <- rep(NA_real_, P)
  I2[ok] <- pmax(0, (Q[ok] - df[ok]) / pmax(Q[ok], 1e-12))
  out <- list(r_global = as.numeric(r_global), k_valid = k_valid, max_abs = max_abs, I2 = I2)
  for (eps in eps_grid) {
    tag <- gsub("\\.", "p", formatC(eps, format = "f", digits = 2))
    out[[paste0("n_pos_eps", tag)]] <- rowSums(valid & (R >= eps))
    out[[paste0("n_neg_eps", tag)]] <- rowSums(valid & (R <= -eps))
  }
  out
}

eps_tag <- function(e) gsub("\\.", "p", formatC(e, format = "f", digits = 2))

# Point-estimate summary of one condition.
rate_summary <- function(pm, eps = 0.05, t = 0.2, min_k = 10) {
  valid <- pm$k_valid >= min_k & is.finite(pm$r_global)
  rg <- pm$r_global
  sg <- sign(rg)
  rr <- pm$r_resid
  one_eps <- function(e) {
    pos <- pm[[paste0("n_pos_eps", eps_tag(e))]]
    neg <- pm[[paste0("n_neg_eps", eps_tag(e))]]
    mixed <- pos > 0 & neg > 0
    maj <- sign(pos - neg)
    elig_t <- valid & abs(rg) >= t & maj != 0 & sg != 0
    elig_0 <- valid & maj != 0 & sg != 0
    list(frac_mixed = mean(mixed[valid]),
         frac_mixed_given_maxabs = if (any(valid & pm$max_abs >= t)) mean(mixed[valid & pm$max_abs >= t]) else NA_real_,
         n_elig_simpson = sum(elig_t),
         frac_simpson = if (sum(elig_t) > 0) mean(sg[elig_t] != maj[elig_t]) else NA_real_,
         n_elig_simpson_t0 = sum(elig_0),
         frac_simpson_t0 = if (sum(elig_0) > 0) mean(sg[elig_0] != maj[elig_0]) else NA_real_,
         frac_tie = mean((maj == 0)[valid]))
  }
  m <- one_eps(eps)
  m0 <- one_eps(0)
  m10 <- one_eps(0.10)
  flip_den <- valid & sg != 0 & is.finite(rr) & sign(rr) != 0
  flip_den_t <- flip_den & abs(rg) >= t
  data.frame(
    n_pairs = length(rg), n_valid = sum(valid),
    frac_mixed_eps0p05 = m$frac_mixed,
    frac_maxabs_ge0p2 = mean((pm$max_abs >= t)[valid]),
    frac_mixed_eps0p05_given_maxabs_ge0p2 = m$frac_mixed_given_maxabs,
    frac_rglobal_ge0p2 = mean((abs(rg) >= t)[valid]),
    median_abs_rglobal = stats::median(abs(rg[valid])),
    n_elig_simpson_t0p2 = m$n_elig_simpson,
    frac_simpson_eps0p05_t0p2 = m$frac_simpson,
    n_elig_simpson_t0 = m$n_elig_simpson_t0,
    frac_simpson_eps0p05_t0 = m$frac_simpson_t0,
    frac_majority_tie_eps0p05 = m$frac_tie,
    n_flip_den = sum(flip_den),
    frac_resid_flip = mean((sg != sign(rr))[flip_den]),
    n_flip_den_t0p2 = sum(flip_den_t),
    frac_resid_flip_t0p2 = if (sum(flip_den_t) > 0) mean((sg != sign(rr))[flip_den_t]) else NA_real_,
    median_I2 = stats::median(pm$I2[valid], na.rm = TRUE),
    mean_I2 = mean(pm$I2[valid], na.rm = TRUE),
    I2_q25 = as.numeric(stats::quantile(pm$I2[valid], 0.25, na.rm = TRUE)),
    I2_q75 = as.numeric(stats::quantile(pm$I2[valid], 0.75, na.rm = TRUE)),
    frac_I2_gt0p75 = mean(pm$I2[valid] > 0.75, na.rm = TRUE),
    frac_mixed_eps0 = m0$frac_mixed,
    frac_simpson_eps0_t0p2 = m0$frac_simpson,
    frac_mixed_eps0p10 = m10$frac_mixed,
    frac_simpson_eps0p10_t0p2 = m10$frac_simpson,
    stringsAsFactors = FALSE)
}

# Between-cohort variance of cohort means (N_c-weighted and unweighted) and eta^2, per feature.
between_stats <- function(X, ctx) {
  lev <- levels(ctx)
  n_c <- as.numeric(table(ctx)[lev])
  Mc <- vapply(lev, function(l) rowMeans(X[, ctx == l, drop = FALSE]), numeric(nrow(X)))
  gm_w <- as.numeric(Mc %*% n_c) / sum(n_c)
  var_between_w <- as.numeric(((Mc - gm_w)^2) %*% n_c) / sum(n_c)
  var_between_u <- apply(Mc, 1, stats::var)
  ss_total <- rowSums((X - rowMeans(X))^2)
  eta2 <- var_between_w * sum(n_c) / ss_total
  list(var_between_w = var_between_w, var_between_u = var_between_u, eta2 = eta2, Mc = Mc)
}

# Exemplar pair statistics on a (possibly adjusted) matrix.
exemplar_stats <- function(A, ctx, n_eff_by_ctx, gene, mir, eps = 0.05, t = 0.2) {
  x <- A[gene, ]
  y <- A[mir, ]
  cc <- context_correlations(x, y, ctx)
  r_ctx <- cc$r_by_context
  n_eff <- n_eff_by_ctx[names(r_ctx)]
  het <- heterogeneity_Q_I2(r_ctx, n_eff)
  sf <- simpson_flag(cc$r_global, r_ctx, eps = eps)
  rres <- mean_residualized_cor(x, y, ctx)
  cl <- classify_pair(cc$r_global, r_ctx, n_eff, eps = eps)
  data.frame(gene = gene, mir = mir, r_global = cc$r_global, r_resid = rres,
             k_cohorts = sum(is.finite(r_ctx)), n_neg_any = sum(r_ctx < 0, na.rm = TRUE),
             n_pos_eps = sf$n_pos, n_neg_eps = sf$n_neg, majority_sign = sf$majority_sign,
             simpson_reversal = sf$flag, eligible_abs_rglobal_ge0p2 = abs(cc$r_global) >= t,
             resid_flip = sign(cc$r_global) != sign(rres),
             Q = het$Q, I2 = het$I2, classification = cl$label,
             r_c_min = min(r_ctx, na.rm = TRUE), r_c_max = max(r_ctx, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Run all conditions (forked workers; reference BLAS is single-threaded)
# ---------------------------------------------------------------------------
t0 <- tic()
covmat_all <- cbind(as.matrix(cov[, c("LF", "purity", "marker_immune", "marker_stromal", "est_immune", "est_stromal")]), cs_abs)

run_condition <- function(cn) {
  bl <- blocks[[cn$block]]
  s <- bl$samples
  ctx <- droplevels(factor(cohort_all[s], levels = levels(cancerType)))
  lev <- levels(ctx)
  n_c <- as.numeric(table(ctx)[lev])
  names(n_c) <- lev
  E <- E_all[, s, drop = FALSE]
  if (cn$covset == "none") {
    A <- E
    p_within <- 0L
    R2 <- rep(NA_real_, nrow(E))
  } else {
    W <- covmat_all[s, covsets[[cn$covset]]$cols, drop = FALSE]
    stopifnot(all(is.finite(W)))
    adj <- adjust_matrix(E, ctx, W, cn$mode)
    A <- adj$A
    p_within <- adj$p_within
    R2 <- adj$R2
  }
  n_eff <- n_c - p_within
  names(n_eff) <- lev
  Ag <- A[row_gene, , drop = FALSE]
  Am <- A[row_mir, , drop = FALSE]
  P <- G * M
  R <- matrix(NA_real_, nrow = P, ncol = length(lev))
  for (j in seq_along(lev)) {
    idx <- which(ctx == lev[j])
    R[, j] <- as.numeric(cor_rows(Ag, Am, idx))
  }
  r_global <- as.numeric(cor_rows(Ag, Am, seq_along(s)))
  Ares <- residualize_matrix(A, ctx)
  r_resid <- as.numeric(cor_rows(Ares[row_gene, , drop = FALSE], Ares[row_mir, , drop = FALSE], seq_along(s)))
  pm <- pair_metrics(r_global, R, n_eff, eps_grid)
  pm$r_resid <- r_resid
  rm(R)
  summ <- rate_summary(pm, eps = eps_main, t = t_main, min_k = min_k)
  summ$n_samples <- length(s)
  summ$n_cohorts <- length(lev)
  summ$min_cohort_n <- min(n_c)
  summ$p_covariates <- if (cn$covset == "none") 0L else length(covsets[[cn$covset]]$cols)
  summ$p_within <- p_within
  summ$n_cohort_level_residual_df <- if (identical(cn$mode, "within_between")) length(lev) - 1L - summ$p_covariates else NA_integer_
  summ$min_n_eff <- min(n_eff)
  summ$median_R2_covariates <- stats::median(R2[c(row_gene, row_mir)], na.rm = TRUE)
  summ$median_R2_covariates_mRNA <- stats::median(R2[row_gene], na.rm = TRUE)
  summ$median_R2_covariates_miRNA <- stats::median(R2[row_mir], na.rm = TRUE)
  # Between-cohort mean-shift statistics before/after (domain features only).
  dom <- c(row_gene, row_mir)
  bs_before <- between_stats(E[dom, , drop = FALSE], ctx)
  bs_after <- between_stats(A[dom, , drop = FALSE], ctx)
  ms <- list(ratio_w = bs_after$var_between_w / bs_before$var_between_w,
             ratio_u = bs_after$var_between_u / bs_before$var_between_u,
             eta2_before = bs_before$eta2, eta2_after = bs_after$eta2,
             feature_type = feature_type[dom])
  # Exemplars.
  ex <- do.call(rbind, lapply(seq_len(nrow(exemplar_pairs)), function(i) {
    cbind(set = exemplar_pairs$set[i],
          exemplar_stats(A, ctx, n_eff, exemplar_pairs$gene[i], exemplar_pairs$mir[i], eps = eps_main, t = t_main),
          stringsAsFactors = FALSE)
  }))
  list(summary = summ, pm = pm, meanshift = ms, exemplars = ex, n_c = n_c, n_eff = n_eff)
}

res <- parallel::mclapply(conditions, run_condition, mc.cores = min(n_cores, length(conditions)), mc.preschedule = FALSE)
failed <- vapply(res, function(r) inherits(r, "try-error") || is.null(r$summary), logical(1))
if (any(failed)) stop("Conditions failed: ", paste(cond_id[failed], collapse = ", "), "\n", res[[which(failed)[1]]])
toc("conditions", t0)

# ---------------------------------------------------------------------------
# Pair-level bootstrap (B replicates, same resampled pairs for every condition)
# ---------------------------------------------------------------------------
t0 <- tic()
P <- G * M
boot_metrics <- c("frac_mixed_eps0p05", "frac_mixed_eps0p05_given_maxabs_ge0p2", "frac_rglobal_ge0p2",
                  "frac_simpson_eps0p05_t0p2", "frac_simpson_eps0p05_t0", "frac_resid_flip", "frac_resid_flip_t0p2",
                  "median_I2", "frac_I2_gt0p75", "frac_simpson_eps0_t0p2", "frac_simpson_eps0p10_t0p2")

# Precompute index sets per condition so a replicate only needs weighted sums.
prep_boot <- function(pm) {
  valid <- pm$k_valid >= min_k & is.finite(pm$r_global)
  rg <- pm$r_global
  sg <- sign(rg)
  rr <- pm$r_resid
  one <- function(e) {
    pos <- pm[[paste0("n_pos_eps", eps_tag(e))]]
    neg <- pm[[paste0("n_neg_eps", eps_tag(e))]]
    maj <- sign(pos - neg)
    elig <- valid & abs(rg) >= t_main & maj != 0 & sg != 0
    list(elig = which(elig), rev = which(elig & sg != maj), mixed = which(valid & pos > 0 & neg > 0),
         elig0 = which(valid & maj != 0 & sg != 0), rev0 = which(valid & maj != 0 & sg != 0 & sg != maj))
  }
  e05 <- one(0.05); e0 <- one(0); e10 <- one(0.10)
  flip_den <- valid & sg != 0 & is.finite(rr) & sign(rr) != 0
  vI <- which(valid & is.finite(pm$I2))
  o <- order(pm$I2[vI])
  list(valid = which(valid), mixed = e05$mixed, maxabs = which(valid & pm$max_abs >= t_main),
       mixed_maxabs = which(valid & pm$max_abs >= t_main & pm[["n_pos_eps0p05"]] > 0 & pm[["n_neg_eps0p05"]] > 0),
       rg_t = which(valid & abs(rg) >= t_main),
       elig = e05$elig, rev = e05$rev, elig0 = e05$elig0, rev0 = e05$rev0,
       flip_den = which(flip_den), flip = which(flip_den & sg != sign(rr)),
       flip_den_t = which(flip_den & abs(rg) >= t_main), flip_t = which(flip_den & abs(rg) >= t_main & sg != sign(rr)),
       I2_idx = vI[o], I2_sorted = pm$I2[vI][o], I2_gt = which(valid & is.finite(pm$I2) & pm$I2 > 0.75),
       elig_e0 = e0$elig, rev_e0 = e0$rev, elig_e10 = e10$elig, rev_e10 = e10$rev)
}
bprep <- lapply(res, function(r) prep_boot(r$pm))

boot_one <- function(b) {
  set.seed(seed_boot + b)
  w <- tabulate(sample.int(P, P, replace = TRUE), nbins = P)
  out <- matrix(NA_real_, nrow = length(bprep), ncol = length(boot_metrics), dimnames = list(names(bprep), boot_metrics))
  for (i in seq_along(bprep)) {
    bp <- bprep[[i]]
    sw <- function(ix) sum(w[ix])
    nv <- sw(bp$valid)
    ww <- w[bp$I2_idx]
    cs <- cumsum(ww)
    med <- if (length(cs) > 0 && cs[length(cs)] > 0) bp$I2_sorted[which(cs >= cs[length(cs)] / 2)[1]] else NA_real_
    out[i, ] <- c(sw(bp$mixed) / nv, sw(bp$mixed_maxabs) / sw(bp$maxabs), sw(bp$rg_t) / nv,
                  sw(bp$rev) / sw(bp$elig), sw(bp$rev0) / sw(bp$elig0),
                  sw(bp$flip) / sw(bp$flip_den), sw(bp$flip_t) / sw(bp$flip_den_t),
                  med, sw(bp$I2_gt) / sw(bp$I2_idx),
                  sw(bp$rev_e0) / sw(bp$elig_e0), sw(bp$rev_e10) / sw(bp$elig_e10))
  }
  out
}
boot_list <- parallel::mclapply(seq_len(B), boot_one, mc.cores = n_cores)
boot_arr <- simplify2array(boot_list)       # conditions x metrics x B
toc("bootstrap", t0)

# ---------------------------------------------------------------------------
# adjusted_rates.csv
# ---------------------------------------------------------------------------
t0 <- tic()
baseline_id <- function(b) paste(b, "none", "none", sep = "|")
rows <- list()
for (i in seq_along(conditions)) {
  cn <- conditions[[i]]
  sm <- res[[i]]$summary
  base <- res[[baseline_id(cn$block)]]$summary
  ib <- which(cond_id == baseline_id(cn$block))
  row <- data.frame(condition = cond_id[i], block = cn$block, block_label = blocks[[cn$block]]$label,
                    covset = cn$covset, covset_label = covsets[[cn$covset]]$label, mode = cn$mode,
                    stringsAsFactors = FALSE)
  row <- cbind(row, sm)
  for (mname in boot_metrics) {
    v <- boot_arr[i, mname, ]
    row[[paste0(mname, "_ci_lo")]] <- as.numeric(stats::quantile(v, 0.025, na.rm = TRUE))
    row[[paste0(mname, "_ci_hi")]] <- as.numeric(stats::quantile(v, 0.975, na.rm = TRUE))
  }
  # Unadjusted values on the same block and the same pairs, with paired-difference CIs.
  for (mname in c("frac_simpson_eps0p05_t0p2", "frac_mixed_eps0p05", "frac_mixed_eps0p05_given_maxabs_ge0p2", "median_I2",
                  "frac_resid_flip", "frac_resid_flip_t0p2", "frac_rglobal_ge0p2")) {
    row[[paste0("unadj_", mname)]] <- base[[mname]]
    d <- boot_arr[i, mname, ] - boot_arr[ib, mname, ]
    row[[paste0("diff_", mname)]] <- sm[[mname]] - base[[mname]]
    row[[paste0("diff_", mname, "_ci_lo")]] <- as.numeric(stats::quantile(d, 0.025, na.rm = TRUE))
    row[[paste0("diff_", mname, "_ci_hi")]] <- as.numeric(stats::quantile(d, 0.975, na.rm = TRUE))
  }
  rows[[i]] <- row
}
rates_df <- do.call(rbind, rows)
rates_df$B <- B
rates_df$eps <- eps_main
rates_df$t_rglobal <- t_main
rates_df$min_k <- min_k
rates_df$G <- G
rates_df$M <- M
utils::write.csv(rates_df, file.path(out_dir, "adjusted_rates.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# reversal_transitions.csv
# ---------------------------------------------------------------------------
rev_sets <- function(pm) {
  valid <- pm$k_valid >= min_k & is.finite(pm$r_global)
  rg <- pm$r_global
  sg <- sign(rg)
  pos <- pm[["n_pos_eps0p05"]]
  neg <- pm[["n_neg_eps0p05"]]
  maj <- sign(pos - neg)
  elig <- valid & abs(rg) >= t_main & maj != 0 & sg != 0
  list(elig = elig, rev = elig & sg != maj, rg = rg, maj = maj, sg = sg, valid = valid)
}
trans_rows <- list()
for (i in seq_along(conditions)) {
  cn <- conditions[[i]]
  if (cn$covset == "none") next
  u <- rev_sets(res[[baseline_id(cn$block)]]$pm)
  a <- rev_sets(res[[i]]$pm)
  n_u <- sum(u$rev)
  n_a <- sum(a$rev)
  persist <- u$rev & a$rev
  vanish <- u$rev & !a$rev
  vanish_inelig <- vanish & !a$elig
  vanish_inelig_small <- vanish & !a$elig & abs(a$rg) < t_main
  vanish_elig <- vanish & a$elig
  newr <- a$rev & !u$rev
  new_from_inelig <- newr & !u$elig
  new_from_elig <- newr & u$elig
  # Sign changes of r_global among unadjusted reversals.
  sign_change_rg <- u$rev & (a$sg != u$sg)
  maj_change <- u$rev & (a$maj != u$maj)
  trans_rows[[length(trans_rows) + 1]] <- data.frame(
    condition = cond_id[i], block = cn$block, covset = cn$covset, mode = cn$mode,
    n_unadjusted_reversals = n_u, n_adjusted_reversals = n_a,
    n_persist = sum(persist), frac_persist = sum(persist) / n_u,
    n_vanish = sum(vanish), frac_vanish = sum(vanish) / n_u,
    n_vanish_no_longer_eligible = sum(vanish_inelig), frac_vanish_no_longer_eligible = sum(vanish_inelig) / n_u,
    n_vanish_abs_rglobal_below_t = sum(vanish_inelig_small), frac_vanish_abs_rglobal_below_t = sum(vanish_inelig_small) / n_u,
    n_vanish_eligible_not_reversed = sum(vanish_elig), frac_vanish_eligible_not_reversed = sum(vanish_elig) / n_u,
    n_new = sum(newr), frac_new_relative_to_unadjusted = sum(newr) / n_u, frac_new_among_adjusted = sum(newr) / n_a,
    n_new_from_ineligible = sum(new_from_inelig), n_new_from_eligible_nonreversed = sum(new_from_elig),
    n_unadj_rev_rglobal_sign_changed = sum(sign_change_rg), frac_unadj_rev_rglobal_sign_changed = sum(sign_change_rg) / n_u,
    n_unadj_rev_majority_changed = sum(maj_change), frac_unadj_rev_majority_changed = sum(maj_change) / n_u,
    n_eligible_unadjusted = sum(u$elig), n_eligible_adjusted = sum(a$elig),
    n_eligible_both = sum(u$elig & a$elig),
    frac_reversal_among_eligible_both_unadjusted = sum(u$rev & u$elig & a$elig) / sum(u$elig & a$elig),
    frac_reversal_among_eligible_both_adjusted = sum(a$rev & u$elig & a$elig) / sum(u$elig & a$elig),
    stringsAsFactors = FALSE)
}
trans_df <- do.call(rbind, trans_rows)
utils::write.csv(trans_df, file.path(out_dir, "reversal_transitions.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# meanshift_shrinkage.csv
# ---------------------------------------------------------------------------
ms_rows <- list()
qs <- function(x) as.numeric(stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE))
for (i in seq_along(conditions)) {
  cn <- conditions[[i]]
  ms <- res[[i]]$meanshift
  cohort_df_i <- res[[i]]$summary$n_cohort_level_residual_df
  for (ft in c("all", "mRNA", "miRNA")) {
    sel <- if (ft == "all") rep(TRUE, length(ms$ratio_w)) else ms$feature_type == ft
    qw <- qs(ms$ratio_w[sel]); qu <- qs(ms$ratio_u[sel])
    ms_rows[[length(ms_rows) + 1]] <- data.frame(
      condition = cond_id[i], block = cn$block, covset = cn$covset, mode = cn$mode, feature_type = ft, n_features = sum(sel),
      n_cohort_level_residual_df = cohort_df_i,
      ratio_between_var_weighted_q25 = qw[1], ratio_between_var_weighted_median = qw[2], ratio_between_var_weighted_q75 = qw[3],
      ratio_between_var_unweighted_q25 = qu[1], ratio_between_var_unweighted_median = qu[2], ratio_between_var_unweighted_q75 = qu[3],
      median_eta2_before = stats::median(ms$eta2_before[sel], na.rm = TRUE),
      median_eta2_after = stats::median(ms$eta2_after[sel], na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
}
ms_df <- do.call(rbind, ms_rows)
utils::write.csv(ms_df, file.path(out_dir, "meanshift_shrinkage.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# exemplars.csv
# ---------------------------------------------------------------------------
ex_rows <- list()
for (i in seq_along(conditions)) {
  cn <- conditions[[i]]
  ex <- res[[i]]$exemplars
  base_ex <- res[[baseline_id(cn$block)]]$exemplars
  ex$condition <- cond_id[i]
  ex$block <- cn$block
  ex$covset <- cn$covset
  ex$mode <- cn$mode
  key <- paste(ex$gene, ex$mir)
  bkey <- paste(base_ex$gene, base_ex$mir)
  mb <- match(key, bkey)
  ex$unadj_r_global <- base_ex$r_global[mb]
  ex$unadj_simpson_reversal <- base_ex$simpson_reversal[mb]
  ex$unadj_classification <- base_ex$classification[mb]
  ex$unadj_I2 <- base_ex$I2[mb]
  ex$unadj_n_neg_any <- base_ex$n_neg_any[mb]
  ex$classification_changed <- ex$classification != ex$unadj_classification
  ex$reversal_flag_changed <- ex$simpson_reversal != ex$unadj_simpson_reversal
  ex_rows[[i]] <- ex
}
ex_df <- do.call(rbind, ex_rows)
ex_df <- ex_df[, c("condition", "block", "covset", "mode", "set", "gene", "mir", setdiff(colnames(ex_df), c("condition", "block", "covset", "mode", "set", "gene", "mir")))]
utils::write.csv(ex_df, file.path(out_dir, "exemplars.csv"), row.names = FALSE)

# Per-pair data for the baselines and the headline conditions (for the figure and for the report).
key_conditions <- c(baseline_id("all"), baseline_id("matched"),
                    "matched|leuk_purity|pooled", "matched|leuk_purity|within", "matched|leuk_purity|within_between",
                    "matched|cibersort_purity|pooled", "matched|cibersort_purity|within", "matched|cibersort_purity|within_between",
                    "all|estimate|pooled", "all|estimate|within", "all|marker|pooled", "all|marker|within")
key_conditions <- intersect(key_conditions, cond_id)
pair_level <- lapply(key_conditions, function(id) {
  pm <- res[[id]]$pm
  rs <- rev_sets(pm)
  list(r_global = pm$r_global, r_resid = pm$r_resid, I2 = pm$I2, k_valid = pm$k_valid,
       n_pos_eps0p05 = pm$n_pos_eps0p05, n_neg_eps0p05 = pm$n_neg_eps0p05, eligible = rs$elig, reversal = rs$rev)
})
names(pair_level) <- key_conditions
saveRDS(list(genes = genes, mirs = mirs, pair_gene = rep(genes, times = M), pair_mir = rep(mirs, each = G),
             conditions = pair_level, params = list(G = G, M = M, eps = eps_main, t = t_main, min_k = min_k)),
        file.path(out_dir, "pair_level_key_conditions.rds"))
toc("tables", t0)

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
t0 <- tic()
mode_show <- c("pooled", "within", "within_between")
mode_col <- c(pooled = "#d95f02", within = "#1b9e77", within_between = "#7570b3")
mode_pch <- c(pooled = 16, within = 17, within_between = 15)
mode_lab <- c(pooled = "pooled regression", within = "within-cohort (means kept)", within_between = "within + between cohort")
x_layout <- rbind(
  data.frame(block = "all", covset = c("none", "marker", "estimate"), x = 1:3, stringsAsFactors = FALSE),
  data.frame(block = "matched", covset = c("none", "leuk", "leuk_purity", "cibersort", "cibersort_purity", "marker", "estimate"), x = 5:11, stringsAsFactors = FALSE)
)
x_labels <- c(none = "unadjusted", marker = "marker scores", estimate = "ESTIMATE", leuk = "LF",
              leuk_purity = "LF + purity", cibersort = "CIBERSORT", cibersort_purity = "CIBERSORT + purity")

panel_rates <- function(metric, ylab, panel_label) {
  lo_col <- paste0(metric, "_ci_lo"); hi_col <- paste0(metric, "_ci_hi")
  yr <- range(c(rates_df[[metric]][rates_df$mode %in% c("none", mode_show)],
                rates_df[[lo_col]][rates_df$mode %in% c("none", mode_show)],
                rates_df[[hi_col]][rates_df$mode %in% c("none", mode_show)]), na.rm = TRUE)
  yr <- c(max(0, yr[1] - 0.05 * diff(yr)), yr[2] + 0.12 * diff(yr))
  plot(NA, xlim = c(0.4, 11.6), ylim = yr, xaxt = "n", xlab = "", ylab = ylab, las = 1)
  axis(1, at = x_layout$x, labels = x_labels[x_layout$covset], las = 2, cex.axis = 0.8)
  abline(v = 4, col = "grey70", lty = 3)
  for (b in c("all", "matched")) {
    base <- rates_df[rates_df$block == b & rates_df$covset == "none", ]
    xr <- range(x_layout$x[x_layout$block == b])
    segments(xr[1] - 0.4, base[[metric]], xr[2] + 0.4, base[[metric]], lty = 2, col = "grey40")
    rect(xr[1] - 0.4, base[[lo_col]], xr[2] + 0.4, base[[hi_col]], col = grDevices::adjustcolor("grey60", 0.25), border = NA)
  }
  for (i in seq_len(nrow(x_layout))) {
    b <- x_layout$block[i]; cs_name <- x_layout$covset[i]; x0 <- x_layout$x[i]
    if (cs_name == "none") {
      r <- rates_df[rates_df$block == b & rates_df$covset == "none", ]
      arrows(x0, r[[lo_col]], x0, r[[hi_col]], angle = 90, code = 3, length = 0.02, col = "black")
      points(x0, r[[metric]], pch = 22, bg = "black", col = "black", cex = 1.1)
    } else {
      for (m in seq_along(mode_show)) {
        r <- rates_df[rates_df$block == b & rates_df$covset == cs_name & rates_df$mode == mode_show[m], ]
        xx <- x0 + (m - 2) * 0.25
        arrows(xx, r[[lo_col]], xx, r[[hi_col]], angle = 90, code = 3, length = 0.02, col = mode_col[mode_show[m]])
        points(xx, r[[metric]], pch = mode_pch[mode_show[m]], col = mode_col[mode_show[m]], cex = 1.1)
      }
    }
  }
  usr <- par("usr")
  text(2, usr[4] - 0.04 * diff(usr[3:4]), sprintf("all samples (N = %s)", format(length(blocks$all$samples), big.mark = ",")), cex = 0.8)
  text(8, usr[4] - 0.04 * diff(usr[3:4]), sprintf("matched samples (N = %s)", format(length(blocks$matched$samples), big.mark = ",")), cex = 0.8)
  mtext(panel_label, side = 3, adj = 0, line = 0.6, font = 2, cex = 1.1)
}

set.seed(seed_fig)
grDevices::cairo_pdf(file.path(out_dir, "Figure_CompositionAdjusted.pdf"), width = 11, height = 8.5)
par(mfrow = c(2, 2), mar = c(7.8, 4.8, 2.4, 0.8), oma = c(0, 0, 0.6, 0), mgp = c(2.9, 0.7, 0), cex.lab = 0.95)
panel_rates("frac_simpson_eps0p05_t0p2", expression("Simpson reversal rate at " * abs(r[global]) >= 0.2), "A.")
legend("bottomleft", legend = c("unadjusted", mode_lab[mode_show]), pch = c(22, mode_pch[mode_show]),
       col = c("black", mode_col[mode_show]), pt.bg = c("black", NA, NA, NA), bty = "n", cex = 0.8)
panel_rates("median_I2", expression("Median " * I^2 * " across pairs"), "B.")

# Panel C: unadjusted vs adjusted r_global for the fullest matched adjustment (pooled mode).
cid <- "matched|cibersort_purity|pooled"
bid <- baseline_id("matched")
u <- rev_sets(res[[bid]]$pm); a <- rev_sets(res[[cid]]$pm)
cat_pair <- rep("other", P)
cat_pair[u$rev & a$rev] <- "persist"
cat_pair[u$rev & !a$rev] <- "vanish"
cat_pair[!u$rev & a$rev] <- "new"
show_n <- min(30000L, P)
idx_show <- sample.int(P, show_n)
# Make sure the rarer categories are visible: draw "other" first, then the rest.
cat_col <- c(other = grDevices::adjustcolor("grey55", 0.25), persist = "#d62728", vanish = "#1f77b4", new = "#ff7f0e")
plot(u$rg[idx_show], a$rg[idx_show], type = "n", las = 1,
     xlab = expression("Unadjusted " * r[global]), ylab = expression("Adjusted " * r[global]))
o_idx <- idx_show[cat_pair[idx_show] == "other"]
points(u$rg[o_idx], a$rg[o_idx], pch = 16, cex = 0.3, col = cat_col["other"])
for (cc in c("vanish", "new", "persist")) {
  ii <- idx_show[cat_pair[idx_show] == cc]
  if (length(ii) > 0) points(u$rg[ii], a$rg[ii], pch = 16, cex = 0.45, col = cat_col[cc])
}
abline(0, 1, col = "grey30", lty = 1)
abline(v = c(-t_main, t_main), h = c(-t_main, t_main), col = "grey60", lty = 3)
legend("topleft", bty = "n", cex = 0.8, pch = 16, col = cat_col[c("persist", "vanish", "new", "other")],
       legend = c("reversal persists", "reversal vanishes", "new reversal", "other pairs"))
mtext("C.", side = 3, adj = 0, line = 0.6, font = 2, cex = 1.1)

# Panel D: between-cohort variance of cohort means, after/before, matched block.
d_sets <- c("leuk", "leuk_purity", "cibersort", "cibersort_purity", "marker", "estimate")
d_modes <- c("pooled", "within_between")
box_data <- list(); box_col <- character(0); box_at <- numeric(0)
for (j in seq_along(d_sets)) for (m in seq_along(d_modes)) {
  id <- paste("matched", d_sets[j], d_modes[m], sep = "|")
  box_data[[length(box_data) + 1]] <- pmin(res[[id]]$meanshift$ratio_w, 1.5)
  box_col <- c(box_col, mode_col[d_modes[m]])
  box_at <- c(box_at, j + (m - 1.5) * 0.36)
}
boxplot(box_data, at = box_at, col = grDevices::adjustcolor(box_col, 0.5), border = box_col, outline = FALSE,
        xaxt = "n", las = 1, boxwex = 0.3, ylim = c(0, 1.2),
        ylab = expression(Var[between] * " of cohort means, adjusted / unadjusted"))
axis(1, at = seq_along(d_sets), labels = x_labels[d_sets], las = 2, cex.axis = 0.8)
abline(h = 1, lty = 2, col = "grey40")
legend("bottomleft", legend = mode_lab[d_modes], fill = grDevices::adjustcolor(mode_col[d_modes], 0.5), border = mode_col[d_modes], bty = "n", cex = 0.8)
mtext("D.", side = 3, adj = 0, line = 0.6, font = 2, cex = 1.1)
# overall title lives in the LaTeX caption, not in the image
dev.off()
toc("figure", t0)

# ---------------------------------------------------------------------------
# Runtime and console summary
# ---------------------------------------------------------------------------
total_s <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
rt <- c(sprintf("total_seconds\t%.1f", total_s),
        sprintf("%s\t%.1f", names(section_times), unlist(section_times)),
        sprintf("n_cores\t%d", n_cores), sprintf("B\t%d", B), sprintf("G\t%d", G), sprintf("M\t%d", M),
        sprintf("quick\t%s", QUICK), sprintf("R\t%s", R.version.string),
        sprintf("estimate_pkg\t%s", estimate_version), sprintf("date\t%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
writeLines(rt, file.path(out_dir, "runtime.txt"))

show <- rates_df[, c("condition", "n_samples", "n_cohorts", "frac_simpson_eps0p05_t0p2", "frac_simpson_eps0p05_t0p2_ci_lo",
                     "frac_simpson_eps0p05_t0p2_ci_hi", "median_I2", "frac_mixed_eps0p05", "frac_resid_flip", "median_R2_covariates")]
print(show, digits = 4, row.names = FALSE)
message(sprintf("Done in %.1f min. Outputs in %s", total_s / 60, normalizePath(out_dir)))
