#!/usr/bin/env Rscript
# Supplement to tcga_checks.R:
#   (i) all-sample vs primary-only comparison with SKCM (353 of 451 samples metastatic) excluded from both;
#  (ii) repeated TCGA participants (same 12-character barcode) within the analyzed matrix.
# Outputs: results/revision/sampletype_fisher_checks/tcga_sampletype_metrics_noSKCM.csv, tcga_repeated_participants.csv
suppressPackageStartupMessages({ source("00-paths.R"); source(file.path("R", "contextcorr.R")); library(parallel) })
N_CORES <- 16L; EPS <- 0.05; THR <- 0.2; K_MIN <- 2L; BOOT_CHUNKS <- 20L; SEED <- 20261111L
out_dir <- file.path("results", "revision", "sampletype_fisher_checks")
ex <- parse(file.path("scripts", "revision", "uniform_table.R"))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
                  identical(e[[3]][[1]], as.name("function"))) eval(e, envir = globalenv())
load(file.path(paths$mirtcga, "matchedNoGBM.rda"))
code <- substr(colnames(mRNAdata), 14, 15); cancerType <- as.character(cancerType)
G <- 5000; M <- 200
top_by_var <- function(A, n) names(sort(apply(A, 1, var), decreasing = TRUE))[seq_len(n)]
genes <- top_by_var(mRNAdata, G); mirs <- top_by_var(miRdata, M)
run <- function(keep, label, seed) {
  X <- mRNAdata[genes, keep]; Y <- miRdata[mirs, keep]; ct <- cancerType[keep]
  tab <- table(ct); retained <- names(tab)[tab >= 20]
  pop <- which(ct %in% retained); ctx <- droplevels(factor(ct[pop], levels = retained))
  rg <- as.vector(cor_rows_product(X[, pop], Y[, pop]))
  R <- sapply(retained, function(cc) { idx <- pop[ctx == cc]; as.vector(cor_rows_product(X[, idx, drop = FALSE], Y[, idx, drop = FALSE])) })
  rr <- as.vector(cor_rows_product(residualize_by_group(X[, pop], ctx), residualize_by_group(Y[, pop], ctx)))
  n_ctx <- as.integer(tab[retained]); names(n_ctx) <- retained
  res <- compute_metrics(rg, rr, R, n_ctx, boot_B = 1000L, boot_seed = seed)
  ci <- res$ci
  data.frame(variant = label, n_samples = length(pop), n_cohorts = length(retained), t(res$metrics), check.names = FALSE,
             simpson_ci_lo = ci$ci95_lo_pct[1], simpson_ci_hi = ci$ci95_hi_pct[1])
}
noskcm <- cancerType != "SKCM"
out <- rbind(run(which(noskcm), "all_samples_noSKCM", SEED + 51L),
             run(which(noskcm & code == "01"), "primary_tumor_only_noSKCM", SEED + 61L))
write.csv(out, file.path(out_dir, "tcga_sampletype_metrics_noSKCM.csv"), row.names = FALSE)
pid <- substr(colnames(mRNAdata), 1, 12)
rep_all <- table(pid); rep_p <- table(pid[code == "01"])
part <- data.frame(analysis = c("all samples", "primary tumors only"),
                   n_samples = c(length(pid), sum(code == "01")),
                   n_participants = c(length(rep_all), length(rep_p)),
                   n_participants_with_2plus_samples = c(sum(rep_all >= 2), sum(rep_p >= 2)),
                   n_samples_from_repeated_participants = c(sum(rep_all[rep_all >= 2]), sum(rep_p[rep_p >= 2])))
by_type <- table(code[pid %in% names(rep_all)[rep_all >= 2]])
write.csv(part, file.path(out_dir, "tcga_repeated_participants.csv"), row.names = FALSE)
print(out[, c("variant", "n_samples", "m5_pct_simpson_reversal_ge_0.2", "m3_median_I2", "simpson_ci_lo", "simpson_ci_hi")]); print(part); print(by_type)
