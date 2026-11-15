# Regenerate the GTEx pooled-versus-residualized panel using the 20 retained
# SMTSD tissues for both feature selection and correlation estimation.

source("00-paths.R")
suppressPackageStartupMessages({library(parallel)})
N_CORES <- max(1L, min(12L, parallel::detectCores() - 70L))
out_png <- "figures/fig_gtex_resid.png"
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

rowcor_pairs <- function(X, Y, ix, iy, chunk = 20000L) {
  zx <- suppressWarnings(t(scale(t(X))))
  zy <- if (identical(X, Y)) zx else suppressWarnings(t(scale(t(Y))))
  denom <- ncol(X) - 1
  out <- numeric(length(ix))
  for (a in seq(1L, length(ix), by = chunk)) {
    b <- min(length(ix), a + chunk - 1L); ii <- a:b
    out[ii] <- rowSums(zx[ix[ii], , drop = FALSE] * zy[iy[ii], , drop = FALSE]) / denom
  }
  out[!is.finite(out)] <- NA_real_
  out
}
sample_pairs <- function(G, P) {
  i <- sample.int(G, size = 2 * P, replace = TRUE)
  j <- sample.int(G, size = 2 * P, replace = TRUE)
  ok <- i != j; i <- i[ok]; j <- j[ok]
  a <- pmin(i, j); b <- pmax(i, j)
  u <- !duplicated(paste(a, b, sep = "_")); a <- a[u]; b <- b[u]
  if (length(a) < P) stop("Not enough unique pairs sampled.")
  list(i = a[1:P], j = b[1:P])
}

f_gtex <- file.path(paths$mirtcga, "xena-gtex.Rda")
load(f_gtex)
Xg <- as.matrix(gtex$mrna); rm(gtex)
attr <- utils::read.delim(Sys.getenv("GTEX_ATTRIBUTES", file.path(paths$mirtcga, "GTEx_SampleAttributes_v8.txt")),
                          stringsAsFactors = FALSE, quote = "")
smtsd <- attr$SMTSD[match(gsub("\\.", "-", colnames(Xg)), attr$SAMPID)]
tis_tab <- sort(table(smtsd), decreasing = TRUE)
keep_tissues <- names(tis_tab)[tis_tab >= 150L]
keep_tissues <- keep_tissues[seq_len(min(20, length(keep_tissues)))]
pop <- which(smtsd %in% keep_tissues)
message(sprintf("tissues: %d, samples in population: %d", length(keep_tissues), length(pop)))

G_gtex <- 1500L; P_gtex <- 200000L
gvar_20 <- apply(Xg[, pop, drop = FALSE], 1, var)
genes <- names(sort(gvar_20, decreasing = TRUE))[1:G_gtex]
set.seed(1)
pr <- sample_pairs(G_gtex, P_gtex); ii <- pr$i; jj <- pr$j

Xs <- Xg[genes, , drop = FALSE]; rm(Xg); gc()
# residualize within each retained tissue; population = the 20 tissues
Xr <- Xs
for (tt in keep_tissues) {
  idx <- which(smtsd == tt)
  Xr[, idx] <- Xs[, idx, drop = FALSE] - rowMeans(Xs[, idx, drop = FALSE])
}
r_global <- rowcor_pairs(Xs[, pop, drop = FALSE], Xs[, pop, drop = FALSE], ii, jj)
r_resid  <- rowcor_pairs(Xr[, pop, drop = FALSE], Xr[, pop, drop = FALSE], ii, jj)

flip <- sign(r_global) != sign(r_resid) & r_global != 0 & r_resid != 0
rate_all <- 100 * mean(flip, na.rm = TRUE)
rate_02  <- 100 * mean(flip[abs(r_global) >= 0.2], na.rm = TRUE)
message(sprintf("resid flip: all %.4f%%  |r_global|>=0.2 %.4f%%", rate_all, rate_02))

set.seed(1)
ok <- is.finite(r_global) & is.finite(r_resid)
ix <- which(ok); if (length(ix) > 200000) ix <- sample(ix, 200000)

grDevices::png(out_png, width = 7, height = 6, units = "in", res = 300)
par(mar = c(4.6, 4.8, 3.0, 1.0), font.main = 2, mgp = c(2.9, 0.8, 0))
plot(r_global[ix], r_resid[ix], pch = 16, col = rgb(0, 0, 0, 0.05),
     cex.lab = 1.4, cex.axis = 1.2, cex.main = 1.4,
     xlab = expression("Pooled correlation " * italic(r)[global]),
     ylab = expression("Tissue-residualized correlation " * italic(r)[resid]),
     main = "GTEx: Pooled vs Tissue-Residualized Correlation")
abline(h = 0, v = 0, lwd = 2, col = "grey70")
abline(a = 0, b = 1, lwd = 2, col = "red")
dev.off()
message("wrote ", out_png)
write.csv(data.frame(quantity = c("resid_flip_all_pct", "resid_flip_ge_0.2_pct",
                                  "n_pairs", "n_tissues", "n_samples_population"),
                     value = c(rate_all, rate_02, P_gtex, length(keep_tissues), length(pop))),
          "results/revision/uniform_table/fig4B_regeneration.csv", row.names = FALSE)
