# Which pairs pass BH in the within-cohort permutation null (permutation 1 of tcga_checks.R).
# Run from the repo root: Rscript scripts/revision/null_tail.R
suppressPackageStartupMessages({ source("00-paths.R"); source("R/contextcorr.R") })
ex <- parse("scripts/revision/uniform_table.R")
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) eval(e)
load(file.path(paths$mirtcga, "matchedNoGBM.rda")); cancerType <- as.character(cancerType)
G <- 5000; M <- 200
top <- function(A, n) names(sort(apply(A, 1, var), decreasing = TRUE))[seq_len(n)]
X <- mRNAdata[top(mRNAdata, G), ]; Y <- miRdata[top(miRdata, M), ]
coh <- sort(unique(cancerType)); n_ctx <- as.integer(table(cancerType)[coh])
set.seed(20261111L + 100L)   # same seed as permutation 1 in tcga_checks.R
Rp <- matrix(NA_real_, G * M, length(coh)); Zf <- matrix(NA_real_, G * M, length(coh))
for (j in seq_along(coh)) { idx <- which(cancerType == coh[j]); Yc <- t(apply(Y[, idx], 1, function(v) v[sample.int(length(v))]))
  Rp[, j] <- as.vector(cor_rows_product(X[, idx], Yc))
  zx <- rowMeans(X[, idx] <= log2(10) + 1e-9); zy <- rowMeans(Y[, idx] <= 1e-12)
  Zf[, j] <- pmax(rep(zx, times = M), rep(zy, each = G)) }
valid <- is.finite(Rp) & abs(Rp) < 1; Rv <- Rp; Rv[!valid] <- 0; z <- atanh(Rv); z[!valid] <- 0
W <- matrix(rep(n_ctx - 3, each = G * M), G * M) * valid
Q <- rowSums(W * z^2) - rowSums(W * z)^2 / rowSums(W); df <- rowSums(valid) - 1
p <- pchisq(Q, df, lower.tail = FALSE); q <- p.adjust(p, "BH"); hit <- which(q < 0.05)
cat("BH hits:", length(hit), "\n")
# For each hit: the cohort contributing most to Q, its |r_c|, its zero fraction
contrib <- W * (z - rowSums(W * z) / rowSums(W))^2
jmax <- apply(contrib[hit, , drop = FALSE], 1, which.max)
zf <- Zf[cbind(hit, jmax)]; rr <- abs(Rp[cbind(hit, jmax)])
cat("top-cohort zero fraction quantiles:", round(quantile(zf, c(.1,.25,.5,.75,.9)), 3), "\n")
cat("share with top-cohort zero fraction >= 0.9:", mean(zf >= 0.9), " >=0.5:", mean(zf >= 0.5), "\n")
cat("top-cohort |r_c| quantiles:", round(quantile(rr, c(.1,.5,.9)), 3), "\n")
cat("background share of cells with zero fraction >= 0.9:", mean(Zf[valid] >= 0.9), "\n")
