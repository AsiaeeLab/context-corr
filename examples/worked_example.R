# Worked example for the context-corr reporting checklist (base R only).
#
# Run from the repo root:
#   Rscript examples/worked_example.R
#
# 1. Simulates 4 contexts x 150 samples and four feature pairs of known type:
#    a genuine Simpson reversal built from opposing context means, a globally
#    robust pair, a context-switching pair, and a null pair.
# 2. Calls context_report() on each and prints the five checklist items plus
#    the classification.
# 3. Shows residualize_matrix() and random_effects_adjacency() on a 20-feature
#    toy matrix whose pooled correlations are dominated by context means.

find_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  cand <- c(getwd(), if (length(f) == 1) file.path(dirname(normalizePath(f)), ".."))
  for (d in cand) {
    if (file.exists(file.path(d, "R", "contextcorr.R"))) return(normalizePath(d))
  }
  stop("Could not find R/contextcorr.R; run from the repo root.", call. = FALSE)
}
root <- find_root()
source(file.path(root, "R", "contextcorr.R"))

set.seed(20260903)
k <- 4
n_per <- 150
context <- factor(rep(paste0("C", seq_len(k)), each = n_per))
n <- length(context)

# One feature pair with chosen within-context correlations r_within[c] and
# chosen context means: x = mu_x[c] + e1,  y = mu_y[c] + r e1 + sqrt(1 - r^2) e2.
sim_pair <- function(r_within, mu_x, mu_y, context) {
  x <- y <- numeric(length(context))
  for (i in seq_along(levels(context))) {
    idx <- which(context == levels(context)[i])
    e1 <- rnorm(length(idx))
    e2 <- rnorm(length(idx))
    x[idx] <- mu_x[i] + e1
    y[idx] <- mu_y[i] + r_within[i] * e1 + sqrt(1 - r_within[i]^2) * e2
  }
  list(x = x, y = y)
}

pairs <- list(
  # Within every context r = -0.4, but both means rise together across contexts,
  # so the pooled correlation is strongly positive: a Simpson reversal.
  compositional = sim_pair(rep(-0.4, k), mu_x = c(0, 1.5, 3, 4.5), mu_y = c(0, 1.5, 3, 4.5), context),
  # The same r = +0.5 in every context, no mean shifts.
  globally_robust = sim_pair(rep(0.5, k), mu_x = rep(0, k), mu_y = rep(0, k), context),
  # r = +0.5 in two contexts and -0.5 in the other two: the pooled r is near 0.
  context_switching = sim_pair(c(0.5, 0.5, -0.5, -0.5), mu_x = rep(0, k), mu_y = rep(0, k), context),
  # Independent noise everywhere.
  indeterminate = sim_pair(rep(0, k), mu_x = rep(0, k), mu_y = rep(0, k), context)
)

for (nm in names(pairs)) {
  cat("\n==== Pair simulated as:", nm, "====\n")
  rpt <- context_report(pairs[[nm]]$x, pairs[[nm]]$y, context, eps = 0.05, method = "pearson")
  print_context_report(rpt)
}

# ---------------------------------------------------------------------------
# 20-feature toy matrix: two co-regulated modules (g1-g5 share latent factor f1,
# g6-g10 share f2) on top of context mean shifts that increase for g1-g3 and
# decrease for g4-g5. Within every context cor(g1, g4) is about +0.64; pooled
# over contexts it is strongly negative.
# ---------------------------------------------------------------------------
cat("\n==== 20-feature toy matrix (features x samples) ====\n")
p <- 20
ctx_i <- as.integer(context)
shift <- c(0, 1.5, 3, 4.5)
f1 <- rnorm(n)
f2 <- rnorm(n)
expr <- matrix(rnorm(p * n, sd = 0.6), nrow = p,
               dimnames = list(paste0("g", seq_len(p)), paste0("s", seq_len(n))))
expr[1:5, ] <- expr[1:5, ] + 0.8 * rep(f1, each = 5)
expr[6:10, ] <- expr[6:10, ] + 0.8 * rep(f2, each = 5)
expr[1:3, ] <- expr[1:3, ] + rep(shift[ctx_i], each = 3)
expr[4:5, ] <- expr[4:5, ] - rep(shift[ctx_i], each = 2)

res <- residualize_matrix(expr, context)
ctx_means <- sapply(levels(context), function(cx) rowMeans(res[, context == cx, drop = FALSE]))
cat("dim(res):", dim(res), "; max |per-context mean| after residualization:",
    format(max(abs(ctx_means)), digits = 3), "\n")

adj <- random_effects_adjacency(expr, context, min_n = 4)
show <- c("g1", "g4", "g6", "g11")
cat("\nPooled Pearson correlation (all samples):\n")
print(round(stats::cor(t(expr))[show, show], 2))
cat("\nRandom-effects pooled within-context correlation r_RE:\n")
print(round(adj$r_RE[show, show], 2))
cat("\nCorrelation of the mean-residualized matrix:\n")
print(round(stats::cor(t(res))[show, show], 2))
cat("\nI2 (%) per pair:\n")
print(round(100 * adj$I2[show, show]))

cat("\ncontext_report() for g1 vs g4:\n")
print_context_report(context_report(expr["g1", ], expr["g4", ], context))

# adj$r_RE is a plain symmetric matrix with unit diagonal, so it can replace a
# pooled correlation matrix directly, e.g. a WGCNA-style soft threshold:
#   A <- abs(adj$r_RE)^6
# and `res` can replace the raw matrix in any correlation-based network tool.
