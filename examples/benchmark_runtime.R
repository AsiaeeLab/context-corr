# Wall-clock timing of the reporting functions at the paper's scale
# (31 contexts, 8,890 samples), on synthetic data. Base R only.
#
# Run from the repo root:
#   Rscript examples/benchmark_runtime.R

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

set.seed(1)
n <- 8890
k <- 31
# Unequal context sizes (smallest 36, like the smallest TCGA cohort), summing to n.
raw <- stats::rgamma(k, shape = 2)
sizes <- pmax(36L, as.integer(floor(raw / sum(raw) * (n - 36 * k)))) 
sizes[k] <- sizes[k] + (n - sum(sizes))
context <- factor(rep(paste0("ctx", seq_len(k)), times = sizes))
shift <- stats::rnorm(k)

x <- stats::rnorm(n) + shift[as.integer(context)]
y <- 0.3 * x + stats::rnorm(n)

cat(R.version.string, "\n")
cat("platform:", R.version$platform, "\n")
cat("n =", n, "samples, k =", k, "contexts (sizes", min(sizes), "to", max(sizes), ")\n\n")

t_pair <- system.time(rpt <- context_report(x, y, context))
cat(sprintf("context_report, one pair:                 %.3f s elapsed\n", t_pair[["elapsed"]]))

n_rep <- 100
t_many <- system.time(for (i in seq_len(n_rep)) context_report(x, y, context))
cat(sprintf("context_report, %d pairs in a loop:       %.3f s elapsed (%.4f s per pair)\n",
            n_rep, t_many[["elapsed"]], t_many[["elapsed"]] / n_rep))

p <- 200
expr <- matrix(stats::rnorm(p * n), nrow = p, dimnames = list(paste0("g", seq_len(p)), NULL))
expr <- expr + outer(stats::rnorm(p), shift[as.integer(context)])
t_res <- system.time(res <- residualize_matrix(expr, context))
cat(sprintf("residualize_matrix, %d x %d:             %.3f s elapsed\n", p, n, t_res[["elapsed"]]))

t_adj <- system.time(adj <- random_effects_adjacency(expr, context, min_n = 4))
cat(sprintf("random_effects_adjacency, %d features:    %.3f s elapsed (%d pairs)\n",
            p, t_adj[["elapsed"]], p * (p - 1) / 2))
