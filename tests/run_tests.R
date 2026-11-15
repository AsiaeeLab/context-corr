# Plain-R tests for R/contextcorr.R (base R only; testthat is not needed).
#
# Run from the repo root:
#   Rscript tests/run_tests.R
# Exits with status 1 if any check fails.

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
source(file.path(root, "R", "sc_context_pipeline.R"))

n_pass <- 0L
n_fail <- 0L
check <- function(name, expr) {
  ok <- tryCatch({ expr; TRUE }, error = function(e) {
    cat("        ", conditionMessage(e), "\n")
    FALSE
  })
  cat(if (ok) "[PASS]" else "[FAIL]", name, "\n")
  if (ok) n_pass <<- n_pass + 1L else n_fail <<- n_fail + 1L
  invisible(ok)
}
near <- function(a, b, tol = 1e-10) isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))

# ---------------------------------------------------------------------------
# Shared simulation (same design as examples/worked_example.R).
# ---------------------------------------------------------------------------
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
set.seed(20260903)
k <- 4
context <- factor(rep(paste0("C", seq_len(k)), each = 150))
toy <- list(
  compositional = sim_pair(rep(-0.4, k), c(0, 1.5, 3, 4.5), c(0, 1.5, 3, 4.5), context),
  globally_robust = sim_pair(rep(0.5, k), rep(0, k), rep(0, k), context),
  context_switching = sim_pair(c(0.5, 0.5, -0.5, -0.5), rep(0, k), rep(0, k), context),
  indeterminate = sim_pair(rep(0, k), rep(0, k), rep(0, k), context)
)

# ---------------------------------------------------------------------------
# (a) Toy pairs receive the intended classification.
# ---------------------------------------------------------------------------
for (nm in names(toy)) {
  check(sprintf("(a) toy pair '%s' is classified as %s", nm, nm), {
    rpt <- context_report(toy[[nm]]$x, toy[[nm]]$y, context)
    stopifnot(identical(rpt$classification$label, nm))
  })
}
check("(a) compositional toy pair: pooled r > 0, every within-context r < 0, simpson flag TRUE, r_RE and r_resid < 0", {
  rpt <- context_report(toy$compositional$x, toy$compositional$y, context)
  stopifnot(rpt$r_global > 0.3, all(rpt$r_by_context < -0.2), isTRUE(rpt$simpson$flag),
            rpt$random_effects$r_RE < 0, rpt$r_resid < 0)
})
check("(a) context_switching toy pair: |r_global| small, I2 large, tau2 > 0", {
  rpt <- context_report(toy$context_switching$x, toy$context_switching$y, context)
  stopifnot(abs(rpt$r_global) < 0.15, rpt$heterogeneity$I2 > 0.9, rpt$random_effects$tau2 > 0)
})
check("(a) classify_pair precedence: a switching pair that also reverses the pooled sign is context_switching", {
  r <- c(A = 0.6, B = 0.6, C = -0.6, D = -0.6, E = -0.6)
  n <- rep(200, 5)
  cl <- classify_pair(r_global = 0.5, r, n)
  stopifnot(identical(cl$label, "context_switching"), isTRUE(cl$simpson$flag),
            cl$n_switch_pos == 2, cl$n_switch_neg == 3)
})
check("(a) classify_pair: compositional when pooled sign opposes the majority vote without switching", {
  r <- c(A = -0.2, B = -0.25, C = -0.15, D = -0.3)
  cl <- classify_pair(r_global = 0.4, r, rep(200, 4))
  stopifnot(identical(cl$label, "compositional"), cl$n_switch_pos == 0)
})
check("(a) classify_pair: globally_robust needs every significant context to share the pooled sign", {
  r <- c(A = 0.4, B = 0.35, C = 0.02, D = 0.5)
  cl <- classify_pair(r_global = 0.4, r, rep(200, 4))
  stopifnot(identical(cl$label, "globally_robust"), cl$n_sig == 3, cl$n_sig_disagree == 0)
  # one significant dissenter below switch_thresh -> indeterminate
  r2 <- c(A = 0.4, B = 0.35, C = -0.2, D = 0.5)
  cl2 <- classify_pair(r_global = 0.3, r2, rep(200, 4))
  stopifnot(identical(cl2$label, "indeterminate"), cl2$n_sig_disagree == 1)
})
check("(a) classify_pair: indeterminate when no context is individually significant or none is informative", {
  cl <- classify_pair(r_global = 0.06, c(0.1, 0.05, -0.02), c(20, 20, 20))
  stopifnot(identical(cl$label, "indeterminate"), cl$n_sig == 0)
  cl2 <- classify_pair(r_global = 0.5, c(0.9, 0.8), c(3, 2))
  stopifnot(identical(cl2$label, "indeterminate"), cl2$k_informative == 0)
  # A perfect correlation is invalid on the Fisher-z scale. With only one
  # other valid context, no substantive context label may be assigned.
  cl3 <- classify_pair(r_global = 0.5, c(perfect = 1, valid = -0.4), c(200, 200))
  stopifnot(identical(cl3$label, "indeterminate"), cl3$k_informative == 1,
            is.na(cl3$contexts["perfect", "p"]))
})
check("(a) contextcorr_fisher_z_ci matches z +/- 1.96/sqrt(n - 3) and the two-sided normal p", {
  ci <- contextcorr_fisher_z_ci(c(a = 0.3, b = 0.5, c = 0.9), c(50, 3, 10))
  crit <- stats::qnorm(0.975)
  stopifnot(near(ci$ci_lo[1], tanh(atanh(0.3) - crit / sqrt(47))),
            near(ci$ci_hi[1], tanh(atanh(0.3) + crit / sqrt(47))),
            near(ci$p[1], 2 * stats::pnorm(-atanh(0.3) * sqrt(47))),
            is.na(ci$p[2]), !ci$excludes_zero[2],
            ci$excludes_zero[3], identical(rownames(ci), c("a", "b", "c")))
})

# ---------------------------------------------------------------------------
# (b) random_effects_pool reproduces the literal DerSimonian-Laird formulas.
# ---------------------------------------------------------------------------
check("(b) random_effects_pool matches hand-coded DerSimonian-Laird (heterogeneous case, tau2 > 0)", {
  r <- c(0.30, 0.10, -0.20, 0.50, 0.05)
  n <- c(50, 120, 80, 30, 200)
  w <- n - 3
  z <- atanh(r)
  zbar <- sum(w * z) / sum(w)
  Q <- sum(w * (z - zbar)^2)
  df <- length(r) - 1
  denom <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / denom)
  wRE <- 1 / (1 / w + tau2)
  zRE <- sum(wRE * z) / sum(wRE)
  rRE <- tanh(zRE)
  re <- random_effects_pool(r, n)
  stopifnot(re$k == 5, tau2 > 0, near(re$tau2, tau2), near(re$z_RE, zRE), near(re$r_RE, rRE))
  het <- heterogeneity_Q_I2(r, n)
  stopifnot(near(het$Q, Q), het$df == df, near(het$I2, max(0, (Q - df) / Q)),
            near(het$p, stats::pchisq(Q, df, lower.tail = FALSE)))
})
check("(b) random_effects_pool: homogeneous contexts give tau2 = 0 and the fixed-effect pooled z", {
  r <- c(0.30, 0.31, 0.29)
  n <- rep(100, 3)
  w <- n - 3
  z <- atanh(r)
  re <- random_effects_pool(r, n)
  stopifnot(re$tau2 == 0, near(re$z_RE, sum(w * z) / sum(w)), near(re$r_RE, tanh(sum(w * z) / sum(w))))
})
check("(b) random_effects_pool ignores contexts with NA r, n < 4, or |r| = 1", {
  r <- c(0.30, 0.10, -0.20, 0.50, 0.05)
  n <- c(50, 120, 80, 30, 200)
  a <- random_effects_pool(r, n)
  b <- random_effects_pool(c(r, NA, 0.9, 1), c(n, 100, 3, 50))
  stopifnot(identical(a, b))
})

# ---------------------------------------------------------------------------
# (c) residualize_matrix gives per-context means of ~0.
# ---------------------------------------------------------------------------
check("(c) residualize_matrix: per-context means ~0, dimnames kept, agrees with mean_residualized_cor", {
  set.seed(1)
  ctx <- factor(rep(c("a", "b", "c"), times = c(20, 25, 15)))
  m <- matrix(rnorm(10 * 60), 10, 60, dimnames = list(paste0("f", 1:10), paste0("s", 1:60)))
  m <- m + rep(c(a = 5, b = -3, c = 1)[as.integer(ctx)], each = 10)
  res <- residualize_matrix(m, ctx)
  means <- sapply(levels(ctx), function(cx) rowMeans(res[, ctx == cx, drop = FALSE]))
  stopifnot(max(abs(means)) < 1e-12, identical(dim(res), dim(m)), identical(dimnames(res), dimnames(m)),
            near(stats::cor(res[1, ], res[2, ]), mean_residualized_cor(m[1, ], m[2, ], ctx)))
})
check("(c) residualize_matrix: drops NA-context samples, singleton contexts residualize to 0", {
  m <- matrix(1:12, 3, 4)
  ctx <- c("a", "a", NA, "b")
  res <- suppressMessages(residualize_matrix(m, ctx))
  stopifnot(ncol(res) == 3, all(res[, 3] == 0), all(abs(rowMeans(res[, 1:2])) < 1e-12))
})

# ---------------------------------------------------------------------------
# (d) random_effects_adjacency is symmetric with unit diagonal.
# ---------------------------------------------------------------------------
check("(d) random_effects_adjacency: symmetric, unit diagonal, dimnames, symmetric I2", {
  set.seed(2)
  ctx <- factor(rep(paste0("c", 1:5), each = 40))
  m <- matrix(rnorm(8 * 200), 8, 200, dimnames = list(paste0("g", 1:8), NULL))
  adj <- random_effects_adjacency(m, ctx)
  A <- adj$r_RE
  stopifnot(identical(dim(A), c(8L, 8L)), identical(rownames(A), paste0("g", 1:8)),
            identical(colnames(A), paste0("g", 1:8)), all(diag(A) == 1),
            max(abs(A - t(A))) == 0, all(is.finite(A[upper.tri(A)])),
            max(abs(adj$I2 - t(adj$I2)), na.rm = TRUE) == 0, all(is.na(diag(adj$I2))),
            all(adj$k[upper.tri(adj$k)] == 5))
})
check("(d) random_effects_adjacency entries equal random_effects_pool / heterogeneity_Q_I2 on per-context cor", {
  set.seed(3)
  ctx <- factor(rep(paste0("c", 1:4), times = c(30, 50, 3, 70)))  # one context below min_n
  m <- matrix(rnorm(6 * 153), 6, 153)
  adj <- random_effects_adjacency(m, ctx, min_n = 4)
  stopifnot(identical(names(adj$n_by_context), c("c1", "c2", "c4")))
  use <- levels(ctx)[table(ctx) >= 4]
  r_c <- sapply(use, function(cx) stats::cor(m[2, ctx == cx], m[5, ctx == cx]))
  n_c <- as.integer(table(ctx)[use])
  stopifnot(near(adj$r_RE["f2", "f5"], random_effects_pool(r_c, n_c)$r_RE),
            near(adj$tau2["f2", "f5"], random_effects_pool(r_c, n_c)$tau2),
            near(adj$I2["f2", "f5"], heterogeneity_Q_I2(r_c, n_c)$I2))
  sub <- random_effects_adjacency(m, ctx, features = c(2, 5))
  stopifnot(identical(dim(sub$r_RE), c(2L, 2L)), near(sub$r_RE[1, 2], adj$r_RE["f2", "f5"]))
})

# ---------------------------------------------------------------------------
# context_report structure and internal consistency.
# ---------------------------------------------------------------------------
check("context_report returns the five checklist items consistent with the component functions", {
  x <- toy$compositional$x
  y <- toy$compositional$y
  ctx3 <- as.character(context)
  ctx3[1:3] <- "tiny"  # a context with n = 3: r computed but excluded (min_n = 4)
  rpt <- context_report(x, y, ctx3)
  stopifnot(all(c("r_global", "r_by_context", "n_by_context", "heterogeneity", "random_effects",
                  "r_resid", "simpson", "classification") %in% names(rpt)),
            near(rpt$r_global, stats::cor(x[ctx3 != "tiny"], y[ctx3 != "tiny"])),
            rpt$n == sum(ctx3 != "tiny"),
            rpt$n_by_context[["tiny"]] == 3, is.na(rpt$r_by_context[["tiny"]]),
            !("tiny" %in% rpt$contexts_used), rpt$k == 4,
            identical(rpt$heterogeneity, heterogeneity_Q_I2(rpt$r_by_context, rpt$n_by_context)),
            identical(rpt$random_effects, random_effects_pool(rpt$r_by_context, rpt$n_by_context)),
            near(rpt$r_resid, mean_residualized_cor(x[ctx3 != "tiny"], y[ctx3 != "tiny"], ctx3[ctx3 != "tiny"])),
            identical(rpt$simpson$flag, simpson_flag(rpt$r_global, rpt$r_by_context)$flag),
            identical(rpt$classification$simpson$flag, rpt$simpson$flag),
            rpt$classification$label %in% c("globally_robust", "compositional", "context_switching", "indeterminate"))
  out <- capture.output(print_context_report(rpt))
  stopifnot(any(grepl("^Classification: compositional", out)))

  # Regression guard for the retained-population convention used by the
  # single-cell pipelines: with exactly one retained context, centering cannot
  # change a correlation and the residualized-flip rate must be zero.
  ctx_one <- factor(c(rep("retained", 8), rep("dropped", 3)))
  m_one <- rbind(
    x = seq_len(11),
    y = c(1, 2, 4, 3, 5, 7, 6, 8, 11, 9, 10)
  )
  p_one <- data.frame(pair_id = 1L, idx_x = 1L, idx_y = 2L,
                      feature_x = "x", feature_y = "y")
  one <- suppressMessages(sc_compute_pair_statistics_fast_pearson(
    m_one, m_one, p_one, ctx_one, n_min = 4, validate_n = 1
  ))$summary_pairs
  stopifnot(near(one$r_global, one$r_resid, tol = 1e-10),
            mean(sign(one$r_global) != sign(one$r_resid)) == 0)
})

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
