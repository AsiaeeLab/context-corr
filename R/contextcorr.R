# Context-aware correlation helpers
#
# This repo is not packaged as an R package; this file provides a small, reusable
# API that scripts and Rmds can `source()` to implement “correlation in context”
# reporting consistently.

contextcorr_fisher_z <- function(r) {
  r <- pmin(pmax(r, -0.999999), 0.999999)
  0.5 * log((1 + r) / (1 - r))
}

contextcorr_inv_fisher_z <- function(z) {
  e <- exp(2 * z)
  (e - 1) / (e + 1)
}

context_correlations <- function(expr_x, expr_y, context, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  stopifnot(length(expr_x) == length(expr_y), length(expr_x) == length(context))

  ok <- is.finite(expr_x) & is.finite(expr_y) & !is.na(context)
  expr_x <- expr_x[ok]
  expr_y <- expr_y[ok]
  context <- context[ok]

  context <- as.factor(context)
  levels_ctx <- levels(context)

  r_global <- suppressWarnings(stats::cor(expr_x, expr_y, method = method))

  r_by_context <- setNames(rep(NA_real_, length(levels_ctx)), levels_ctx)
  n_by_context <- setNames(integer(length(levels_ctx)), levels_ctx)

  for (ctx in levels_ctx) {
    idx <- which(context == ctx)
    n_by_context[ctx] <- length(idx)
    if (length(idx) < 3) next
    r_by_context[ctx] <- suppressWarnings(stats::cor(expr_x[idx], expr_y[idx], method = method))
  }

  list(
    r_global = r_global,
    r_by_context = r_by_context,
    n_by_context = n_by_context,
    k = sum(is.finite(r_by_context)),
    df = max(sum(is.finite(r_by_context)) - 1, 0)
  )
}

heterogeneity_Q_I2 <- function(r_by_context, n_by_context) {
  stopifnot(length(r_by_context) == length(n_by_context))

  r <- as.numeric(r_by_context)
  n <- as.numeric(n_by_context)

  valid <- is.finite(r) & is.finite(n) & n >= 4 & abs(r) < 1
  r <- r[valid]
  n <- n[valid]

  k <- length(r)
  df <- max(k - 1, 0)
  if (k < 2) {
    return(list(k = k, df = df, Q = NA_real_, p = NA_real_, I2 = NA_real_))
  }

  w <- n - 3
  z <- contextcorr_fisher_z(r)
  z_bar <- sum(w * z) / sum(w)
  Q <- sum(w * (z - z_bar) ^ 2)

  p <- stats::pchisq(Q, df = df, lower.tail = FALSE)
  I2 <- if (is.finite(Q) && Q > 0) max(0, (Q - df) / Q) else NA_real_

  list(k = k, df = df, Q = Q, p = p, I2 = I2)
}

random_effects_pool <- function(r_by_context, n_by_context) {
  stopifnot(length(r_by_context) == length(n_by_context))

  r <- as.numeric(r_by_context)
  n <- as.numeric(n_by_context)

  valid <- is.finite(r) & is.finite(n) & n >= 4 & abs(r) < 1
  r <- r[valid]
  n <- n[valid]

  k <- length(r)
  if (k < 2) {
    return(list(k = k, tau2 = NA_real_, z_RE = NA_real_, r_RE = NA_real_))
  }

  w <- n - 3
  z <- contextcorr_fisher_z(r)
  z_bar <- sum(w * z) / sum(w)
  Q <- sum(w * (z - z_bar) ^ 2)
  df <- k - 1

  denom <- sum(w) - sum(w ^ 2) / sum(w)
  tau2 <- max(0, (Q - df) / denom)

  wRE <- 1 / (1 / w + tau2)
  zRE <- sum(wRE * z) / sum(wRE)
  rRE <- contextcorr_inv_fisher_z(zRE)

  list(k = k, tau2 = tau2, z_RE = zRE, r_RE = rRE)
}

mean_residualized_cor <- function(expr_x, expr_y, context, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  stopifnot(length(expr_x) == length(expr_y), length(expr_x) == length(context))

  ok <- is.finite(expr_x) & is.finite(expr_y) & !is.na(context)
  expr_x <- expr_x[ok]
  expr_y <- expr_y[ok]
  context <- as.factor(context[ok])

  x_res <- expr_x - ave(expr_x, context, FUN = mean)
  y_res <- expr_y - ave(expr_y, context, FUN = mean)

  suppressWarnings(stats::cor(x_res, y_res, method = method))
}

simpson_flag <- function(r_global, r_by_context, eps = 0.05, majority = c("count", "median")) {
  majority <- match.arg(majority)

  r <- as.numeric(r_by_context)
  valid <- is.finite(r) & abs(r) >= eps
  r_use <- r[valid]

  n_pos <- sum(r_use > 0)
  n_neg <- sum(r_use < 0)
  n_used <- length(r_use)

  majority_sign <- 0L
  if (n_used > 0) {
    if (majority == "count") {
      majority_sign <- sign(n_pos - n_neg)
    } else {
      majority_sign <- sign(stats::median(r_use))
    }
  }

  global_sign <- sign(r_global)
  flag <- is.finite(r_global) && global_sign != 0 && majority_sign != 0 && (global_sign != majority_sign)

  list(
    flag = flag,
    r_global = r_global,
    global_sign = global_sign,
    majority_sign = majority_sign,
    n_pos = n_pos,
    n_neg = n_neg,
    n_used = n_used,
    eps = eps,
    majority = majority
  )
}

