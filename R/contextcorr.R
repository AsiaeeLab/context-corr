# Context-aware correlation helpers
#
# This repository is not packaged as an R package; this file provides a small,
# reusable interface that scripts and Rmds can `source()` to implement
# “correlation in context”
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


# ---------------------------------------------------------------------------
# Reporting-checklist extensions (BMC Genomics revision).
#
# Conventions shared by everything below (they match the paper):
# - Sign tolerance eps = 0.05: a context votes +1 if r_c >= eps, -1 if
#   r_c <= -eps, and 0 otherwise. Mixed-sign = both votes present.
# - Simpson reversal = sign(r_global) differs from the majority context vote
#   (simpson_flag with majority = "count").
# - Residualized correlation = correlation after subtracting per-context means
#   from both variables (mean_residualized_cor / residualize_matrix).
# - Fisher-z interval for r_c: z_c +/- qnorm(1 - alpha / 2) / sqrt(n_c - 3),
#   mapped back to the r scale; with alpha = 0.05 this is the 1.96 interval.
# - Random-effects pooling = DerSimonian-Laird on Fisher z with w = n - 3
#   (random_effects_pool).
# ---------------------------------------------------------------------------

contextcorr_fisher_z_ci <- function(r, n, alpha = 0.05) {
  # Per-context Fisher-z interval and two-sided test.
  # Returns a data.frame with one row per context: r, n, z, se, ci_lo, ci_hi
  # (r scale), p, and `excludes_zero` (equivalent to p < alpha).
  # Contexts with n < 4 have no finite standard error and get NA.
  stopifnot(length(r) == length(n))
  nm <- names(r)
  r <- as.numeric(r)
  n <- as.numeric(n)

  ok <- is.finite(r) & is.finite(n) & n >= 4 & abs(r) < 1
  z <- se <- lo <- hi <- p <- rep(NA_real_, length(r))
  z[ok] <- contextcorr_fisher_z(r[ok])
  se[ok] <- 1 / sqrt(n[ok] - 3)
  crit <- stats::qnorm(1 - alpha / 2)
  lo[ok] <- contextcorr_inv_fisher_z(z[ok] - crit * se[ok])
  hi[ok] <- contextcorr_inv_fisher_z(z[ok] + crit * se[ok])
  p[ok] <- 2 * stats::pnorm(-abs(z[ok]) / se[ok])

  out <- data.frame(
    r = r, n = n, z = z, se = se, ci_lo = lo, ci_hi = hi, p = p,
    excludes_zero = ok & (lo > 0 | hi < 0),
    stringsAsFactors = FALSE
  )
  out$excludes_zero[is.na(out$excludes_zero)] <- FALSE
  if (!is.null(nm)) rownames(out) <- nm
  out
}

classify_pair <- function(r_global, r_by_context, n_by_context, eps = 0.05,
                          switch_thresh = 0.3, alpha = 0.05) {
  # Classify one feature pair from its pooled and within-context correlations.
  #
  # Labels (checked in this order; the first that applies wins):
  #   "context_switching": at least one context with r_c >= +switch_thresh and
  #       another with r_c <= -switch_thresh, both with Fisher-z p < alpha.
  #   "compositional": operational Simpson-reversal label, i.e. the pooled sign
  #       disagrees with the eps-tolerant majority of context signs
  #       (simpson_flag(..., majority = "count")). This label alone does not
  #       identify the mechanism; inspect r_resid/the covariance decomposition.
  #   "globally_robust": |r_global| >= eps, at least one context has a Fisher-z
  #       CI excluding 0, every such context has the pooled sign, and there is
  #       no reversal.
  #   "indeterminate": none of the above (no informative context, no context
  #       individually distinguishable from 0, |r_global| < eps, or a
  #       significant dissenting context that does not reach switch_thresh).
  #
  # Precedence: context_switching is tested before compositional because a
  # switching pair (strong opposite signs) frequently also produces a majority
  # reversal of the pooled sign; the switching label is the more specific
  # description. A switching pair can never be globally_robust, since it has
  # significant contexts of both signs. compositional is tested before
  # globally_robust because globally_robust requires no reversal.
  #
  # Only "informative" contexts (finite r_c, |r_c| < 1, and n_c >= 4, so a
  # finite Fisher-z standard error exists) enter the vote and the tests. At
  # least two informative contexts are required for a substantive label.
  #
  # Returns a list whose first element `label` is the class; the remaining
  # elements are the supporting facts (per-context table, counts, the
  # simpson_flag output, a one-line `reason`, and the thresholds used).
  stopifnot(length(r_by_context) == length(n_by_context))

  ci <- contextcorr_fisher_z_ci(r_by_context, n_by_context, alpha = alpha)
  informative <- is.finite(ci$r) & abs(ci$r) < 1 & is.finite(ci$n) & ci$n >= 4
  significant <- informative & is.finite(ci$p) & ci$p < alpha

  r_inf <- as.numeric(r_by_context)
  r_inf[!informative] <- NA_real_
  names(r_inf) <- names(r_by_context)

  global_sign <- if (is.finite(r_global) && abs(r_global) >= eps) sign(r_global) else 0
  sign_eps <- ifelse(informative & abs(ci$r) >= eps, sign(ci$r), 0)

  n_sig <- sum(significant)
  n_sig_agree <- sum(significant & sign(ci$r) == global_sign)
  n_sig_disagree <- n_sig - n_sig_agree
  n_switch_pos <- sum(significant & ci$r >= switch_thresh)
  n_switch_neg <- sum(significant & ci$r <= -switch_thresh)

  sf <- simpson_flag(r_global, r_inf, eps = eps, majority = "count")

  enough_contexts <- sum(informative) >= 2
  is_switching <- enough_contexts && n_switch_pos >= 1 && n_switch_neg >= 1
  is_compositional <- enough_contexts && isTRUE(sf$flag)
  is_robust <- enough_contexts && global_sign != 0 && n_sig >= 1 && n_sig_disagree == 0 && !is_compositional

  sign_chr <- function(s) c("-", "0", "+")[s + 2]
  if (is_switching) {
    label <- "context_switching"
    reason <- sprintf(
      "%d context(s) with r >= +%.2f and %d with r <= -%.2f, all with Fisher-z p < %g",
      n_switch_pos, switch_thresh, n_switch_neg, switch_thresh, alpha
    )
  } else if (is_compositional) {
    label <- "compositional"
    reason <- sprintf(
      "pooled sign (%s) disagrees with the majority of context votes (+%d / -%d at eps = %g)",
      sign_chr(sf$global_sign), sf$n_pos, sf$n_neg, eps
    )
  } else if (is_robust) {
    label <- "globally_robust"
    reason <- sprintf(
      "pooled sign (%s) agrees with all %d context(s) whose Fisher-z CI excludes 0 (%d informative); no reversal",
      sign_chr(global_sign), n_sig, sum(informative)
    )
  } else {
    label <- "indeterminate"
    reason <- if (sum(informative) < 2) {
      "fewer than two informative contexts (need finite |r_c| < 1 and n_c >= 4)"
    } else if (n_sig == 0) {
      sprintf("no context has a Fisher-z CI excluding 0 at alpha = %g", alpha)
    } else if (global_sign == 0) {
      sprintf("|r_global| < eps = %g, so there is no pooled sign to test, and no switching", eps)
    } else {
      sprintf(
        "%d significant context(s) oppose the pooled sign without a majority reversal, and no direction reaches |r| >= %.2f",
        n_sig_disagree, switch_thresh
      )
    }
  }

  contexts <- ci
  contexts$sign_eps <- as.integer(sign_eps)
  contexts$informative <- informative
  contexts$significant <- significant

  list(
    label = label,
    reason = reason,
    r_global = r_global,
    global_sign = as.integer(global_sign),
    k_informative = sum(informative),
    n_sig = n_sig,
    n_sig_agree = n_sig_agree,
    n_sig_disagree = n_sig_disagree,
    n_switch_pos = n_switch_pos,
    n_switch_neg = n_switch_neg,
    simpson = sf,
    contexts = contexts,
    eps = eps,
    switch_thresh = switch_thresh,
    alpha = alpha
  )
}

residualize_matrix <- function(expr, context) {
  # Subtract each feature's per-context mean.
  # expr: numeric matrix, features x samples; context: one label per sample.
  # Samples with a missing context are dropped (with a message), so the result
  # aligns with context[!is.na(context)]. Dimnames are preserved. A singleton
  # context residualizes to exactly 0. The output is what WGCNA / ARACNe /
  # GENIE3-style tools can consume in place of the raw matrix.
  expr <- as.matrix(expr)
  if (!is.numeric(expr)) stop("`expr` must be a numeric matrix (features x samples).", call. = FALSE)
  if (ncol(expr) != length(context)) stop("`context` must have one entry per column of `expr`.", call. = FALSE)

  keep <- !is.na(context)
  if (!all(keep)) {
    message("residualize_matrix: dropping ", sum(!keep), " sample(s) with missing context.")
    expr <- expr[, keep, drop = FALSE]
    context <- context[keep]
  }
  context <- droplevels(as.factor(context))

  out <- expr
  for (ctx in levels(context)) {
    idx <- which(context == ctx)
    mu <- rowMeans(expr[, idx, drop = FALSE], na.rm = TRUE)
    out[, idx] <- expr[, idx, drop = FALSE] - mu
  }
  out
}

random_effects_adjacency <- function(expr, context, features = NULL, min_n = 4,
                                     method = c("pearson", "spearman")) {
  # Random-effects (DerSimonian-Laird) pooled within-context correlation for
  # every feature pair, as a drop-in similarity matrix for network tools.
  #
  # expr: numeric matrix, features x samples; context: one label per sample;
  # features: optional row subset (names or indices); min_n: contexts with
  # fewer samples are excluded (random_effects_pool additionally requires
  # n_c >= 4 and |r_c| < 1 per pair, so values of min_n below 4 change nothing).
  #
  # Per-context correlations are computed once per context with stats::cor on
  # the whole feature block; pooling then loops over the p(p-1)/2 pairs and
  # calls random_effects_pool and heterogeneity_Q_I2, so cost is O(p^2 k).
  #
  # Returns a list: r_RE (p x p, symmetric, diagonal 1), I2, tau2 and k (p x p,
  # symmetric, diagonal NA), n_by_context (sizes of the contexts used),
  # method, min_n.
  method <- match.arg(method)
  expr <- as.matrix(expr)
  if (!is.numeric(expr)) stop("`expr` must be a numeric matrix (features x samples).", call. = FALSE)
  if (ncol(expr) != length(context)) stop("`context` must have one entry per column of `expr`.", call. = FALSE)
  if (!is.null(features)) {
    if (is.character(features) && !all(features %in% rownames(expr))) {
      stop("Some `features` are not row names of `expr`.", call. = FALSE)
    }
    expr <- expr[features, , drop = FALSE]
  }
  if (is.null(rownames(expr))) rownames(expr) <- paste0("f", seq_len(nrow(expr)))
  feats <- rownames(expr)
  p <- length(feats)
  if (p < 1) stop("`expr` has no features.", call. = FALSE)

  keep <- !is.na(context)
  expr <- expr[, keep, drop = FALSE]
  context <- droplevels(as.factor(context[keep]))
  n_all <- table(context)
  use <- names(n_all)[n_all >= min_n]
  k <- length(use)
  if (k < 2) stop("Fewer than two contexts with n >= min_n; nothing to pool.", call. = FALSE)
  n_use <- as.integer(n_all[use])
  names(n_use) <- use

  # Lower-triangle pair index (i > j), in the same order as R[lower.tri(R)].
  lt <- which(lower.tri(matrix(0, p, p)), arr.ind = TRUE)
  n_pairs <- nrow(lt)

  # k x n_pairs matrix of within-context correlations.
  r_ctx <- matrix(NA_real_, nrow = k, ncol = n_pairs)
  for (ci in seq_len(k)) {
    idx <- which(context == use[ci])
    R_c <- suppressWarnings(stats::cor(t(expr[, idx, drop = FALSE]), method = method))
    r_ctx[ci, ] <- R_c[lower.tri(R_c)]
  }

  r_RE <- matrix(NA_real_, p, p, dimnames = list(feats, feats))
  I2 <- tau2 <- r_RE
  k_used <- matrix(NA_integer_, p, p, dimnames = list(feats, feats))
  for (m in seq_len(n_pairs)) {
    i <- lt[m, 1]
    j <- lt[m, 2]
    r_m <- r_ctx[, m]
    re <- random_effects_pool(r_m, n_use)
    het <- heterogeneity_Q_I2(r_m, n_use)
    r_RE[i, j] <- r_RE[j, i] <- re$r_RE
    tau2[i, j] <- tau2[j, i] <- re$tau2
    I2[i, j] <- I2[j, i] <- het$I2
    k_used[i, j] <- k_used[j, i] <- as.integer(re$k)
  }
  diag(r_RE) <- 1

  list(
    r_RE = r_RE,
    I2 = I2,
    tau2 = tau2,
    k = k_used,
    n_by_context = n_use,
    method = method,
    min_n = min_n
  )
}

context_report <- function(expr_x, expr_y, context, eps = 0.05,
                           method = c("pearson", "spearman"), min_n = 4,
                           switch_thresh = 0.3, alpha = 0.05) {
  # The five-item reporting checklist for one feature pair, in one call:
  #   1. r_global: pooled correlation over samples in retained contexts.
  #   2. r_by_context, n_by_context: within-context correlations and sizes
  #      (contexts with n_c < min_n keep their n but have r set to NA and are
  #      excluded from every pooled or conditional estimand, mirroring the
  #      paper's retained-context analysis population).
  #   3. heterogeneity: Cochran's Q, df, p, I2 (heterogeneity_Q_I2).
  #   4. random_effects: DerSimonian-Laird k, tau2, z_RE, r_RE (random_effects_pool).
  #   5. r_resid: mean-residualized correlation over the retained contexts.
  # Plus simpson = simpson_flag(...) and classification = classify_pair(...).
  method <- match.arg(method)
  stopifnot(length(expr_x) == length(expr_y), length(expr_x) == length(context))

  ok <- is.finite(expr_x) & is.finite(expr_y) & !is.na(context)
  expr_x <- expr_x[ok]
  expr_y <- expr_y[ok]
  context <- droplevels(as.factor(context[ok]))

  cc <- context_correlations(expr_x, expr_y, context, method = method)
  r_ctx <- cc$r_by_context
  n_ctx <- cc$n_by_context

  use_ctx <- n_ctx >= min_n
  r_ctx_use <- r_ctx
  r_ctx_use[!use_ctx] <- NA_real_

  retained <- as.character(context) %in% names(n_ctx)[use_ctx]
  r_global <- NA_real_
  if (sum(retained) >= 2L) {
    r_global <- suppressWarnings(stats::cor(expr_x[retained], expr_y[retained],
                                             method = method))
  }

  r_ctx_valid <- r_ctx_use
  r_ctx_valid[!is.finite(r_ctx_valid) | abs(r_ctx_valid) >= 1] <- NA_real_

  het <- heterogeneity_Q_I2(r_ctx_valid, n_ctx)
  re <- random_effects_pool(r_ctx_valid, n_ctx)

  r_resid <- NA_real_
  if (any(use_ctx)) {
    ctx_ok <- as.character(context)
    ctx_ok[!(ctx_ok %in% names(n_ctx)[use_ctx])] <- NA_character_
    r_resid <- mean_residualized_cor(expr_x, expr_y, context = ctx_ok, method = method)
  }

  sf <- simpson_flag(r_global, r_ctx_valid, eps = eps, majority = "count")
  cl <- classify_pair(r_global, r_ctx_valid, n_ctx, eps = eps,
                      switch_thresh = switch_thresh, alpha = alpha)

  list(
    r_global = r_global,
    n = sum(retained),
    r_by_context = r_ctx_valid,
    n_by_context = n_ctx,
    k = sum(is.finite(r_ctx_valid)),
    contexts_used = names(n_ctx)[use_ctx],
    heterogeneity = het,
    random_effects = re,
    r_resid = r_resid,
    simpson = sf,
    classification = cl,
    settings = list(eps = eps, method = method, min_n = min_n,
                    switch_thresh = switch_thresh, alpha = alpha)
  )
}

print_context_report <- function(report, digits = 3) {
  # Compact console view of a context_report() result. Returns it invisibly.
  fmt <- function(x) ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
  s <- report$settings
  cat("Context-aware correlation report (", s$method, "; eps = ", s$eps,
      ", min_n = ", s$min_n, ")\n", sep = "")
  cat("1. Pooled correlation: r_global = ", fmt(report$r_global),
      " (n = ", report$n, ")\n", sep = "")
  cat("2. Within-context correlations (k = ", report$k, " contexts used):\n", sep = "")
  tab <- report$classification$contexts
  show <- data.frame(
    context = rownames(tab),
    n = tab$n,
    r = round(tab$r, digits),
    ci_lo = round(tab$ci_lo, digits),
    ci_hi = round(tab$ci_hi, digits),
    p = signif(tab$p, 2),
    vote = tab$sign_eps,
    stringsAsFactors = FALSE
  )
  print(show, row.names = FALSE)
  het <- report$heterogeneity
  cat("3. Heterogeneity: Q = ", fmt(het$Q), " on ", het$df, " df, p = ",
      if (is.finite(het$p)) format.pval(het$p, digits = 2) else "NA",
      ", I2 = ", if (is.finite(het$I2)) paste0(formatC(100 * het$I2, digits = 1, format = "f"), "%") else "NA",
      "\n", sep = "")
  re <- report$random_effects
  cat("4. Random-effects pooled correlation: r_RE = ", fmt(re$r_RE),
      " (tau2 = ", fmt(re$tau2), ", k = ", re$k, ")\n", sep = "")
  cat("5. Mean-residualized correlation: r_resid = ", fmt(report$r_resid), "\n", sep = "")
  sf <- report$simpson
  cat("Simpson reversal flag: ", sf$flag, " (pooled sign ", sf$global_sign,
      ", majority sign ", sf$majority_sign, "; votes +", sf$n_pos, " / -", sf$n_neg, ")\n", sep = "")
  cat("Classification: ", report$classification$label, "\n  ",
      report$classification$reason, "\n", sep = "")
  invisible(report)
}
