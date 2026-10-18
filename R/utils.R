# Internal helpers shared across the package. None of these are exported.

# Resolve a variance from `variance` and/or `sd` arguments.
#
# If both are supplied and disagree, `variance` wins and a warning is issued
# (the behavior of the d/p/q/rBetaMS() family).
.resolve_variance <- function(variance, sd) {
  if (!is.null(variance) && !is.null(sd)) {
    if (!isTRUE(all.equal(variance, sd^2))) {
      warning("Nonequivalent values of variance and sd specified. Using variance.")
    }
    return(variance)
  }
  if (!is.null(variance)) {
    return(variance)
  }
  if (!is.null(sd)) {
    return(sd^2)
  }
  stop("Either `variance` or `sd` must be specified.")
}

# Assemble the list-of-lists structure returned by the moment functions.
# `raw`, `central` and `standardized` are numeric vectors (orders 1..k);
# only the types requested are included, always in the order
# raw -> central -> standardized.
.moment_list <- function(types, raw, central, standardized) {
  out <- list()
  if (any(types == "raw")) out[["raw"]] <- as.list(raw())
  if (any(types == "central")) out[["central"]] <- as.list(central())
  if (any(types == "standardized")) out[["standardized"]] <- as.list(standardized())
  out
}

# Raw, central and standardized moments of a discrete distribution given its
# support and probability weights.
.pmf_moments <- function(support, weights, types, orders) {
  i <- seq_len(orders)
  mu <- sum(weights * support)
  central <- function() vapply(i, function(r) sum(weights * (support - mu)^r), numeric(1))
  .moment_list(
    types,
    raw = function() vapply(i, function(r) sum(weights * support^r), numeric(1)),
    central = central,
    standardized = function() central() / sqrt(sum(weights * (support - mu)^2))^i
  )
}

# Integrate `f` over [lower, upper] intersected with the support [l, u] of a
# four-parameter Beta density. Integrating only where the density is non-zero
# is both faster and more reliable than handing `integrate()` a range in which
# most of the integrand is identically zero.
.integrate_support <- function(f, l, u, lower = 0, upper = 1) {
  lo <- max(lower, l)
  hi <- min(upper, u)
  if (!(lo < hi)) {
    return(0)
  }
  stats::integrate(f, lower = lo, upper = hi)$value
}

# Probability mass function of the (four-parameter) Beta-Binomial distribution
# at the values `x`.
.betabinom_pmf <- function(x, N, l, u, alpha, beta) {
  valid <- is.finite(x) & x >= 0 & x <= N & x == round(x)
  out <- numeric(length(x))
  if (!any(valid)) {
    return(out)
  }
  if (l == 0 && u == 1) {
    # Closed form for the standard (two-parameter) Beta-Binomial.
    xv <- x[valid]
    out[valid] <- exp(lchoose(N, xv) + lbeta(xv + alpha, N - xv + beta) - lbeta(alpha, beta))
  } else {
    out[valid] <- vapply(x[valid], function(xi) {
      .integrate_support(function(p) stats::dbinom(xi, N, p) * dBeta.4P(p, l, u, alpha, beta),
                         l, u, 0, 1)
    }, numeric(1))
  }
  out
}

# Probability that a discrete variable on 0:N with probability masses `pmf`
# falls below `q` (lower.tail = TRUE, i.e. P[X < q]) or at or above `q`
# (lower.tail = FALSE, i.e. P[X >= q]). This is the convention used throughout
# the package, which differs from base R's P[X <= q].
.discrete_tail <- function(q, pmf, lower.tail = TRUE) {
  N <- length(pmf) - 1
  idx <- pmin(pmax(ceiling(q), 0), N + 1)
  if (lower.tail) {
    c(0, cumsum(pmf))[idx + 1]
  } else {
    c(rev(cumsum(rev(pmf))), 0)[idx + 1]
  }
}

# P[X < c] (lower.tail = TRUE) or P[X >= c] (lower.tail = FALSE) under Lord's
# two-term approximation to the compound Binomial distribution, in closed form.
#
# Summing the two-term pmf over 0..(c - 1) turns every dbinom() into a pbinom():
#   F(m) = pbinom(m, N, p) - k p (1 - p) [pbinom(m, N - 2, p)
#            - 2 pbinom(m - 1, N - 2, p) + pbinom(m - 2, N - 2, p)]
# with m = ceiling(c) - 1. The bracketed term sums to zero over the full
# support, so the same expression with lower.tail = FALSE gives the upper tail.
.pcbinom_tail <- function(c, N, k, p, lower.tail = TRUE) {
  m <- ceiling(c) - 1
  stats::pbinom(m, N, p, lower.tail = lower.tail) -
    k * p * (1 - p) * (stats::pbinom(m, N - 2, p, lower.tail = lower.tail) -
                         2 * stats::pbinom(m - 1, N - 2, p, lower.tail = lower.tail) +
                         stats::pbinom(m - 2, N - 2, p, lower.tail = lower.tail))
}

# Make a vector of (possibly slightly negative) probability weights usable by
# sample(): shift negative weights up as documented, then normalize.
.fix_weights <- function(weights) {
  if (any(weights < 0)) {
    warning("Density function returned negative probabilities.")
    weights <- weights + abs(min(weights))
  }
  weights / sum(weights)
}

# Normalize the `true.model` argument to "4P" or "2P".
.true_model <- function(true.model) {
  model <- substr(toupper(as.character(true.model)[1]), 1, 1)
  if (identical(model, "4")) return("4P")
  if (identical(model, "2")) return("2P")
  stop("`true.model` must be either \"4P\" or \"2P\".")
}

# Whether a set of fitted Beta parameters is impermissible (out of the [0, 1]
# interval, non-positive shape parameters, or not finite).
.impermissible <- function(params) {
  vals <- c(params$l, params$u, params$alpha, params$beta)
  any(!is.finite(vals)) || params$l < 0 || params$u > 1 || params$alpha <= 0 || params$beta <= 0
}

# Fit the true-score Beta distribution from its first four raw moments `m`
# (on the proportion scale). Returns the parameter list and a flag telling
# whether the two-parameter fail-safe was engaged.
.fit_true_score <- function(m, true.model, failsafe, l, u) {
  s2 <- m[2] - m[1]^2
  g3 <- (m[3] - 3 * m[1] * m[2] + 2 * m[1]^3) / sqrt(s2)^3
  g4 <- (m[4] - 4 * m[1] * m[3] + 6 * m[1]^2 * m[2] - 3 * m[1]^4) / sqrt(s2)^4
  engaged <- FALSE
  if (.true_model(true.model) == "4P") {
    params <- Beta.4p.fit(mean = m[1], variance = s2, skewness = g3, kurtosis = g4)
    if (isTRUE(failsafe) && .impermissible(params)) {
      warning("Fail-safe engaged: l = ", params$l, ", u = ", params$u,
              ", alpha = ", params$alpha, ", beta = ", params$beta,
              ". \n  Finding permissible solution for the true-score distribution in accordance with specifications.",
              call. = FALSE)
      params <- Beta.2p.fit(mean = m[1], variance = s2, l = l, u = u)
      engaged <- TRUE
    }
  } else {
    params <- Beta.2p.fit(mean = m[1], variance = s2, l = l, u = u)
  }
  list(params = params, engaged = engaged)
}

# Warn if fitted true-score parameters are outside the [0, 1] interval.
.warn_out_of_bounds <- function(params) {
  if (isTRUE(params$l < 0) || isTRUE(params$u > 1)) {
    warning("Parameter out of bounds: l = ", round(params$l, 4), ", u = ", round(params$u, 4),
            ", alpha = ", round(params$alpha, 4), ", beta = ", round(params$beta, 4),
            ". Consider constraining the fitting procedure further (e.g., set the location-parameters).",
            call. = FALSE)
  }
}

# Warn if observed values fall outside the declared [min, max] scale.
.warn_outside_scale <- function(x, min, max) {
  if (base::min(x) < min || base::max(x) > max) {
    warning("Observed values not within the specified [", min, ", ", max,
            "] bounds (observed min = ", base::min(x), ", observed max = ", base::max(x), ").",
            call. = FALSE)
  }
}

# Chi-square model-fit test from per-bin expected and observed counts.
#
# Bins whose expected count is below `minbin` are merged forward into the next
# bin; if the last bin is still too small it is merged into the one before it.
# Degrees of freedom follow Lord (1965, p. 244): number of bins minus one
# (for the fixed total) minus the number of moments matched when fitting the
# true-score distribution (2 for a two-parameter solution, 4 otherwise). No
# degree of freedom is deducted for the effective test length or Lord's k,
# which are estimated from the reliability coefficient rather than from the
# score distribution.
.chisq_modelfit <- function(expected, observed, minbin, n_params, rownames = NULL) {
  tab <- rbind(expected, observed)
  dimnames(tab) <- list(rownames, NULL)
  keep <- rep(TRUE, ncol(tab))
  for (i in seq_len(ncol(tab) - 1)) {
    if (tab[1, i] < minbin) {
      tab[, i + 1] <- tab[, i + 1] + tab[, i]
      keep[i] <- FALSE
    }
  }
  tab <- tab[, keep, drop = FALSE]
  J <- ncol(tab)
  if (J > 1 && tab[1, J] < minbin) {
    tab[, J - 1] <- tab[, J - 1] + tab[, J]
    tab <- tab[, -J, drop = FALSE]
  }
  chisquared <- sum((tab[2, ] - tab[1, ])^2 / tab[1, ])
  df <- ncol(tab) - 1 - n_params
  list(contingencytable = tab,
       chisquared = chisquared,
       df = df,
       pvalue = stats::pchisq(chisquared, df, lower.tail = FALSE))
}

# Row and column labels for the multiple-category accuracy and consistency
# matrices. `cuts` are the interior cut-points on the reporting scale.
.category_labels <- function(cuts, first_prefix) {
  K <- length(cuts) + 1
  lower <- c(NA, cuts)
  upper <- c(cuts, NA)
  vapply(seq_len(K), function(j) {
    if (j == 1) paste(first_prefix, "<", upper[j])
    else if (j < K) paste(" >=", lower[j], "& <", upper[j])
    else paste(" >=", lower[j])
  }, character(1))
}

# Overall and category-specific accuracy output for multiple-category
# classification, given a (normalized) confusion matrix with observed
# categories in rows and true categories in columns.
.mc_accuracy <- function(camat) {
  specific <- list()
  for (i in seq_len(ncol(camat))) {
    TP <- camat[i, i]
    FP <- sum(camat[i, -i])
    FN <- sum(camat[-i, i])
    TN <- sum(camat[-i, -i])
    specific[[paste0("Category.", i)]] <- list(confusionmatrix = confmat(TP, TN, FP, FN),
                                               statistics = caStats(TP, TN, FP, FN))
  }
  list(overall = list(confusionmatrix = camat, accuracy = sum(diag(camat))),
       specific = specific)
}

# Overall and category-specific consistency output for multiple-category
# classification, given a (normalized) consistency matrix.
.mc_consistency <- function(ccmat) {
  kappa <- function(p, p_c) (p - p_c) / (1 - p_c)
  p <- sum(diag(ccmat))
  p_c <- sum(colSums(ccmat)^2)
  specific <- list()
  for (i in seq_len(ncol(ccmat))) {
    p_i <- ccmat[i, i]
    p_c_i <- sum(ccmat[i, ])^2
    specific[[paste0("Category.", i)]] <- list(statistics = list(p = p_i, p_c = p_c_i, Kappa = kappa(p_i, p_c_i)))
  }
  list(overall = list(consistencymatrix = ccmat,
                      statistics = list(p = p, p_c = p_c, Kappa = kappa(p, p_c))),
       specific = specific)
}

# Two-category consistency output (matrix and statistics) from the four cell
# probabilities.
.binary_consistency <- function(ii, ij, ji, jj) {
  ccmat <- matrix(c(ii, ji, ij, jj), nrow = 2, dimnames = list(c("i", "j"), c("i", "j")))
  ccmat <- ccmat / sum(ccmat)
  list(consistencymatrix = ccmat,
       classification.consistency = ccStats(ccmat["i", "i"], ccmat["i", "j"],
                                            ccmat["j", "i"], ccmat["j", "j"]))
}

# ROC coordinates: evaluate `accuracy_at(cut)` (which returns the list produced
# by caStats()) at each cut-point.
.roc_matrix <- function(cuts, accuracy_at) {
  out <- matrix(nrow = length(cuts), ncol = 7,
                dimnames = list(NULL, c("FPR", "TPR", "Youden.J", "Cutoff", "Accuracy", "PPV", "NPV")))
  out[, "Cutoff"] <- cuts
  for (i in seq_along(cuts)) {
    s <- accuracy_at(cuts[i])
    out[i, c("FPR", "TPR", "Youden.J", "Accuracy", "PPV", "NPV")] <-
      c(1 - s$Specificity, s$Sensitivity, s$Youden.J, s$Accuracy, s$PPV, s$NPV)
  }
  # Undefined rates at the extremes of the scale (division by zero).
  out[is.na(out[, "FPR"]), "FPR"] <- 0
  out[is.na(out[, "TPR"]), "TPR"] <- 1
  out[is.na(out[, "PPV"]), "PPV"] <- 1
  out[is.na(out[, "NPV"]), "NPV"] <- 1
  out
}

# Draw an ROC curve from the output of .roc_matrix().
.roc_plot <- function(rocmat, truecut, AUC = FALSE, maxJ = FALSE, maxAcc = FALSE, locate = NULL) {
  graphics::plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
                 xlab = "False-Positive Rate (1 - Specificity)",
                 ylab = "True-Positive Rate (Sensitivity)",
                 main = paste("ROC curve for true-cut equal to", truecut))
  graphics::abline(h = seq(0, 1, .1), v = seq(0, 1, .1), col = "lightgrey", lty = "dotted")
  graphics::lines(rocmat[, "FPR"], rocmat[, "TPR"], type = "s", lwd = 2)

  mark <- function(rows, label) {
    graphics::points(rocmat[rows, "FPR"], rocmat[rows, "TPR"], cex = 1.5, pch = 19)
    graphics::text(rocmat[rows, "FPR"] + .025, rocmat[rows, "TPR"] - .025, labels = label, adj = c(0, 1))
  }
  mark_max <- function(column, what) {
    best <- max(rocmat[, column])
    rows <- which(rocmat[, column] == best)
    mark(rows, paste0("Maximum ", what, " (", round(best, 3), ") at cut-off = ",
                      round(rocmat[rows[1], "Cutoff"], 3), "."))
  }

  if (AUC) {
    graphics::legend("bottomright", bty = "n", cex = 1.5,
                     legend = paste("AUC =", round(AUC(rocmat[, "FPR"], rocmat[, "TPR"]), 3)))
  }
  if (maxJ) mark_max("Youden.J", "Youden's J.")
  if (maxAcc) mark_max("Accuracy", "Accuracy")
  if (!is.null(locate)) {
    target <- as.numeric(locate[2])
    what <- tolower(locate[1])
    if (startsWith(what, "se")) {
      row <- which(rocmat[, "TPR"] >= target)[1]
      value <- rocmat[row, "TPR"]
      statistic <- "Sensitivity >= "
    } else if (startsWith(what, "sp")) {
      rows <- which(rocmat[, "FPR"] <= 1 - target)
      row <- rows[length(rows)]
      value <- 1 - rocmat[row, "FPR"]
      statistic <- "Specificity >= "
    } else if (startsWith(what, "p")) {
      row <- which(rocmat[, "PPV"] <= target)[1]
      value <- rocmat[row, "PPV"]
      statistic <- "PPV <= "
    } else if (startsWith(what, "n")) {
      row <- which(rocmat[, "NPV"] >= target)[1]
      value <- rocmat[row, "NPV"]
      statistic <- "NPV >= "
    } else {
      stop("`locate[1]` must be one of \"sensitivity\", \"specificity\", \"PPV\", or \"NPV\".")
    }
    if (length(row) == 0 || is.na(row)) {
      warning("No cut-off satisfying `locate` was found.", call. = FALSE)
    } else {
      mark(row, paste0(statistic, target, " (", round(value, 3), ") at cut-off = ",
                       rocmat[row, "Cutoff"], "."))
    }
  }
  invisible(NULL)
}
