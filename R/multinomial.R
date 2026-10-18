#' An Implementation of the Compound Multinomial Approach to Estimate Classification Consistency and Accuracy based on Item Responses.
#'
#' @description An implementation of the multinomial and compound multinomial approach to classification consistency and accuracy of Lee (2007) and Lee, Brennan, and Wan (2009), as implemented in the MULT-CLASS program. Unlike the Livingston and Lewis and Hanson and Brennan approaches, no true-score distribution is fitted. Each examinee's items within a fixed item set are treated as a random sample from a universe of items, so that the examinee's observed proportions of responses at each score point estimate the probabilities of a multinomial error model. The distribution of an examinee's (weighted) total score on a replication is the convolution of these multinomial distributions across item sets ("compound" multinomial). Classification indices are averaged over examinees.
#' @param x A matrix or data frame of item scores, with rows representing examinees and columns representing items. Missing values are not allowed.
#' @param cut The cutoff value on the (weighted) total-score scale for classifying observations into above/below categories. Total scores at or above \code{cut} are classified as above the cut.
#' @param truecut Optional specification of a cutoff for the reference ("true") classification, on the same scale as \code{cut}. Defaults to \code{cut}. Useful for producing ROC curves (see documentation for the \code{MN.ROC()} function).
#' @param sets A list defining the fixed item sets, each element a vector of column indices or column names of \code{x}. Every column must belong to exactly one set, and each set must contain at least two items. Default is \code{NULL}, treating all items as a single set (the multinomial model).
#' @param weights A vector of positive weights, one per item set, applied to the set's summed score when forming the total score. Default is \code{NULL} (all weights equal to 1).
#' @param scorepoints Optional list (one element per item set) of the possible score points of the items in each set. Every observed item score must be one of these. Default is \code{NULL}, which uses the distinct item scores observed in each set. Score points that are never observed have probability zero for every examinee, so they affect only the reported range of possible total scores.
#' @param output Character vector indicating which types of statistics (i.e, accuracy and/or consistency) are to be computed and included in the output. Permissible values are \code{"accuracy"} and \code{"consistency"}.
#' @param correct Logical. Whether to apply the bias correction of Brennan and Lee (2006), in which each examinee's score-point proportions are shrunk toward the group mean proportions by the square root of the reliability of the total score. Default is \code{TRUE}.
#' @param conditional Logical. Whether to include the examinee-level (conditional) results in the output. Default is \code{FALSE}.
#' @return A list containing the model parameters (the model-based reliability of the total score, the observed total-score variance and average error variance, the number of items, and the range of possible total scores), the confusion matrix containing estimated proportions of true/false pass/fail categorizations for a test, diagnostic performance statistics, and / or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa). If \code{conditional = TRUE}, also a data frame of examinee-level results: the observed total score, the reference category, the probabilities of consistent (\code{phi}) and accurate (\code{gamma}) classification, the probabilities of classification above (\code{overclassification}) and below (\code{underclassification}) the reference category, the estimated error variance, and the probability of each category.
#' @note As in the LL and HB functions, the "positive" category in the confusion matrix is the category below the cut.
#' @note The reliability coefficient is model-based: one minus the ratio of the average estimated error variance to the observed total-score variance, where the error variance of each set's summed score is estimated within examinees with the \eqn{n / (n - 1)} correction, and the total-score variance uses the denominator \eqn{N} (the number of examinees), as in MULT-CLASS. It is computed from the data, so no reliability coefficient is supplied.
#' @note The reference ("true") category used for accuracy is the category of each examinee's observed total score with respect to \code{truecut}. This follows Lee, Brennan, and Wan (2009) and MULT-CLASS, and is not an independently known true classification. With \code{correct = TRUE}, only the category probabilities are shrunk; the reference category is unchanged.
#' @note Item sets should contain enough items for the within-examinee proportions to be informative; with few items per set the error model is coarse.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a 20-item test.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
#' rawdata <- sapply(1:20, function(i) rbinom(1000, 1, p.success))
#'
#' # Suppose the cutoff value for attaining a pass is 10 items correct. To
#' # estimate classification accuracy and consistency under the multinomial
#' # model (all items in one set) with bias correction:
#' MN.CA(rawdata, cut = 10)
#'
#' # Mixed-format test: 8 polytomous items scored 2/4/6/8 with weight 1, and
#' # 12 dichotomous items with weight 2.
#' poly <- matrix(sample(c(2, 4, 6, 8), 8000, replace = TRUE), ncol = 8)
#' mixed <- cbind(poly, rawdata[, 1:12])
#' MN.CA(mixed, cut = 50, sets = list(1:8, 9:20), weights = c(1, 2),
#' scorepoints = list(c(2, 4, 6, 8), c(0, 1)))
#' @references Lee, Won-Chan (2007). Multinomial and Compound Multinomial Error Models for Tests With Complex Item Scoring. Applied Psychological Measurement, 31(4), 255--274. doi: 10.1177/0146621606294206
#' @references Lee, Won-Chan, Brennan, Robert L., and Wan, Lei (2009). Classification Consistency and Accuracy for Complex Assessments Under the Compound Multinomial Model. Applied Psychological Measurement, 33(5), 374--390. doi: 10.1177/0146621608321759
#' @references Brennan, Robert L. and Lee, Won-Chan (2006). Correcting for Bias in Single-Administration Decision Consistency Indexes. CASMA Research Report 18. University of Iowa.
#' @export
MN.CA <- function(x, cut, truecut = NULL, sets = NULL, weights = NULL, scorepoints = NULL, output = c("accuracy", "consistency"), correct = TRUE, conditional = FALSE) {
  if (length(cut) != 1) {
    stop("`cut` must be a single value. For several cut-points, use MN.CA.MC().")
  }
  model <- .mn_model(x, sets, weights, scorepoints, correct)
  fit <- .mn_classify(model, cut, if (is.null(truecut)) cut else truecut)
  out <- list(parameters = model$parameters)
  output <- tolower(output)
  if (any(startsWith(output, "a"))) {
    a <- fit$accuracy
    # Category 1 (below the cut) is the "positive" category, as in LL.CA().
    camat <- confmat(a[1, 1], a[2, 2], a[1, 2], a[2, 1], "prop")
    out[["confusionmatrix"]] <- camat
    out[["classification.accuracy"]] <- caStats(camat[1, 1], camat[1, 2], camat[2, 1], camat[2, 2])
  }
  if (any(startsWith(output, "c"))) {
    b <- fit$consistency
    out <- c(out, .binary_consistency(b[1, 1], b[1, 2], b[2, 1], b[2, 2]))
  }
  if (conditional) {
    out[["conditional"]] <- fit$conditional
  }
  out
}

#' An Implementation of the Compound Multinomial Approach to Estimate Classification Consistency and Accuracy for Multiple Classifications based on Item Responses.
#'
#' @description An implementation of the multinomial and compound multinomial approach to classification consistency and accuracy of Lee (2007) and Lee, Brennan, and Wan (2009), as implemented in the MULT-CLASS program, for two or more categories. See the documentation for \code{MN.CA()} for a description of the model.
#' @param x A matrix or data frame of item scores, with rows representing examinees and columns representing items. Missing values are not allowed.
#' @param cut A vector of cut-off values on the (weighted) total-score scale for classifying observations into two or more categories. Total scores at or above a cut-off are classified into the category above it.
#' @param truecut Optional vector of cut-off values for the reference ("true") classification, of the same length as \code{cut}. Defaults to \code{cut}.
#' @param sets A list defining the fixed item sets, each element a vector of column indices or column names of \code{x}. Every column must belong to exactly one set, and each set must contain at least two items. Default is \code{NULL}, treating all items as a single set (the multinomial model).
#' @param weights A vector of positive weights, one per item set, applied to the set's summed score when forming the total score. Default is \code{NULL} (all weights equal to 1).
#' @param scorepoints Optional list (one element per item set) of the possible score points of the items in each set. Default is \code{NULL}, which uses the distinct item scores observed in each set.
#' @param correct Logical. Whether to apply the bias correction of Brennan and Lee (2006). Default is \code{TRUE}.
#' @param conditional Logical. Whether to include the examinee-level (conditional) results in the output. Default is \code{FALSE}.
#' @return A list containing the model parameters (see \code{MN.CA()}), the overall and category-specific confusion matrices and diagnostic performance statistics, and the overall and category-specific consistency matrices and statistics, in the same format as the output of \code{LL.CA.MC()} and \code{HB.CA.MC()}. The overall accuracy output additionally includes the ordered over- and underclassification rates of Lee, Brennan, and Wan (2009): the expected proportions of examinees classified above (\code{overclassification}) and below (\code{underclassification}) their reference category. MULT-CLASS reports these as "false positive" and "false negative" error. If \code{conditional = TRUE}, also a data frame of examinee-level results (see \code{MN.CA()}).
#' @note See the notes in the documentation for \code{MN.CA()} concerning the model-based reliability, the reference classification, and the bias correction.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a 20-item test.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
#' rawdata <- sapply(1:20, function(i) rbinom(1000, 1, p.success))
#'
#' # Estimate accuracy and consistency for three categories, with cut-points
#' # at 8 and 12 items correct:
#' (output <- MN.CA.MC(rawdata, cut = c(8, 12)))
#'
#' # The output can be organized in tabular format with MC.out.tabular():
#' MC.out.tabular(output)
#' @references Lee, Won-Chan (2007). Multinomial and Compound Multinomial Error Models for Tests With Complex Item Scoring. Applied Psychological Measurement, 31(4), 255--274. doi: 10.1177/0146621606294206
#' @references Lee, Won-Chan, Brennan, Robert L., and Wan, Lei (2009). Classification Consistency and Accuracy for Complex Assessments Under the Compound Multinomial Model. Applied Psychological Measurement, 33(5), 374--390. doi: 10.1177/0146621608321759
#' @references Brennan, Robert L. and Lee, Won-Chan (2006). Correcting for Bias in Single-Administration Decision Consistency Indexes. CASMA Research Report 18. University of Iowa.
#' @export
MN.CA.MC <- function(x, cut, truecut = NULL, sets = NULL, weights = NULL, scorepoints = NULL, correct = TRUE, conditional = FALSE) {
  if (is.null(truecut)) {
    truecut <- cut
  }
  if (length(truecut) != length(cut)) {
    stop("`cut` and `truecut` must have the same length.")
  }
  model <- .mn_model(x, sets, weights, scorepoints, correct)
  fit <- .mn_classify(model, cut, truecut)
  out <- list(parameters = model$parameters)

  camat <- fit$accuracy
  dimnames(camat) <- list(.category_labels(cut, "Observed"), .category_labels(truecut, "True    "))
  out[["accuracy"]] <- .mc_accuracy(camat)
  out$accuracy$overall[["overclassification"]] <- mean(fit$conditional$overclassification)
  out$accuracy$overall[["underclassification"]] <- mean(fit$conditional$underclassification)

  ccmat <- fit$consistency
  dimnames(ccmat) <- list(.category_labels(cut, "2nd adm."), .category_labels(cut, "1st adm."))
  out[["consistency"]] <- .mc_consistency(ccmat)
  if (conditional) {
    out[["conditional"]] <- fit$conditional
  }
  out
}

#' ROC curves for the compound multinomial approach.
#'
#' @description Generate a ROC curve plotting the false-positive rate against the true-positive rate at different cut-off values across the observed-score scale, under the multinomial or compound multinomial model (see documentation for the \code{MN.CA()} function).
#' @param x A matrix or data frame of item scores, with rows representing examinees and columns representing items.
#' @param truecut The point along the total-score scale that marks the reference classification.
#' @param sets A list defining the fixed item sets (see \code{MN.CA()}). Default is \code{NULL} (a single set).
#' @param weights A vector of positive weights, one per item set. Default is \code{NULL} (all equal to 1).
#' @param scorepoints Optional list of the possible score points of the items in each set (see \code{MN.CA()}).
#' @param correct Logical. Whether to apply the bias correction of Brennan and Lee (2006). Default is \code{TRUE}.
#' @param AUC Logical. Calculate and include the area under the curve? Default is \code{FALSE}.
#' @param maxJ Logical. Mark the point along the curve where Youden's J statistic is maximized? Default is \code{FALSE}.
#' @param maxAcc Logical. Mark the point along the curve where the Accuracy statistic is maximized? Default is \code{FALSE}.
#' @param locate Ask the function to locate the cut-point at which sensitivity or NPV is greater than or equal to some value, or specificity or PPV is lesser than or equal to some value. Take as input a character-vector of length 2, with the first argument being which index is to be found (e.g., "sensitivity"), and the second argument the value to locate (e.g., "0.75"). For example: c("sensitivity", "0.75").
#' @param raw.out Give raw coordinates as output rather than plot? Default is \code{FALSE}
#' @param grainsize Specify the number of cutoff-points for which the ROC curve is to be calculated, evenly spaced across the range of possible total scores. Default is 100.
#' @return A plot tracing the ROC curve for the test, or matrix of coordinates if raw.out is \code{TRUE}.
#' @note The examinee-level score distributions do not depend on the cut-off, so they are computed once and reused for every point on the curve.
#' @examples
#' # Generate some fictional data. Say, 1000 individuals take a 20-item test.
#' set.seed(1234)
#' p.success <- rBeta.4P(1000, 0.15, 0.85, 6, 4)
#' rawdata <- sapply(1:20, function(i) rbinom(1000, 1, p.success))
#'
#' # To draw the ROC curve for a true cut-off at 10 items correct and locate
#' # the points at which Youden's J and Accuracy are maximized:
#' MN.ROC(rawdata, truecut = 10, grainsize = 20, maxJ = TRUE, maxAcc = TRUE)
#' @references Lee, Won-Chan, Brennan, Robert L., and Wan, Lei (2009). Classification Consistency and Accuracy for Complex Assessments Under the Compound Multinomial Model. Applied Psychological Measurement, 33(5), 374--390. doi: 10.1177/0146621608321759
#' @export
MN.ROC <- function(x, truecut, sets = NULL, weights = NULL, scorepoints = NULL, correct = TRUE, AUC = FALSE, maxJ = FALSE, maxAcc = FALSE, locate = NULL, raw.out = FALSE, grainsize = 100) {
  model <- .mn_model(x, sets, weights, scorepoints, correct)
  range <- c(model$parameters$min, model$parameters$max)
  rocmat <- .roc_matrix(seq(range[1], range[2], diff(range) / grainsize), function(cut) {
    a <- .mn_classify(model, cut, truecut)$accuracy
    caStats(a[1, 1], a[2, 2], a[1, 2], a[2, 1])
  })
  if (raw.out) {
    return(rocmat)
  }
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  .roc_plot(rocmat, truecut, AUC = AUC, maxJ = maxJ, maxAcc = maxAcc, locate = locate)
}

# Validate the item-set design and compute everything about the compound
# multinomial model that does not depend on the cut-points: observed totals,
# error variances, reliability, (corrected) score-point probabilities, and each
# distinct examinee profile's distribution over the possible total scores.
.mn_model <- function(x, sets, weights, scorepoints, correct) {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("Item scores in `x` must be numeric.")
  }
  if (nrow(x) < 2) {
    stop("At least two examinees are required.")
  }
  if (anyNA(x)) {
    stop("Missing item scores are not allowed. Remove incomplete rows first.")
  }
  if (is.null(sets)) {
    sets <- list(seq_len(ncol(x)))
  }
  sets <- lapply(sets, function(s) if (is.character(s)) match(s, colnames(x)) else s)
  members <- unlist(sets)
  if (anyNA(members) || anyDuplicated(members) || !setequal(members, seq_len(ncol(x)))) {
    stop("`sets` must assign every column of `x` to exactly one item set.")
  }
  if (any(lengths(sets) < 2)) {
    stop("Each item set must contain at least two items (to estimate the within-set error variance).")
  }
  if (is.null(weights)) {
    weights <- rep(1, length(sets))
  }
  if (length(weights) != length(sets) || any(!is.finite(weights) | weights <= 0)) {
    stop("`weights` must contain one positive weight per item set.")
  }
  if (is.null(scorepoints)) {
    scorepoints <- lapply(sets, function(s) sort(unique(as.vector(x[, s]))))
  }
  if (!is.list(scorepoints) || length(scorepoints) != length(sets)) {
    stop("`scorepoints` must be a list with one vector of score points per item set.")
  }
  scorepoints <- lapply(scorepoints, function(p) sort(unique(p)))

  N <- nrow(x)
  total <- numeric(N)
  error <- numeric(N)
  props <- vector("list", length(sets))
  counts <- vector("list", length(sets))
  for (s in seq_along(sets)) {
    X <- x[, sets[[s]], drop = FALSE]
    points <- scorepoints[[s]]
    if (!all(X %in% points)) {
      stop("Item set ", s, " contains item scores that are not among its declared score points.")
    }
    n <- ncol(X)
    counts[[s]] <- vapply(points, function(p) rowSums(X == p), numeric(N))
    props[[s]] <- counts[[s]] / n
    # Error variance of the set's summed score: n times the within-examinee
    # item-score variance (denominator n - 1). Lee (2007).
    error <- error + weights[s]^2 * rowSums((X - rowMeans(X))^2) * n / (n - 1)
    total <- total + weights[s] * rowSums(X)
  }
  variance <- mean((total - mean(total))^2)  # denominator N, as in MULT-CLASS
  if (variance <= 0) {
    stop("The total score is constant across examinees; reliability and classification indices cannot be estimated.")
  }
  reliability <- 1 - mean(error) / variance

  if (correct) {
    if (!(reliability >= 0 && reliability <= 1)) {
      stop("The model-based reliability (", round(reliability, 4), ") is outside [0, 1], so the bias correction cannot be applied. Inspect the item-set design or set `correct = FALSE`.")
    }
    # Brennan and Lee (2006): shrink each examinee's proportions toward the
    # group mean by the square root of the total-score reliability, using the
    # same coefficient in every set.
    root <- sqrt(reliability)
    props <- lapply(props, function(p) {
      root * p + matrix((1 - root) * colMeans(p), nrow = nrow(p), ncol = ncol(p), byrow = TRUE)
    })
  } else if (!(reliability >= 0 && reliability <= 1)) {
    warning("The model-based reliability is outside [0, 1].", call. = FALSE)
  }

  # Examinees with identical score-point counts share a distribution.
  profile <- do.call(paste, c(as.data.frame(do.call(cbind, counts)), sep = ","))
  first <- !duplicated(profile)
  profile_id <- match(profile, profile[first])
  rows <- which(first)

  # Distribution of the total score for each distinct profile, over the full
  # set of attainable totals (identical for every profile, so a matrix).
  support <- 0
  pmf <- matrix(1, nrow = length(rows), ncol = 1)
  for (s in seq_along(sets)) {
    points <- weights[s] * scorepoints[[s]]
    for (item in seq_along(sets[[s]])) {
      step <- .mn_convolve(support, pmf, points, props[[s]][rows, , drop = FALSE])
      support <- step$support
      pmf <- step$pmf
    }
  }
  if (any(abs(rowSums(pmf) - 1) > 1e-8)) {
    stop("The score distribution failed the probability check.")
  }

  list(parameters = list(reliability = reliability,
                         observed.variance = variance,
                         error.variance = mean(error),
                         items = ncol(x),
                         sets = length(sets),
                         min = min(support),
                         max = max(support),
                         bias.correction = correct),
       total = total,
       error = error,
       support = support,
       pmf = pmf,
       profile = profile_id)
}

# Add one item with score points `points` and per-profile probabilities
# `probs` (profiles x points) to the distribution `pmf` (profiles x support).
.mn_convolve <- function(support, pmf, points, probs) {
  values <- as.vector(outer(support, points, "+"))
  key <- round(values, 9)
  new_support <- sort(unique(key))
  group <- match(key, new_support)
  new_pmf <- matrix(0, nrow = nrow(pmf), ncol = length(new_support))
  for (j in seq_along(points)) {
    idx <- group[(j - 1) * length(support) + seq_along(support)]
    new_pmf[, idx] <- new_pmf[, idx] + pmf * probs[, j]
  }
  list(support = new_support, pmf = new_pmf)
}

# Category of each value given increasing cut-points: values at or above a
# cut-point fall in the category above it.
.mn_category <- function(values, cuts) {
  tolerance <- 1e-9 * pmax(1, abs(cuts))
  1 + rowSums(outer(values, cuts - tolerance, ">="))
}

# Accuracy and consistency matrices and examinee-level results for given
# observed and reference cut-points.
.mn_classify <- function(model, cut, truecut) {
  if (any(diff(cut) <= 0) || any(diff(truecut) <= 0)) {
    stop("Cut-points must be strictly increasing.")
  }
  k <- length(cut) + 1
  # Probability of each observed category, per examinee.
  category <- .mn_category(model$support, cut)
  indicator <- outer(category, seq_len(k), "==") * 1
  q <- (model$pmf %*% indicator)[model$profile, , drop = FALSE]
  reference <- .mn_category(model$total, truecut)
  N <- nrow(q)

  accuracy <- matrix(0, k, k)  # rows observed, columns reference
  for (j in seq_len(k)) {
    accuracy[, j] <- colSums(q[reference == j, , drop = FALSE]) / N
  }
  consistency <- crossprod(q) / N
  above <- vapply(seq_len(N), function(i) sum(q[i, seq_len(k) > reference[i]]), numeric(1))
  below <- vapply(seq_len(N), function(i) sum(q[i, seq_len(k) < reference[i]]), numeric(1))
  conditional <- data.frame(total = model$total,
                            reference = reference,
                            phi = rowSums(q^2),
                            gamma = q[cbind(seq_len(N), reference)],
                            overclassification = above,
                            underclassification = below,
                            error.variance = model$error)
  colnames(q) <- paste0("P.Category.", seq_len(k))
  list(accuracy = accuracy, consistency = consistency, conditional = cbind(conditional, q))
}
